from common import str2bool
from vcf import load_vcf
from mutation_report import normalize_hla_allele_name, group_epitopes
from load_file import expand_transcripts
from immunogenicity import ImmunogenicityPredictor
import mhc_random

from flask import Flask
from flask import (redirect, request, render_template, url_for,
    send_from_directory)
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.user import (current_user, login_required, UserManager,
    UserMixin, SQLAlchemyAdapter)
from flask_mail import Mail, Message
from flask.ext.wtf import Form
from flask.ext.wtf.file import FileField, FileRequired, FileAllowed
from werkzeug import secure_filename
from os import environ, getcwd
from os.path import exists, join
from vcf import load_vcf
from hla_file import read_hla_file
from wtforms import SubmitField, TextField, TextAreaField, validators
from jinja2 import ChoiceLoader, FileSystemLoader
from pandas import DataFrame, Series, concat

class ConfigClass(object):
    # Custom config
    DEBUG = str2bool(environ.get('IMMUNO_DEBUG', 'False'))
    PORT = int(environ.get('IMMUNO_PORT', 5000))
    USE_RELOADER = str2bool(environ.get('IMMUNO_USE_RELOADER', False))
    UPLOAD_FOLDER = join(getcwd(), 'uploads')

    # Flask config
    SECRET_KEY = environ.get('IMMUNO_SECRET_KEY')
    SQLALCHEMY_DATABASE_URI = environ.get('IMMUNO_DB')

    # Flask-User config
    USER_PRODUCT_NAME = 'immuno'
    USER_ENABLE_EMAIL = True
    USER_ENABLE_CHANGE_PASSWORD = True
    USER_ENABLE_CHANGE_USERNAME = False
    USER_ENABLE_CONFIRM_EMAIL = True
    USER_ENABLE_FORGOT_PASSWORD = True
    USER_ENABLE_MULTIPLE_EMAILS = False
    USER_ENABLE_REGISTRATION = True
    USER_ENABLE_RETYPE_PASSWORD = True
    USER_ENABLE_USERNAME = False
    USER_CONFIRM_EMAIL_EXPIRATION = 2 * 24 * 3600
    USER_PASSWORD_HASH = 'bcrypt'
    USER_PASSWORD_HASH_MODE = 'passlib'
    USER_REQUIRE_INVITATION = False
    USER_RESET_PASSWORD_EXPIRATION = 2 * 24 * 3600
    USER_SEND_PASSWORD_CHANGED_EMAIL = True
    USER_SEND_REGISTERED_EMAIL = True
    USER_SEND_USERNAME_CHANGED_EMAIL = False

    # Flask-Mail config
    MAIL_SERVER = environ.get('IMMUNO_MAIL_SERVER')
    MAIL_PORT = int(environ.get('IMMUNO_MAIL_PORT', 5000))
    MAIL_USE_SSL = str2bool(environ.get('IMMUNO_MAIL_USE_SSL', 'False'))
    MAIL_USE_TLS = str2bool(environ.get('IMMUNO_MAIL_USE_TLS', 'False'))
    MAIL_USERNAME = environ.get('IMMUNO_MAIL_USERNAME')
    MAIL_PASSWORD = environ.get('IMMUNO_MAIL_PASSWORD')
    MAIL_DEFAULT_SENDER = environ.get('IMMUNO_MAIL_DEFAULT_SENDER')

app = Flask(__name__)
app.config.from_object(__name__ + '.ConfigClass')
mail = Mail()
mail.init_app(app)
db = SQLAlchemy(app)

class User(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key=True)
    active = db.Column(db.Boolean(), nullable=False, default=False)
    password = db.Column(db.String(255), nullable=False, default='')
    email = db.Column(db.String(255), nullable=False, unique=True)
    confirmed_at = db.Column(db.DateTime())
    reset_password_token = db.Column(db.String(100), nullable=False, default='')

class Patient(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    name = db.Column(db.String(255), nullable=False, unique=True)
    notes = db.Column(db.Text, nullable=True)

class Variant(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    patient_id = db.Column(db.Integer, db.ForeignKey('patient.id'))
    chr = db.Column(db.String(10), nullable=False)
    pos = db.Column(db.Integer, nullable=False)
    ref = db.Column(db.String(1000), nullable=True)
    alt = db.Column(db.String(1000), nullable=False)

    def __init__(self, patient_id, chr, pos, ref, alt):
        self.patient_id = patient_id
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt

class HLAType(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    patient_id = db.Column(db.Integer, db.ForeignKey('patient.id'))
    allele = db.Column(db.String(15), nullable=False)
    mhc_class = db.Column(db.SmallInteger, nullable=False)

    def __init__(self, patient_id, allele, mhc_class):
        self.patient_id = patient_id
        self.allele = allele
        self.mhc_class = mhc_class

db_adapter = SQLAlchemyAdapter(db, User)
user_manager = UserManager(db_adapter, app)

@app.route('/')
def patients():
    if current_user.is_authenticated():
        return redirect(url_for('profile'))
    else:
        return redirect(url_for('user.login'))

@app.route('/profile')
@login_required
def profile():
    patients = Patient.query.with_entities(Patient.name).filter_by(
        user_id=current_user.id).all()
    return render_template('profile.html', patients=patients)

def get_vcf_df(patient_id):
    variants = Variant.query.with_entities(Variant.chr, Variant.pos,
        Variant.ref, Variant.alt).filter_by(patient_id=patient_id).all()
    vcf_df = DataFrame(variants, columns=['chr', 'pos', 'ref', 'alt'])

    # TODO: I added this because downstream functions expect 'info', but this
    # is silly and hacky.
    vcf_df['info'] = Series([None] * len(vcf_df))
    return vcf_df

@app.route('/patient/<name>')
@login_required
def patient(name):
    """This is version of the patient page that literally runs the mutation report
    pipeline synchronously. This is not good.
    """
    id, name, notes = Patient.query.with_entities(Patient.id, Patient.name,
        Patient.notes).filter_by(name=name).one()
    variants = Variant.query.with_entities(Variant.chr, Variant.pos,
        Variant.ref, Variant.alt).filter_by(patient_id=id).all()
    hla_types = HLAType.query.with_entities(HLAType.allele,
        HLAType.mhc_class).filter_by(patient_id=id).all()

    peptide_length = 9
    alleles = [normalize_hla_allele_name(
        allele) for allele, mhc_class in hla_types]

    vcf_df = get_vcf_df(id)
    transcripts_df, vcf_df, variant_report = expand_transcripts(
        vcf_df,
        id,
        min_peptide_length = peptide_length,
        max_peptide_length = peptide_length)

    # TODO: Don't use MHC random
    scored_epitopes = mhc_random.generate_scored_epitopes(transcripts_df, alleles)
    imm = ImmunogenicityPredictor(alleles = alleles)
    scored_epitopes = imm.predict(scored_epitopes)
    peptides = group_epitopes(scored_epitopes)

    return render_template('patient.html',
        name=name,
        notes=notes,
        variants=variants,
        hla_types=hla_types,
        peptides=peptides)

class NewPatientForm(Form):
    name = TextField('Patient Name',
        validators=[validators.required(), validators.length(max=10)])
    notes = TextAreaField('Patient Notes',
        validators=[validators.optional(), validators.length(max=200)])
    vcf_file = FileField('VCF File',
        validators=[FileRequired(), FileAllowed(['vcf'], 'VCF Only')])
    hla_file = FileField('HLA File',
        validators=[FileRequired(), FileAllowed(['hla'], 'HLA Only')])
    submit = SubmitField("Send")

@app.route('/patient/new', methods=['GET', 'POST'])
@login_required
def new_patient():
    """TODO: Implement streaming, or at least delete the file when done."""
    form = NewPatientForm(csrf_enabled=False)
    if form.validate_on_submit():
        patient = create_patient(request=request, user_id=current_user.id)
        db.session.add(patient)
        db.session.commit()

        variants = create_variants(file=request.files['vcf_file'],
            patient_id=patient.id)
        db.session.add_all(variants)
        db.session.commit()

        hla_types = create_hla_types(file=request.files['hla_file'],
            patient_id=patient.id)
        db.session.add_all(hla_types)
        db.session.commit()
        return redirect(url_for('patient', name=patient.name))
    return render_template('upload.html', form=form)

def create_patient(request, user_id):
    name = request.form['name']
    notes = request.form['notes']
    patient = Patient(user_id=user_id, name=name, notes=notes)
    return patient

def create_variants(file, patient_id):
    filename = secure_filename(file.filename)
    filepath = join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)
    vcf_df = load_vcf(filepath)
    variants = []
    for index, row in vcf_df.iterrows():
        chr = row['chr']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        variant = Variant(patient_id=patient_id, chr=chr, pos=pos,
            ref=ref, alt=alt)
        variants.append(variant)
    return variants

def create_hla_types(file, patient_id):
    filename = secure_filename(file.filename)
    filepath = join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)
    alleles = read_hla_file(filepath)
    hla_types = []
    for allele in alleles:
        hla_type = HLAType(patient_id=patient_id, allele=allele, mhc_class=1)
        hla_types.append(hla_type)
    return hla_types

@app.route('/uploads/<filename>')
@login_required
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)
