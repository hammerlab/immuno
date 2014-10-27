from os import environ, getcwd
from os.path import exists, join

from common import str2bool, env_var
from hla_file import read_hla_file
from immunogenicity import ImmunogenicityPredictor
from mhc_common import (normalize_hla_allele_name,
    mhc_class_from_normalized_allele_name)
from mhc_netmhcpan import PanBindingPredictor
from mhc_netmhccons import ConsensusBindingPredictor
import mhc_random
from mutation_report import group_epitopes
from load_file import expand_transcripts
from load_file import load_variants
from vcf import load_vcf

from flask import Flask
from flask import (redirect, request, render_template, url_for,
    send_from_directory, flash, make_response)
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.user import (current_user, login_required, UserManager,
    UserMixin, SQLAlchemyAdapter)
from flask_mail import Mail, Message
from flask.ext.wtf import Form
from flask.ext.wtf.file import FileField, FileRequired, FileAllowed
from jinja2 import ChoiceLoader, FileSystemLoader
from natsort import natsorted
from pandas import DataFrame, Series, concat, merge, read_json
from werkzeug import secure_filename
from wtforms import (SubmitField, TextField, TextAreaField, SelectField,
    validators)

class ConfigClass(object):
    # Custom config
    DEBUG = env_var('IMMUNO_DEBUG', str2bool, False)
    PORT = env_var('IMMUNO_PORT', int, 5000)
    USE_RELOADER = env_var('IMMUNO_USE_RELOADER', str2bool, False)
    UPLOAD_FOLDER = join(getcwd(), 'uploads')

    # Flask config
    SECRET_KEY = environ.get('IMMUNO_SECRET_KEY')
    assert SECRET_KEY, \
        "Environment variable IMMUNO_SECRET_KEY must be set"
    SQLALCHEMY_DATABASE_URI = environ.get('IMMUNO_DB')
    assert SQLALCHEMY_DATABASE_URI, \
        "Environment variable IMMUNO_DB must be set"
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
    assert MAIL_SERVER, \
        "Environment variable IMMUNO_MAIL_SERVER must be set"
    MAIL_PORT = env_var('IMMUNO_MAIL_PORT', int, 5000)
    MAIL_USE_SSL = env_var('IMMUNO_MAIL_USE_SSL', str2bool, False)
    MAIL_USE_TLS = env_var('IMMUNO_MAIL_USE_TLS', str2bool, False)
    MAIL_USERNAME = environ.get('IMMUNO_MAIL_USERNAME')
    assert MAIL_USERNAME, \
        "Environment variable IMMUNO_MAIL_USERNAME must be set"
    MAIL_PASSWORD = environ.get('IMMUNO_MAIL_PASSWORD')
    assert MAIL_PASSWORD, \
        "Environment variable IMMUNO_MAIL_PASSWORD must be set"
    MAIL_DEFAULT_SENDER = environ.get('IMMUNO_MAIL_DEFAULT_SENDER')
    assert MAIL_DEFAULT_SENDER, \
        "Environment variable IMMUNO_MAIL_DEFAULT_SENDER must be set"

app = Flask(__name__)
app.config.from_object(__name__ + '.ConfigClass')
mail = Mail()
mail.init_app(app)
db = SQLAlchemy(app)
print 'DB: %s' % db

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
    display_id = db.Column(db.String(1000), nullable=False, unique=True)

class Variant(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    patient_id = db.Column(db.Integer, db.ForeignKey('patient.id'))
    chr = db.Column(db.String(255), nullable=False)
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

class Run(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    patient_id = db.Column(db.Integer, db.ForeignKey('patient.id'))
    output = db.Column(db.Text, nullable=False)

    def __init__(self, patient_id, output):
        self.patient_id = patient_id
        self.output = output

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
    patients = Patient.query.with_entities(Patient.display_id).filter_by(
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

def run_pipeline(patient_id, score_epitopes):
    """Run the pipeline for this patient, and save the output to the DB as a
    Run."""
    hla_types = HLAType.query.with_entities(HLAType.allele,
        HLAType.mhc_class).filter_by(patient_id=patient_id).all()

    peptide_length = 31
    alleles = [normalize_hla_allele_name(
        allele) for allele, mhc_class in hla_types]

    vcf_df = get_vcf_df(patient_id)
    transcripts_df, vcf_df, variant_report = expand_transcripts(
        vcf_df,
        patient_id,
        min_peptide_length = peptide_length,
        max_peptide_length = peptide_length)

    scored_epitopes = score_epitopes(transcripts_df, alleles)
    imm = ImmunogenicityPredictor(alleles=alleles)
    scored_epitopes = imm.predict(scored_epitopes)

    # TODO(tavi) Make this expansion more robust. It breaks the IEDB predictor,
    # for example.
    short_transcripts_df = transcripts_df[['chr', 'pos', 'ref',
        'alt', 'TranscriptId']]
    scored_epitopes = merge(scored_epitopes, short_transcripts_df,
        on='TranscriptId', how='left')

    run = Run(patient_id=patient_id, output=scored_epitopes.to_json())
    db.session.add(run)

@app.route('/patient/<display_id>')
@login_required
def patient(display_id):
    patient_id, display_id = Patient.query.with_entities(Patient.id,
        Patient.display_id).filter_by(display_id=display_id).one()
    output = Run.query.with_entities(Run.output).filter_by(
        patient_id=patient_id).one()
    scored_epitopes = read_json(output[0])

    peptides = group_epitopes(scored_epitopes, use_transcript_name=True)
    peptides_no_thymic_deletion = group_epitopes(scored_epitopes,
        use_transcript_name=True, no_thymic_deletion=True)

    return render_template('patient.html',
        display_id=display_id,
        peptides=peptides,
        peptides_no_thymic_deletion=peptides_no_thymic_deletion)

@app.route('/patient/hla_types/<display_id>')
@login_required
def hla_types(display_id):
    patient_id, display_id = Patient.query.with_entities(Patient.id,
        Patient.display_id).filter_by(display_id=display_id).one()
    hla_types = HLAType.query.with_entities(HLAType.allele,
        HLAType.mhc_class).filter_by(patient_id=patient_id).all()

    return render_template('hla.html',
        display_id=display_id,
        hla_types=hla_types)

@app.route('/patient/variants/<display_id>')
@login_required
def variants(display_id):
    patient_id, display_id = Patient.query.with_entities(Patient.id,
        Patient.display_id).filter_by(display_id=display_id).one()
    variants = Variant.query.with_entities(Variant.chr, Variant.pos,
        Variant.ref, Variant.alt).filter_by(patient_id=patient_id).all()

    # Sort variants with a natural alphanumeric sort on chr name + position
    # TODO(tavi) Make this two-column sort less hacky
    variants = natsorted(variants,
        key=lambda variant: variant[0].lower() + ' ' + str(variant[1]))
    return render_template('variants.html',
        display_id=display_id,
        variants=variants)

@app.route('/patient/export/<display_id>')
@login_required
def export(display_id):
    patient_id, display_id = Patient.query.with_entities(Patient.id,
        Patient.display_id).filter_by(display_id=display_id).one()
    output = Run.query.with_entities(Run.output).filter_by(
        patient_id=patient_id).one()
    response = make_response(output[0])
    response.headers['Content-Disposition'] = ('attachment;' +
        'filename=epitopes.json')
    return response

class NewPatientForm(Form):
    binding_predictor_choices = {'random': 'Random',
        'netmhccons': 'NetMHCCons',
        'netmhcpan': 'NetMHCPan'}

    display_id = TextField('Patient ID',
        validators=[validators.required(), validators.length(max=1000)])
    vcf_file = FileField('VCF/MAF File',
        validators=[
            FileRequired(),
            FileAllowed(['vcf', 'maf'],
            'Variant file must be VCF or MAF')]
    )
    hla_file = FileField('HLA File',
        validators=[FileRequired(), FileAllowed(['hla'], 'HLA Only')])
    rna_file = FileField('RNA Expression Quantification File',
        validators=[FileAllowed(['txt'], 'TXT Only')])
    binding_predictor = SelectField('MHC Binding Prediction',
        choices=sorted(zip(binding_predictor_choices.keys(),
            binding_predictor_choices.values())),
        validators=[validators.required()])
    submit = SubmitField('Send')

def flash_errors(form):
    for field, errors in form.errors.items():
        for error in errors:
            print error
            flash('Error in the %s field. %s' % (
                getattr(form, field).label.text,
                error
            ), 'error')

def score_epitopes_random(transcripts_df, alleles):
    return mhc_random.generate_scored_epitopes(transcripts_df, alleles)

def score_epitopes_mhc_cons(transcripts_df, alleles):
    predictor = ConsensusBindingPredictor(alleles)
    return predictor.predict(transcripts_df)

def score_epitopes_mhc_pan(transcripts_df, alleles):
    predictor = PanBindingPredictor(alleles)
    return predictor.predict(transcripts_df)

@app.route('/patient/new', methods=['GET', 'POST'])
@login_required
def new_patient():
    """Upload files and run the pipeline for a new patient.
    TODO: Don't run the pipeline synchronously.
    TODO: Implement streaming, or at least delete the file when done.
    """
    form = NewPatientForm(csrf_enabled=False)
    if form.validate_on_submit():
        patient = create_patient(request=request, user_id=current_user.id)
        db.session.add(patient)
        db.session.flush()

        variants = create_variants(file=request.files['vcf_file'],
            patient_id=patient.id)
        db.session.add_all(variants)

        hla_types = create_hla_types(file=request.files['hla_file'],
            patient_id=patient.id)
        db.session.add_all(hla_types)

        binding_predictor = request.form['binding_predictor']
        assert binding_predictor in form.binding_predictor_choices, \
            ('Binding predictor %s should have produced a form error') % (
                binding_predictor)

        if binding_predictor == 'random':
            score_epitopes = score_epitopes_random
        elif binding_predictor == 'netmhccons':
            score_epitopes = score_epitopes_mhc_cons
        else:
            score_epitopes = score_epitopes_mhc_pan

        run_pipeline(patient.id, score_epitopes)
        db.session.commit()
        return redirect(url_for('patient', display_id=patient.display_id))
    flash_errors(form)
    return render_template('upload.html', form=form)

def create_patient(request, user_id):
    display_id = request.form['display_id']
    patient = Patient(user_id=user_id, display_id=display_id)
    return patient

def create_variants(file, patient_id):
    filename = secure_filename(file.filename)
    filepath = join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)
    vcf_df = load_variants(filepath)
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
        allele_normalized = normalize_hla_allele_name(allele)
        mhc_class = mhc_class_from_normalized_allele_name(allele_normalized)
        hla_type = HLAType(patient_id=patient_id, allele=allele_normalized,
            mhc_class=mhc_class)
        hla_types.append(hla_type)
    return hla_types
