from common import str2bool
from vcf import load_vcf

from flask import Flask
from flask import redirect, request, render_template, url_for,
    send_from_directory
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.user import (current_user, login_required, UserManager,
    UserMixin, SQLAlchemyAdapter)
from flask_mail import Mail, Message
from werkzeug import secure_filename
from os import environ, getcwd
from os.path import exists, join
from vcf import load_vcf

ALLOWED_EXTENSIONS = set(['vcf', 'hla'])

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
    mhc_class = db.Column(db.String(1), nullable=False)

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
    return render_template('profile.html')

@app.route('/patient/<patient_filename>')
@login_required
def patient(patient_filename):
    file_names = File.query.with_entities(File.name).filter_by(
        user_id=current_user.id).all()
    # TODO: Loop through all files
    filename = file_names[0][0]
    print join(app.config['UPLOAD_FOLDER'], filename)
    vcf_df = load_vcf(join(app.config['UPLOAD_FOLDER'], filename))
    vcf_rows = [row for _, row in vcf_df.iterrows()]
    return render_template('patient.html',
        patient_id = filename,
        variant_filename = filename,
        vcf = vcf_rows)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route('/upload', methods=['GET', 'POST'])
@login_required
def upload_file():
    """TODO: Implement streaming, or at least delete the file when done."""
    if request.method == 'POST':
        filenames = []
        files = request.files.getlist('file')
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                print filename
                filepath = join(app.config['UPLOAD_FOLDER'], filename)
                file.save(filepath)
                variants = create_variants(filepath)
                db.session.add_all(variants)
                db.session.commit()
                filenames.append(filename)
        return render_template('uploaded.html', filenames=filenames)
    return render_template('upload.html')

def create_variants(filepath):
    vcf_df = load_vcf(filepath)
    variants = []
    for index, row in vcf_df.iterrows():
        patient_id = "testpatient"
        chr = row['chr']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        variant = Variant(patient_id=patient_id, chr=chr, pos=pos,
            ref=ref, alt=alt)
        variants.append(variant)
    return variants

@app.route('/uploads/<filename>')
@login_required
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)
