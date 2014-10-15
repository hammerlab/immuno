from common import str2bool
from vcf import load_vcf

from flask import Flask
from flask import redirect, render_template, url_for
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.user import (current_user, login_required, UserManager,
    UserMixin, SQLAlchemyAdapter)
from flask_mail import Mail, Message
from os import environ, getcwd
from os.path import exists, join

class ConfigClass(object):
    # Custom config
    DEBUG = str2bool(environ.get('IMMUNO_DEBUG', 'False'))
    PORT = int(environ.get('IMMUNO_PORT', 5000))
    USE_RELOADER = str2bool(environ.get('IMMUNO_USE_RELOADER', False))
    VARIANT_PATH = environ.get('IMMUNO_VARIANT_PATH', getcwd())
    HLA_PATH = environ.get('IMMUNO_HLA_PATH', getcwd())

    # Flask config
    SECRET_KEY = environ.get('IMMUNO_SECRET_KEY')
    SQLALCHEMY_DATABASE_URI = environ.get('IMMUNO_DB')

    # Flask-User config
    USER_ENABLE_EMAIL = True
    USER_ENABLE_CHANGE_PASSWORD = True
    USER_ENABLE_CHANGE_USERNAME = False
    USER_ENABLE_CONFIRM_EMAIL = True
    USER_ENABLE_FORGOT_PASSWORD = True
    USER_ENABLE_MULTIPLE_EMAILS = False
    USER_ENABLE_REGISTRATION = True
    USER_ENABLE_RETYPE_PASSWORD = True
    USER_ENABLE_USERNAME = False
    USER_CONFIRM_EMAIL_EXPIRATION = 2*24*3600
    USER_PASSWORD_HASH = 'bcrypt'
    USER_PASSWORD_HASH_MODE = 'passlib'
    USER_REQUIRE_INVITATION = False
    USER_RESET_PASSWORD_EXPIRATION = 2*24*3600
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
    variant_file_path = join(app.config['VARIANT_PATH'], patient_filename)
    if not exists(variant_file_path):
        return 'File not found: %s' % variant_file_path
    if not variant_file_path.endswith('.vcf'):
        return 'Not a VCF file: %s' % variant_file_path
    vcf_df = load_vcf(variant_file_path)
    vcf_rows = [row for _, row in vcf_df.iterrows()]
    return render_template('patient.html',
        patient_id = variant_file_path,
        variant_filename = variant_file_path,
        vcf = vcf_rows)
