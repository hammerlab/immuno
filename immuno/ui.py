from vcf import load_vcf
from flask import Flask
from flask import render_template
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.user import (current_user, login_required, UserManager,
    UserMixin, SQLAlchemyAdapter)
import os
from os.path import exists, join, isdir

class ConfigClass(object):
    # Custom config
    VARIANT_PATH = os.environ.get("IMMUNO_VARIANT_PATH", os.getcwd())
    HLA_PATH = os.environ.get("IMMUNO_HLA_PATH", os.getcwd())
    # Flask config
    SECRET_KEY = os.environ.get("IMMUNO_SECRET_KEY")
    SQLALCHEMY_DATABASE_URI = os.environ.get("IMMUNO_DB")
    # Flask-User config
    USER_ENABLE_EMAIL = False

app = Flask(__name__)
app.config.from_object(__name__ + '.ConfigClass')
db = SQLAlchemy(app)

class User(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key = True)
    active = db.Column(db.Boolean(), nullable = False, default = False)
    username = db.Column(db.String(50), nullable = False, unique = True)
    password = db.Column(db.String(255), nullable = False, default = '')

db_adapter = SQLAlchemyAdapter(db, User)
user_manager = UserManager(db_adapter, app)

@app.route("/")
def patients():
    if not current_user.is_authenticated():
        return render_template("home.html")
    return render_template("patients.html", message = "test message")

@app.route("/patient/<patient_filename>")
def patient(patient_filename):
    if not current_user.is_authenticated():
        return render_template("home.html")
    variant_file_path = join(app.config["VARIANT_PATH"], patient_filename)
    if not exists(variant_file_path):
        return "File not found: %s" % variant_file_path
    if not variant_file_path.endswith(".vcf"):
        return "Not a VCF file: %s" % variant_file_path
    vcf_df = load_vcf(variant_file_path)
    vcf_rows = [row for _, row in vcf_df.iterrows()]
    return render_template("patient.html",
        patient_id = variant_file_path,
        variant_filename = variant_file_path,
        vcf = vcf_rows)

if __name__ == '__main__':
    app.debug = True
    assert isdir(app.config["VARIANT_PATH"]), \
        "Variant path %s must be a directory" % app.config["VARIANT_PATH"]
    assert isdir(app.config["HLA_PATH"]), \
        "HLA path %s must be a directory" % app.config["HLA_PATH"]
    app.run()
