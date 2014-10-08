from vcf import load_vcf
from flask import Flask
from flask import render_template
import os
from os.path import exists, join, isdir

app = Flask(__name__)

# Default to the current working directory if these variables are not set
VARIANT_PATH = os.environ.get("VARIANT_PATH", os.getcwd())
HLA_PATH = os.environ.get("HLA_PATH", os.getcwd())

@app.route('/')
def index():
    return render_template("home.html")

@app.route("/patients")
def patients():
    return render_template("patients.html", message = "test message")

@app.route("/patient/<patient_filename>")
def patient(patient_filename):
    variant_file_path = join(VARIANT_PATH, patient_filename)
    if not exists(variant_file_path):
        return "File not found: %s" % variant_file_path
    if not variant_file_path.endswith(".vcf"):
        return "Not a VCF file: %s" % variant_file_path
    vcf_df = load_vcf(variant_file_path)
    return render_template("patient.html",
        patient_id = variant_file_path,
        variant_filename = variant_file_path,
        df_html = vcf_df.to_html())

if __name__ == '__main__':
    app.debug = True
    assert isdir(VARIANT_PATH), \
        "Variant path %s must be a directory" % VARIANT_PATH
    assert isdir(HLA_PATH), \
        "HLA path %s must be a directory" % HLA_PATH
    app.run()
