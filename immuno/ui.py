from vcf import load_vcf
from flask import Flask
from flask import render_template
from os.path import exists, join, isdir

app = Flask(__name__)

@app.route('/')
def index():
    return render_template("home.html")

@app.route("/patients")
def patients():
    return render_template("patients.html", message = "test message")

@app.route("/patient/<patient_filename>")
def patient(patient_filename):
    variant_file_path = join(app.config["VARIANT_PATH"], patient_filename)
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
    app.config.from_object('config')
    assert isdir(app.config["VARIANT_PATH"]), (
        "Variant path %s must be a directory" % app.config["VARIANT_PATH"])
    assert isdir(app.config["HLA_PATH"]), (
        "HLA path %s must be a directory" % app.config["HLA_PATH"])
    app.run()
