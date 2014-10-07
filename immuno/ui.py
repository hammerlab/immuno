import argparse
from os.path import exists, join, isdir
from vcf import load_vcf
from flask import Flask
from flask import render_template

parser = argparse.ArgumentParser(
    description="Web UI for inspecting cancer neoantigens")

parser.add_argument(
    "--variant-path",
    default = ".",
    help="Path to directory with VCF/MAF files with patient mutations")

parser.add_argument(
    "--hla-path",
    default = ".",
    help="Path directory with files containing patient HLA types")

app = Flask(__name__)

@app.route('/')
def index():
    return render_template("home.html")

@app.route("/patients")
def patients():
    return render_template("patients.html", message = "test message")

@app.route("/patient/<patient_filename>")
def patient(patient_filename):
    variant_path = join(args.variant_path, patient_filename)
    if not exists(variant_path):
        return "File not found: %s" % variant_path
    if not variant_path.endswith(".vcf"):
        return "Not a VCF file: %s" % variant_path
    vcf_df = load_vcf(variant_path)
    return render_template("patient.html",
        patient_id = variant_path,
        variant_filename = variant_path,
        df_html = vcf_df.to_html())

if __name__ == '__main__':
    args = parser.parse_args()
    assert isdir(args.variant_path), \
        "Variant path %s must be a directory" % args.variant_path
    assert isdir(args.hla_path), \
        "HLA path %s must be a directory" % args.hla_path
    app.debug = True
    app.run()
