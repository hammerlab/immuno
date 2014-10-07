import argparse
from os.path import exists, join, isdir
from vcf import load_vcf
from flask import Flask

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
    return 'Immuno UI'

@app.route('/patients')
def patients():
    return 'Patients'

@app.route("/patient/<patient_filename>")
def patient(patient_filename):

    variant_path = join(args.variant_path, patient_filename)
    if not exists(variant_path):
        return "File not found: %s" % variant_path
    if not variant_path.endswith(".vcf"):
        return "Not a VCF file: %s" % variant_path
    vcf_df = load_vcf(variant_path)
    return """
    <html>
    <head><title>Patient ID: %s</title></head>
    <body>
    <div>
    Variant filename: %s
    </div>
    <h2>
    Variants
    </h2>
    %s
    </body>
    </html>
    """ % (
        variant_path,
        variant_path,
        vcf_df.to_html()
    )


if __name__ == '__main__':
    args = parser.parse_args()
    assert isdir(args.variant_path), \
        "Variant path %s must be a directory" % args.variant_path
    assert isdir(args.hla_path), \
        "HLA path %s must be a directory" % args.hla_path
    app.debug = True
    app.run()
