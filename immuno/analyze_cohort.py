import argparse 
import logging
from os import listdir
from os.path import join, split 
from glob import glob
from collections import OrderedDict

from common import init_logging
from immunogenicity import ImmunogenicityPredictor
from load_file import load_file
from mhc_common import normalize_hla_allele_name
from mhc_netmhcpan import PanBindingPredictor
from mutation_report import print_mutation_report

parser = argparse.ArgumentParser()


parser.add_argument("--input-dir",
  	type = str, 
  	required = True, 
    help="Directory containing MAF or VCF input files")

parser.add_argument("--hla-input-dir",
	type = str, 
	default = None, 
	help = "Directory containing HLA allele files (with suffix .hla), if omitted assumed to be same as input-dir")

parser.add_argument("--output",
	type = str, 
	help = "Path to output file")


parser.add_argument("--quiet",
	type = str, 
	help = "Suppress INFO log messages")


MUTATION_FILE_EXTENSIONS = ["maf", "vcf"]

def find_mutation_files(input_dir_string):
	"""
	Collect all .vcf/.maf file paths in the dir(s) given as a comma-separated string.
	Returns a dictionary mapping base filenames to full paths. 
	"""
	mutation_files = OrderedDict()
	for dirpath in input_dir_string.split(","):
		for filename in listdir(dirpath):
			path = join(dirpath, filename)
			parts = filename.split(".")
			ext = parts[-1]
			patient_id = ".".join(parts[:-1]) 
			if ext in MUTATION_FILE_EXTENSIONS:
				logging.info("Reading mutation file %s", path)
				assert patient_id not in mutation_files, \
					"Duplicate files for %s: %s and %s" % (patient_id, mutation_files[patient_id], path)
				mutation_files[patient_id] = path
	return mutation_files

def find_hla_files(input_dir_string):
	"""
	Collect all .hla files  in the dir(s) given as a comma-separated string, 
	read in all the HLA alleles and normalize them. 

	Returns a dictionary mapping base filenames to lists of HLA allele names. 
	"""
	
	hla_types = {}
	for dirpath in input_dir_string.split(","):
		for filename in listdir(dirpath):
			parts = filename.split(".")
			ext = parts[-1]
			patient_id = ".".join(parts[:-1]) 
			if ext == 'hla':
				path = join(dirpath, filename)
				logging.info("Reading HLA file %s", path)
				assert patient_id not in hla_types, "Duplicate HLA files for %s" % patient_id 
				if patient_id not in mutation_files:
					logging.info("Skipping %s because don't have corresponding mutation file", path)
				else:
					alleles = []
					with open(path, 'r') as f:
						contents = f.read()
						for line in contents.split("\n"):
							for raw_allele_name in line.split(","):
								alleles.append(normalize_hla_allele_name(raw_allele_name))
					hla_types[patient_id] = alleles
	return hla_types

if __name__ == "__main__":

	args = parser.parse_args()

	init_logging(args.quiet)
	mutation_files = find_mutation_files(args.input_dir)

	# if no HLA input dir is specified then assume .hla files in the same dir
	# as the .maf/.vcf files 
	hla_dir_arg = args.hla_input_dir if args.hla_input_dir else args.input_dir
	hla_types = find_hla_files(hla_dir_arg)

	# make sure we have HLA types for each patient
	for patient_id, path in mutation_files.iteritems():
		assert patient_id in hla_types, "Missing HLA types for %s (%s)" % (patient_id, path)

	# dictionary that maps each patient ID to a tuple with three fields:
	# 	- number of mutated genes
	# 	- number of mutated genes with MHC binding mutated epitope
	# 	- number of mutated genes with immunogenic mutated epitope
	mutation_counts = {}
	for patient_id, path in mutation_files.iteritems():
		hla_allele_names = hla_types[patient_id]
		logging.info("Processing %s with HLA alleles %s", path, hla_allele_names)
		transcripts_df, raw_genomic_mutation_df, variant_report = load_file(path)

        # print each genetic mutation applied to each possible transcript
        # and either why it failed or what protein mutation resulted
        if not args.quiet:
            print_mutation_report(path, variant_report, raw_genomic_mutation_df, transcripts_df)

        print transcripts_df
    	mhc = PanBindingPredictor(hla_allele_names)
    	imm = ImmunogenicityPredictor(alleles = hla_allele_names)

    	scored_epitopes = mhc.predict(transcripts_df)
    	scored_epitopes = imm.predict(scored_epitopes)
    	print scored_epitopes
