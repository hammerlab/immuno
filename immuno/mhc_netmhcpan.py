
import subprocess 
import tempfile
import os

import pandas as pd 

from mhc_common import normalize_hla_allele_name

class PanBindingPredictor(object):
	def __init__(self, hla_alleles):
		valid_alleles_str = subprocess.check_output(["netMHCpan", "-listMHC"])
		valid_alleles = set([])
		for line in valid_alleles_str.split("\n"):
			if not line.startswith("#"):
				valid_alleles.add(line)
		print "Total alleles available for NetMHCpan: %d" % len(valid_alleles)
		self.alleles = []
		for allele in hla_alleles:
			allele = normalize_hla_allele_name(allele.strip().upper())
			# for some reason netMHCpan drop the "*" in names such as "HLA-A*03:01" becomes "HLA-A03:01"
			if  allele.replace("*", "") not in valid_alleles:
				print "Skipping %s (not available in NetMHCpan)" % allele
			else:
				self.alleles.append(allele)


	def predict(self, df):
		"""
		Given a dataframe of mutated amino acid sequences, run each sequence through NetMHCpan
		"""

		results = []

		print "Calling netMHCpan..."
		for mutation_entry in df.to_records():
			input_file = tempfile.NamedTemporaryFile("w", prefix="peptide", delete=False)
			input_file.write(">seq\n")
			input_file.write(mutation_entry['SourceSequence'])
			input_file.close()	
			

			output_file = tempfile.NamedTemporaryFile("w", prefix="netMHCpan_output", delete=False)
			output_file.close()

			print ">", mutation_entry.SourceSequence
			for allele in self.alleles:
				print "--", allele 
				subprocess.check_call(["netMHCpan", "-a", allele.replace("*", ""), "-l", "9", "-xls", "-xlsfile", output_file.name, "-f", input_file.name, ])
				epitopes_df = pd.read_csv(output_file.name, sep='\t', skiprows = 1)
				# columns: Pos    Peptide   ID  1-log50k          nM  Rank
				for epitope_row in epitopes_df.to_records():
					pos = epitope_row['Pos']
					epitope = epitope_row['Peptide']
					ic50 = epitope_row['nM']
					rank = epitope_row['Rank']
					new_row = {}
					# fields shared by all epitopes 
					new_row['SourceSequence'] = mutation_entry.SourceSequence
					new_row['MutationStart'] = mutation_entry.MutationStart
					new_row['MutationEnd'] = mutation_entry.MutationEnd
					new_row['GeneInfo'] = mutation_entry.GeneInfo
					new_row['Gene'] = mutation_entry.Gene
					new_row['MutationInfo'] = mutation_entry.MutationInfo
					new_row['TranscriptId'] = mutation_entry.TranscriptId

					# fields specific to this epitope 
					new_row['Allele'] = allele
					new_row['EpitopeStart'] = pos 
					new_row['EpitopeEnd'] = pos + 9
					new_row['Epitope'] = epitope 
					new_row['MHC_IC50'] = ic50
					new_row['MHC_PercentileRank'] = rank 
					results.append(new_row)
			os.remove(input_file.name)
			os.remove(output_file.name)
		return pd.DataFrame.from_records(results)