"""

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""
import argparse

import pandas as pd
from epitopes import reduced_alphabet 

from pipeline import ImmunoPipeline
from immunogenicity import ImmunogenicityRFModel
from binding import IEDBMHCBinding
from cleavage import ProteasomalCleavage
from Bio import SeqIO
from maf_to_epitopes import get_eptiopes_from_maf

def get_epitopes_from_fasta(fasta_files):
  epitopes = []
  for fasta_file in fasta_files:
    epitope_data = SeqIO.parse(fasta_file, 'fasta')
    epitopes += [pd.DataFrame({'peptide':  pd.Series(list(set([e.seq for e in epitope_data])))})]
  return pd.concat(epitopes)


def create_pipeline():
  pipeline = ImmunoPipeline()
  return pipeline

def add_scoring(pipeline, alleles):
  pipeline.add_scorer(IEDBMHCBinding(name='mhc', alleles=alleles))
  pipeline.add_scorer(ProteasomalCleavage(name='protocleave'))
  #pipeline.add_scorer(SelfScorer(name='selfcheck'))
  pipeline.add_scorer(ImmunogenicityRFModel(name='default RF'))
  pipeline.add_scorer(ImmunogenicityRFModel(name = 'murphy10 RF', reduced_alphabet = reduced_alphabet.murphy10))

  return pipeline

DEFAULT_ALLELE = 'HLA-A*02:01'

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", action="append", default=[], help="input file to process")
  parser.add_argument("--allele_file", help="list of alleles")
  parser.add_argument("--output", help="output file for dataframes", required=True)


  args = parser.parse_args()
  if args.input[0].endswith(".vcf"):
    epitope_data = pipeline.add(Variant2Epitope)
  elif args.input[0].endswith(".maf"):
    epitope_data = get_eptiopes_from_maf(args.input)
  elif args.input[0].endswith(".fasta") or args.input[0].endswith(".fa"):
    epitope_data =  get_epitopes_from_fasta(args.input)

  if args.allele_file:
    alleles = [l.strip() for l in open(allele_file)]
  else:
    alleles = [DEFAULT_ALLELE]
  pipeline = ImmunoPipeline()
  add_scoring(pipeline, alleles)
  data = pipeline.score(epitope_data)

  data.to_csv(args.output, index=False)
