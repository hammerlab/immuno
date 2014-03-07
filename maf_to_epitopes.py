#!/usr/bin/env python

# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""

Given a directory of refseq protein fasta files and a maf file, generate a csv file that
maps each maf file to a list of epitopes.

Example:

./maf_to_epitopes.py --refseq-dir data/refseq \ 
           --maf data/maf/step4_LUSC_Paper_v8.aggregated.tcga.maf2.4.migrated.somatic.maf \
           epitopes_raw.csv

TODO:
 - handle mutations other than simple missense point mutations
 
"""

import gzip, logging, argparse, glob, re, pickle

import pandas
import Bio.SeqIO

import common

def refseq_id_to_sequence(refseq_filenames):
  result = {}
  for filename in refseq_filenames:
    logging.info("loading: %s", filename)
    if filename.endswith(".gz"):
      fd = gzip.GzipFile(filename)
    else:
      fd = open(filename)
    gen = Bio.SeqIO.parse(fd, "fasta")
    for record in gen:
      try:
        name = record.id.split("|")[3]
        if name.startswith("NP_"):
          # TODO(odonnt02): handle multiple entries more intelligently than this.
          before_dot = name.split('.')[0]
          result[before_dot] = record
      except IndexError:
        pass
    fd.close()
  logging.info("loaded %d refseq sequences from %d files", len(result), len(refseq_filenames))
  return result

def open_maf(filename):
  print filename
  logging.info("Opening %s" % filename)
  with open(filename) as fd:
    lines_to_skip = 0
    while next(fd).startswith('#'):
        lines_to_skip += 1
  return pandas.read_csv(filename, skiprows=lines_to_skip, sep="\t", low_memory=False)

SINGLE_AMINO_ACID_SUBSTITUTION = re.compile("p.([A-Z])([0-9]+)([A-Z])")
def extract(maf_df, refseq_map, epitope_max_length):
  half_max = int(epitope_max_length / 2)
  # print maf_df.columns
  # print maf_df
  filtered = maf_df[['Hugo_Symbol', 'Chromosome', 'Start_position', 'Genome_Change', "Refseq_prot_Id", "Protein_Change", 'Tumor_Sample_Barcode']].dropna()
  cols = ['Hugo_Symbol', 'Chromosome', 'Start_position', 'Genome_Change', 'MAF Num', 'Epitope', 'Mutation Position', "Refseq Protein Id", 'SampleId']
  result = dict((col, []) for col in cols)
  result['MAF Num'] = None
  for (index, row) in filtered.iterrows():
    ref_seq = refseq_map.get(row['Refseq_prot_Id'])
    if ref_seq is None:
      continue
    #logging.debug("rfseq match %s", row['Refseq_prot_Id'])
    match = SINGLE_AMINO_ACID_SUBSTITUTION.match(row['Protein_Change'])
    if match is None:
      continue

    (wild_type, str_position, mutation) = match.groups()
    position = int(str_position) - 1

    if wild_type == mutation:
      continue

    if len(ref_seq.seq) <= position:
      logging.warning("Mismatch in protein %s: ref is only %d long, but mutation is at pos %d",
              row['Refseq_prot_Id'],
              len(ref_seq.seq),
              position)
      continue

    if ref_seq.seq[position] != wild_type:
      logging.warning("Mismatch in protein %s at pos %d: ref is %s but expected %s",
              row['Refseq_prot_Id'],
              position,
              ref_seq.seq[position],
              wild_type)
      continue

    # Make mutation
    full_epitope = ref_seq.seq.tomutable()
    full_epitope[position] = mutation
    trimmed_epitope = full_epitope[max(0, position - half_max) :
                                   min(len(full_epitope), position + half_max + 1)]
    position_in_epitope = position - max(0, position - half_max)
    assert(trimmed_epitope[position_in_epitope] == full_epitope[position])
    result['Epitope'].append(trimmed_epitope)
    result['Mutation Position'].append(position_in_epitope)
    result['Refseq Protein Id'].append(row['Refseq_prot_Id'])
    result['SampleId'].append(row['Tumor_Sample_Barcode'])

    result['Hugo_Symbol'].append(row['Hugo_Symbol'])
    result['Chromosome'].append(row['Chromosome'])
    result['Start_position'].append(row['Start_position'])
    result['Genome_Change'].append(row['Genome_Change'])

  return pandas.DataFrame(result, columns = cols)


def get_eptiopes_from_maf(maf_files, refseq_dir=['data/refseq'], epitope_max_length=17):
  refseq_filenames = []
  for dir in refseq_dir:
    print dir
    refseq_filenames.extend(glob.glob("%s/*" % dir))
  refseq_map = refseq_id_to_sequence(refseq_filenames)
  result_dfs = []
  for (maf_num, maf_filename) in enumerate(maf_files):
    df = open_maf(maf_filename)
    result_df = extract(df, refseq_map, epitope_max_length)
    result_df['MAF Num'] = maf_num
    result_dfs.append(result_df)
  result_df = pandas.concat(result_dfs)
  return result_df[:20000]

def go():
  parser = argparse.ArgumentParser(usage = __doc__)
  parser.add_argument("--refseq-file", action="append", default=[])
  parser.add_argument('--refseq-dir', action="append", default=[])
  parser.add_argument('--refseq-num', type=int, metavar="X",
    help="Read only first X refseq files (for debugging quickly).")
  parser.add_argument('--maf', action="append", default=[])
  parser.add_argument('--epitope-max-length', type=int, default=17)
  parser.add_argument("out")

  args = parser.parse_args()

  refseq_filenames = list(args.refseq_file)
  for dir in args.refseq_dir:
    refseq_filenames.extend(glob.glob("%s/*" % dir))
  if args.refseq_num:
    refseq_filenames = refseq_filenames[:args.refseq_num]
  refseq_map = refseq_id_to_sequence(refseq_filenames)

  result_dfs = []
  for (maf_num, maf_filename) in enumerate(args.maf):
    df = open_maf(maf_filename)
    result_df = extract(df, refseq_map, args.epitope_max_length)
    result_df['MAF Num'] = maf_num
    result_dfs.append(result_df)
  result_df = pandas.concat(result_dfs)
    
  logging.info("Matched %d epitopes", len(result_df))
  result_df.to_csv(args.out, index = False)
  logging.info("Wrote: %s", args.out)
    
if __name__ == "__main__":
  go()
  
