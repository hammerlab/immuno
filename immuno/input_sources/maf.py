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


import gzip, logging, argparse, glob, re, pickle

import pandas
import Bio.SeqIO

from .. import common

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

SINGLE_AMINO_ACID_SUBSTITUTION = "p.([A-Z])([0-9]+)([A-Z])"
DELETION = "p.([A-Z])([0-9]+)del"

def extract(maf_df, refseq_map, epitope_max_length):
    half_max = int(epitope_max_length / 2)
    old_cols = [
        'Hugo_Symbol',
        'Chromosome',
        'Start_position',
        'Genome_Change',
        "Refseq_prot_Id",
        "Protein_Change",
        'Tumor_Sample_Barcode'
    ]
    filtered = maf_df[old_cols].dropna()
    new_cols = [
        'Hugo_Symbol',
        'Chromosome',
        'Start_position',
        'Genome_Change',
        'MAF Num',
        'Peptide',
        'Mutation Position',
        "Refseq Protein Id",
        'SampleId'
    ]
    result = dict((col, []) for col in cols)
    result['MAF Num'] = None
    for (index, row) in filtered.iterrows():
        ref_seq = refseq_map.get(row['Refseq_prot_Id'])
        if ref_seq is None:
            continue
        # TODO: Make this work for deletions, multi-residue indels, &c
        match = SINGLE_AMINO_ACID_SUBSTITUTION.match(row['Protein_Change'])
        if match is None:
            continue

        (wild_type, str_position, mutation) = match.groups()
        position = int(str_position) - 1

        if wild_type == mutation:
            continue

        if len(ref_seq.seq) <= position:
            logging.warning(
                "Mismatch in protein %s: ref is %d long, but mutation at %d",
                row['Refseq_prot_Id'],
                len(ref_seq.seq),
                position)
            continue

        if ref_seq.seq[position] != wild_type:
            logging.warning(
                "Mismatch in protein %s at pos %d: ref is %s but expected %s",
                row['Refseq_prot_Id'],
                position,
                ref_seq.seq[position],
                wild_type)
            continue

        # Make mutation
        full_peptide = ref_seq.seq.tomutable()
        full_peptide[position] = mutation
        start = max(0, position - half_max)
        stop = min(len(full_peptide), position + half_max + 1)
        trimmed_peptide = full_epitope[start:stop]
        position_in_epitope = position - max(0, position - half_max)
        assert(trimmed_epitope[position_in_peptide] == full_peptide[position])
        result['Peptide'].append(trimmed_epitope)
        result['Mutation Position'].append(position_in_epitope)
        result['Refseq Protein Id'].append(row['Refseq_prot_Id'])
        result['SampleId'].append(row['Tumor_Sample_Barcode'])
        result['Hugo_Symbol'].append(row['Hugo_Symbol'])
        result['Chromosome'].append(row['Chromosome'])
        result['Start_position'].append(row['Start_position'])
    result['Genome_Change'].append(row['Genome_Change'])
    return pandas.DataFrame(result, columns = new_cols)


def peptides_from_maf(
        maf_files,
        refseq_dir=['data/refseq'],
        epitope_max_length=31,
        nrows = None):
    refseq_filenames = []
    for dir in refseq_dir:
        refseq_filenames.extend(glob.glob("%s/*" % dir))
    refseq_map = refseq_id_to_sequence(refseq_filenames)
    result_dfs = []
    for (maf_num, maf_filename) in enumerate(maf_files):
        df = open_maf(maf_filename)
        result_df = extract(df, refseq_map, epitope_max_length)
        result_df['MAF Num'] = maf_num
        result_dfs.append(result_df)
    result_df = pandas.concat(result_dfs)
    return result_df[:nrows]

