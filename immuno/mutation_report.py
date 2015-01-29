#!/usr/bin/env python2

# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import argparse
import sys

import pandas as pd
from Bio import SeqIO
import numpy as np
from varcode import VariantCollection

from common import peptide_substrings, init_logging, parse_int_list
from epitope_scoring import (
        simple_ic50_epitope_scorer,
        logistic_ic50_epitope_scorer,
)
from group_epitopes import group_epitopes_dataframe
from immunogenicity import (ImmunogenicityPredictor, THYMIC_DELETION_FIELD_NAME)
from load_file import load_file
from mhc_common import normalize_hla_allele_name
from mhc_iedb import IEDB_MHC1
from mhc_netmhcpan import PanBindingPredictor
from mhc_netmhccons import ConsensusBindingPredictor
import mhc_random
from peptide_binding_measure import IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
from strings import load_comma_string
from vaccine_peptides import select_vaccine_peptides

DEFAULT_ALLELE = 'HLA-A*02:01'

parser = argparse.ArgumentParser()

parser.add_argument("--variants-file",
    help="Input file with DNA variants (must be MAF, TAB or VCF format)")


parser.add_argument("--quiet",
    default=False,
    action="store_true",
    help="Suppress verbose output"
)


###
# MHC options
###

mhc_arg_parser = parser.add_argument_group(
    title="MHC",
    description="Which MHC binding predictor to use (default NetMHCpan)")

mhc_arg_parser.add_argument("--hla-file",
    help="File with one HLA allele per line")

mhc_arg_parser.add_argument("--hla",
    help="Comma separated list of allele (default HLA-A*02:01)")

mhc_arg_parser.add_argument("--random-mhc",
    default=False,
    action="store_true",
    help="Random values instead for MHC binding prediction")

mhc_arg_parser.add_argument("--iedb-mhc",
    default=False,
    action="store_true",
    help="Use IEDB's web API for MHC binding")

mhc_arg_parser.add_argument("--netmhc-cons",
    default=False,
    action="store_true",
    help="Use local NetMHCcons binding predictor")


mhc_arg_parser.add_argument("--mhc-epitope-lengths",
    default=[8,9,10,11],
    type=parse_int_list,
    help="Lengths of epitopes to consider for MHC binding prediction")

parser.add_argument("--self-epitope-filter",
    default=False,
    action="store_true",
    help="Filter 'self' epitopes (predicted MHC binders from the proteome)")

# TODO: Add a '--self-epitope-filter-path' option
# for directory with following structure:
# /dir/9/HLA-A0201
# ...
# /dir/11/HLA-C0801
# ....
# /dir/length/allele

parser.add_argument("--output-epitopes-file",
    help="Output CSV file for dataframe containing scored epitopes",
    required=False)

parser.add_argument("--print-epitopes",
    help="Print dataframe with epitope scores",
    default=False,
    action="store_true")


###
# RNA-Seq data
###
rna_group = parser.add_argument_group(
    title="RNA-Seq",
    description="Transcript and gene abundance quantification")

rna_group.add_argument(
    "--rna-gene-fpkm-file",
    help="Tab separated mapping from Ensembl gene IDs to FPKM")

rna_group.add_argument(
    "--rna-transcript-fpkm-file",
    help="Tab separated mapping from Ensembl transcript IDs to FPKM")

###
# Vaccine peptide options
###

vaccine_peptide_arg_parser = parser.add_argument_group(
    title="Vaccine Peptides",
    description="Options affecting selection and display of vaccine peptides")

vaccine_peptide_arg_parser.add_argument("--vaccine-peptide-file",
    default="vaccine-peptides.csv",
    help="Path to CSV file containing predicted vaccine peptides")

vaccine_peptide_arg_parser.add_argument("--vaccine-peptide-count",
    default=None,
    type=int,
    help="How many vaccine peptides do we save? (default: all mutations)")

vaccine_peptide_arg_parser.add_argument(
    "--vaccine-peptide-logistic-epitope-scoring",
    default=False,
    action="store_true",
    help="Use continuous score per epitope (instead of just IC50 <= 500nM)")

vaccine_peptide_arg_parser.add_argument("--vaccine-peptide-length",
    default=31,
    type=int,
    help="Length of vaccine peptides (may contain multiple epitopes)")

vaccine_peptide_arg_parser.add_argument("--vaccine-peptide-padding",
    default=5,
    type=int,
    help="Minimum number of wildtype residues before or after a mutation")

vaccine_peptide_arg_parser.add_argument("--print-vaccine-peptides",
    default=False,
    help="Print vaccine peptides and scores",
    action="store_true")

def load_variants_file(filename):
    variants = VariantCollection(args.variants_file)
    if len(variants) == 0:
        raise ValueError("No variants found in %s" % (filename,))

    logging.info("Loaded variants file %s" % (filename,))

    for variant in variants:
        logging.info(variant)
    return variants

def load_expression_levels(filename, key_field):
    df = pd.read_csv(filename, sep='\t')
    return df[[key_field, 'FPKM']]

def create_vaccine_fasta_lines(vaccine_peptide_records):
    fasta_lines = []
    for i, record in enumerate(vaccine_peptide_records):
        line = ">%d Gene=%s, Transcript=%s, Mut=%s (%d:%d), Score=%0.6f" % (
            i,
            record['Gene'],
            record['TranscriptId'],
            record['PeptideMutationInfo'],
            record['VaccinePeptideMutationStart'],
            record['VaccinePeptideMutationEnd'],
            record['MutantEpitopeScore']
        )
        fasta_lines.append(line)
        fasta_lines.append(record['VaccinePeptide'])
    return fasta_lines

def print_epitopes(source_sequences):
    print
    print "Epitopes"
    print "--------"
    print
    for record in sorted(source_sequences, key=lambda r: r['TranscriptId']):
        mut_start = record['MutationStart']
        mut_end = record['MutationEnd']
        print ">Gene=%s, Transcript=%s, Mut=%s (%d:%d)" % (
            record['Gene'],
            record['TranscriptId'],
            record['PeptideMutationInfo'],
            mut_start,
            mut_end,
        )
        print record['SourceSequence']
        for epitope in sorted(
                record['Epitopes'], key=lambda e: e['EpitopeStart']):
            overlap_start = epitope['EpitopeStart'] < mut_end
            overlap_end = epitope['EpitopeEnd'] > mut_start
            mutant = overlap_start and overlap_end
            print "\t", epitope['Epitope'], ("<-- MUTANT" if mutant else "")
            for prediction in sorted(
                    epitope["MHC_Allele_Scores"],
                    key=lambda p: p['Allele']):
                print "\t\t", "allele = %s, IC50=%0.4f" % (
                    prediction['Allele'],
                    prediction[IC50_FIELD_NAME],
                )

def mhc_binding_prediction(mutated_regions, alleles):
    if args.random_mhc:
        return mhc_random.generate_scored_epitopes(mutated_regions, alleles)
    elif args.iedb_mhc:
        mhc = IEDB_MHC1(alleles=alleles)
        return mhc.predict(mutated_regions)
    elif args.netmhc_cons:
        predictor = ConsensusBindingPredictor(alleles)
        return predictor.predict(mutated_regions)
    else:
        predictor = PanBindingPredictor(alleles)
        return predictor.predict(mutated_regions)

if __name__ == '__main__':
    args = parser.parse_args()
    init_logging(args.quiet)

    vaccine_peptide_length = int(args.vaccine_peptide_length)

    # get rid of gene descriptions if they're in the dataframe
    if args.hla_file:
        alleles = [normalize_hla_allele_name(l) for l in open(args.hla_file)]
    elif args.hla:
        alleles = [normalize_hla_allele_name(l) for l in args.hla.split(",")]
    else:
        alleles = [normalize_hla_allele_name(DEFAULT_ALLELE)]


    variants = load_variants_file(args.variants_file)
    gene_levels = load_expression_levels(args.rna_gene_fpkm_file)
    transcript_levels = load_expression_levels(args.rna_transcript_fpkm_file)

    # dictionary mapping each variant to coding effect
    # some variants may be dropped if they affect insufficiently
    # transcribed genes or have no coding effects
    variant_transcript_effects = select_transcript_effects(
        variants,
        gene_levels,
        transcript_levels,
        gene_fpkm_cutoff=1.0)

    mutated_regions = mutant_amino_acid_sequences(
        variant_transcript_effects,
        vaccine_peptide_length)

    # TODO: generate mutated_regions
    # mutated_regions = pd.concat(mutated_region_dfs)
    scored_epitopes = mhc_binding_prediction(mutated_regions, alleles)

    if args.self_epitope_filter:
        imm = ImmunogenicityPredictor(alleles = alleles)
        scored_epitopes = imm.predict(scored_epitopes)
    else:
        scored_epitopes[THYMIC_DELETION_FIELD_NAME] = False

    if PERCENTILE_RANK_FIELD_NAME in scored_epitopes:
        scored_epitopes = scored_epitopes.sort([PERCENTILE_RANK_FIELD_NAME])

    if args.output_epitopes_file:
        scored_epitopes.to_csv(args.output_epitopes_file, index=False)

    source_sequences = group_epitopes_dataframe(scored_epitopes)

    if args.print_epitopes:
        print_epitopes(source_sequences)

    if args.print_vaccine_peptides or args.vaccine_peptide_file:
        padding = args.vaccine_peptide_padding
        if args.vaccine_peptide_logistic_epitope_scoring:
            epitope_scorer = logistic_ic50_epitope_scorer
        else:
            epitope_scorer = simple_ic50_epitope_scorer

        vaccine_peptide_records = select_vaccine_peptides(
            source_sequences,
            epitope_scorer=epitope_scorer,
            vaccine_peptide_length=peptide_length,
            padding=padding)

        if args.vaccine_peptide_count:
            n = args.vaccine_peptide_count
            vaccine_peptide_records = vaccine_peptide_records[:n]

        fasta_lines = create_vaccine_fasta_lines(vaccine_peptide_records)

        if args.print_vaccine_peptides:
            for line in fasta_lines:
                print line

        if args.vaccine_peptide_file:
            with open(args.vaccine_peptide_file, 'w') as f:
                for line in fasta_lines:
                    f.write(line)
                    f.write("\n")
