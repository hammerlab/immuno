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

"""
This script lets you count the number of immunogenic mutations in a collection
of cancer mutation files (in either MAF or VCF formats). Each mutation file is
expect to have a corresponding .hla file containing the patient's HLA alleles.
The output format is a CSV file with the following fields:

  - patient_id (extracted from the base of each VCF/MAF file)
  - number of coding mutations
  - number of coding mutations which contained MHC epitopes
  - number of coding mutations whose MHC epitopes are expected to be
    immunogenic

Example usage:
  python analyze_cohort.py --input-dir ../canseq/ --hla-input-dir\
          ../canseq-hla/ --output results.csv
"""

import argparse
import logging
from os import listdir
from os.path import join, split, splitext, isfile, abspath
import traceback
from collections import OrderedDict

import pandas as pd
import numpy as np

from common import (init_logging, splitext_permissive, find_paths)
from immunogenicity import ImmunogenicityPredictor
from load_file import (
    load_file, maf_to_vcf, expand_transcripts, load_variants
)
from maf import load_maf, get_patient_id, is_valid_tcga
from mhc_common import normalize_hla_allele_name
from hla_file import read_hla_file
from mhc_netmhcpan import PanBindingPredictor
from mhc_netmhccons import ConsensusBindingPredictor
from mutation_report import print_mutation_report

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--input-dir",
   type=str,
   help="Directory containing MAF or VCF input files")

group.add_argument("--input-file",
   type=str,
   help="Single MAF or VCF input file")

parser.add_argument("--hla-dir",
    type=str,
    default=None,
    help=("Directory containing HLA allele files (with suffix .hla)"))

parser.add_argument("--output",
    default="analyze_cohort_results.csv",
    help="Path to output file")

parser.add_argument("--quiet",
    type=str,
    help="Suppress INFO log messages")

parser.add_argument("--binding-threshold",
    type=int,
    default=500,
    help="Cutoff IC50 score for epitope MHC binding")

parser.add_argument("--combined-maf",
    default=False,
    action="store_true",
    help=("Rather than using filenames to identify patients, "
          "a single MAF file can have multiple tumor barcodes."))

parser.add_argument("--rna-filter-dir",
    type=str,
    default=None,
    help=("Directory containing RNASeq gene expression "
          "levels (one file per patient). If provided, we "
          "filter mutations with no gene expression."))

parser.add_argument("--debug-patient-id",
    type=str,
    default=None,
    help=("If we have a directory or a file containing "
          "multiple patient IDs, limit that collection to "
          "one specific patient ID for debugging."))

parser.add_argument("--debug-scored-epitopes-csv",
    type=str,
    default=None,
    help=("If we have a CSV file representing scored "
          "epitopes, use that instead of running netMHCpan. "
          "If not, generate that CSV file."))

parser.add_argument("--netmhc-cons",
    default=False,
    action="store_true",
    help="Use local NetMHCcons binding predictor (otherwise use NetMHCpan)")

parser.add_argument("--resume",
    default=False,
    action="store_true",
    help="Append to an existing output file"
)


MUTATION_FILE_EXTENSIONS = [".maf", ".vcf"]

def find_mutation_files(
        input_files, combined_maf=False, max_peptide_length=31):
    """
    Collect all .vcf/.maf file paths in the `input_filenames` list.

    Returns a dictionary mapping patient IDs to DataFrames containing basic
    variant information (chr, pos, ref, alt). The patient IDs will be each
    filename without its extension, unless the argument combined_maf is True.
    In this case, patient IDs are derived from the tumor barcode column in
    each MAF file.
    """
    mutation_files = OrderedDict()

    for path in input_filenames:
        _, filename = split(path)
        base, ext = splitext(filename)
        if ext in MUTATION_FILE_EXTENSIONS:
            if ext.endswith('maf') and combined_maf:
                maf_df = load_maf(path)
                file_patients = {}
                for barcode, group_df in (
                        maf_df.groupby(['Tumor_Sample_Barcode'])):
                    vcf_df = maf_to_vcf(group_df)
                    patient_id = get_patient_id(barcode)
                    file_patients[patient_id] = vcf_df
            else:
                patient_id = get_patient_id(base)
                vcf_df = load_variants(path)
                file_patients = {patient_id: vcf_df}

            for patient_id, vcf_df in file_patients.iteritems():
                patient_id = "-".join(patient_id.split("-")[:3])
                if patient_id in mutation_files:
                    logging.warning(
                        "Already processed patient %s before file %s",
                            patient_id,
                            path)
                else:
                    mutation_files[patient_id] = vcf_df
    if args.debug_patient_id:
        patient_id = args.debug_patient_id
        mutation_files = {patient_id: mutation_files[patient_id]}
    return mutation_files


def collect_hla_files(input_dir_string):
    return collect_files(input_dir_string, read_hla_file,
        permissive_parsing = True)


def collect_gene_exp_files(input_dir_string):
    return collect_files(input_dir_string, read_gene_exp_file,
        permissive_parsing = True)


def collect_files(input_dir_string, read_file_fn, permissive_parsing):
    """
    Collect all files in the dir(s) given as a comma-separated string,
    and then perform per-patient ID file_type-specific processing.
    """
    patient_to_data = {}
    for dirpath in input_dir_string.split(","):
        for filename in listdir(dirpath):
            base, ext = splitext_permissive(filename, [".txt"])
            if is_valid_tcga(base):
                patient_id = get_patient_id(base)
                path = join(dirpath, filename)
                result = read_file_fn(path, permissive_parsing)
                if result:
                    patient_to_data[patient_id] = result
    if args.debug_patient_id:
        patient_id = args.debug_patient_id
        patient_to_data = {patient_id: patient_to_data[patient_id]}
    return patient_to_data


def read_gene_exp_file(path, permissive_parsing):
    """
    Read in gene expression counts, returning a set of expressed genes.

    Expects the first column to be the gene name (or "<gene name>|<id>"),
    and the second column to be what we're filtering on.
    """
    logging.info("Reading gene expression file %s", path)
    gene_exp_df = pd.read_csv(path, sep='\t')
    gene_exp_df = gene_exp_df[gene_exp_df.columns[:2]]
    gene_col = gene_exp_df.columns[0]
    count_col = gene_exp_df.columns[1]
    if permissive_parsing:
        gene_exp_df[gene_col] = gene_exp_df[gene_col].str.split('|').map(
                lambda x: x[0])
    gene_exp_df = gene_exp_df[gene_exp_df[count_col] > 0]
    return set(gene_exp_df[gene_col].tolist())

def generate_mutation_counts(
        mutation_files,
        hla_types,
        genes_expressed,
        max_peptide_length=31,
        skip_identifiers = {},
        output_file=None):
    """
    Returns dictionary that maps each patient ID to a tuple with six fields:
        - total number of mutated epitopes across all transcripts
        - number of mutated genes
        - number of mutated genes with MHC binding mutated epitope
        - number of mutated epitopes which are predicted to bind to an MHC
          allele
        - number of mutated genes with at least one immunogenic mutated
          epitope
        - number of mutated epitopes which are predicted to be immunogenic
          (MHC binder + non-self)
    """
    mutation_counts = OrderedDict()
    n = len(mutation_files)
    for i, (patient_id, vcf_df) in enumerate(mutation_files.iteritems()):
        if patient_id in skip_identifiers:
            logging.info("Skipping patient ID %s", patient_id)
            continue
        hla_allele_names = hla_types[patient_id]
        logging.info(
            "Processing %s (#%d/%d) with HLA alleles %s",
            patient_id, i + 1, n, hla_allele_names)

        if not args.quiet:
            print vcf_df

        try:
            transcripts_df, raw_genomic_mutation_df, variant_report = (
                expand_transcripts(
                    vcf_df,
                    patient_id,
                    max_peptide_length=max_peptide_length))
        except KeyboardInterrupt:
            raise
        except:
            logging.warning("Failed to apply mutations for %s", patient_id)
            raise

        # print each genetic mutation applied to each possible transcript
        # and either why it failed or what protein mutation resulted
        if not args.quiet:
            print_mutation_report(
                patient_id,
                variant_report,
                raw_genomic_mutation_df,
                transcripts_df)
            logging.info(
                "Calling MHC binding predictor for %s (#%d/%d)",
                patient_id, i + 1, n)

        def make_mhc_predictor():
            if args.netmhc_cons:
                return ConsensusBindingPredictor(hla_allele_names)
            else:
                return PanBindingPredictor(hla_allele_names)

        # If we want to read scored_epitopes from a CSV file, do that.
        if args.debug_scored_epitopes_csv:
            csv_file = args.debug_scored_epitopes_csv
            if isfile(csv_file):
                scored_epitopes = pd.read_csv(csv_file)
            else:
                mhc = make_mhc_predictor()
                scored_epitopes = mhc.predict(transcripts_df,
                        mutation_window_size=9)
                scored_epitopes.to_csv(csv_file)
        else:
            mhc = make_mhc_predictor()
            scored_epitopes = mhc.predict(transcripts_df,
                    mutation_window_size=9)

        if not args.quiet:
            print scored_epitopes

        imm = ImmunogenicityPredictor(
            alleles=hla_allele_names,
            binding_threshold=args.binding_threshold)
        scored_epitopes = imm.predict(scored_epitopes)
        scored_epitopes.to_csv("scored_epitopes.csv")
        scored_epitopes = pd.read_csv("scored_epitopes.csv")

        grouped = scored_epitopes.groupby(["Gene", "GeneMutationInfo"])
        n_coding_mutations = len(grouped)
        n_epitopes = 0
        n_ligand_mutations = 0
        n_ligands = 0
        n_immunogenic_mutations = 0
        n_immunogenic_epitopes = 0
        for (gene, mut), group in grouped:
            start_mask = group.EpitopeStart < group.MutationEnd
            stop_mask = group.EpitopeEnd > group.MutationStart
            mutated_subset = group[start_mask & stop_mask]
            # we might have duplicate epitopes from multiple transcripts, so
            # drop them
            n_curr_epitopes = len(mutated_subset.groupby(['Epitope']))
            n_epitopes += n_curr_epitopes
            below_threshold_mask = \
                mutated_subset.MHC_IC50 <= args.binding_threshold
            ligands = mutated_subset[below_threshold_mask]
            n_curr_ligands = len(ligands.groupby(['Epitope']))
            n_ligands += n_curr_ligands
            n_ligand_mutations += (n_curr_ligands) > 0
            thymic_deletion_mask = \
                np.array(ligands.ThymicDeletion).astype(bool)
            immunogenic_epitopes = ligands[~thymic_deletion_mask]
            curr_immunogenic_epitopes = immunogenic_epitopes.groupby(['Epitope']).first()
            n_immunogenic_epitopes += len(curr_immunogenic_epitopes)
            n_immunogenic_mutations += len(curr_immunogenic_epitopes) > 0
            logging.info(("%s %s: epitopes %s, ligands %d, imm %d"),
                         gene,
                         mut,
                         n_curr_epitopes,
                         n_curr_ligands,
                         len(curr_immunogenic_epitopes),
                        )
        result_tuple = (
            n_coding_mutations,
            n_epitopes,
            n_ligand_mutations,
            n_ligands,
            n_immunogenic_mutations,
            n_immunogenic_epitopes,
        )
        if output_file:
            data_string = ",".join(str(d) for d in result_tuple)
            output_file.write("%s,%s\n" % (patient_id, data_string))
            output_file.flush()
        mutation_counts[patient_id] = result_tuple
    return mutation_counts

if __name__ == "__main__":
    args = parser.parse_args()

    init_logging(args.quiet)

    # collect input files from commandline arguments specifying either
    # filenames or directories
    input_filenames = []
    if args.input_file:
        for filename in args.input_file.split(","):
            input_filenames.append(abspath(filename))
    if args.input_dir:
        for dirpath in args.input_dir.split(","):
            for filename in listdir(dirpath):
                path = join(dirpath, filename)
                input_filenames.append(path)
    mutation_files = find_mutation_files(input_filenames, args.combined_maf)

    # if no HLA input dir is specified then assume .hla files in the same dir
    # as the .maf/.vcf files
    hla_dir_arg = args.hla_dir if args.hla_dir else args.input_dir
    assert hla_dir_arg, "Specify HLA directory via --hla-dir argument"
    hla_types = collect_hla_files(hla_dir_arg)

    # If no RNA dir is specified, then assume we are not doing RNA filtering
    genes_expressed = None
    if args.rna_filter_dir:
        genes_expressed = collect_gene_exp_files(args.rna_filter_dir)

    missing = set([])
    # make sure we have HLA types for each patient
    for patient_id in mutation_files.iterkeys():
        if patient_id not in hla_types:
            missing.add(patient_id)
    if len(missing) > 0:
        logging.warning("Missing HLA types for %s", list(missing))

    logging.info(
        "Total missing HLA types: %d / %d patients",
        len(missing),
        len(mutation_files))

    for patient_id in missing:
        del mutation_files[patient_id]

    if args.resume:
        with open(args.output, 'r') as f:
            lines = [l for l in f.read().split("\n") if len(l) > 0]
            fields = [l.split(",") for l in lines]
            # expect ID and 6 data fields
            # warning: the number of fields is tightly
            # coupled with the code that actually writes the CSV files
            complete_fields = [l for l in fields if len(l) == 7]
        output_file = open(args.output, 'w')
        finished_identifiers = []
        for l in complete_fields:
            identifier = l[0]
            finished_identifiers.append(identifier)
            output_file.write("%s\n" % ",".join(l))
        output_file.flush()
    else:
        finished_identifiers = []
        output_file = open(args.output, 'w')

    mutation_counts = generate_mutation_counts(
        mutation_files,
        hla_types,
        genes_expressed,
        skip_identifiers = finished_identifiers,
        output_file=output_file)

    output_file.close()

    print
    print "SUMMARY"
    for patient_id, fields in mutation_counts.iteritems():
        (
            n_coding_mutations, n_epitopes,
            n_ligand_mutations, n_ligands,
            n_immunogenic_mutations, n_immunogenic_epitopes,
        ) = fields
        print(
            ("%s: # mutations %d (%d epitopes), # mutations with ligands "
             "%d (%d epitopes), # immunogenic mutations %d (%d epitopes) "
            ) % (
                patient_id,
                n_coding_mutations,
                n_epitopes,
                n_ligand_mutations,
                n_ligands,
                n_immunogenic_mutations,
                n_immunogenic_epitopes,
            ))