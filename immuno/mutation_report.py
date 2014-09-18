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
from mako.template import Template
from mako.lookup import TemplateLookup

from common import peptide_substrings, init_logging
from mhc_iedb import IEDBMHCBinding, normalize_hla_allele_name
from mhc_netmhcpan import PanBindingPredictor
from mhc_netmhccons import ConsensusBindingPredictor
import mhc_random 
from load_file import load_file
from strings import load_comma_string
from immunogenicity import ImmunogenicityPredictor
from vaccine_peptides import build_peptides_dataframe

DEFAULT_ALLELE = 'HLA-A*02:01'

parser = argparse.ArgumentParser()
# must supply either an input file or an amino acid string
input_group = parser.add_argument_group()

input_group.add_argument("--input-file", 
    action="append", 
    default=[],
    help="input file(s) (must be FASTA, MAF, TAB or VCF format)")

input_group.add_argument("--string", 
    default=None,
    help="Literal amino acid string of mutated peptide")

parser.add_argument("--quiet",
    default=False, 
    action="store_true",
    help="Suppress verbose output"
)

parser.add_argument("--peptide-length",
    default=31,
    type=int,
    help="Length of vaccine peptides (may contain multiple epitopes)")

parser.add_argument("--min-peptide-padding", 
    default=0, 
    type=int, 
    help="Minimum number of wildtype residues before or after a mutation")

parser.add_argument("--hla-file",
    help="File with one HLA allele per line")

parser.add_argument(
    "--hla",
    help="Comma separated list of allele (default HLA-A*02:01)")

parser.add_argument("--random-mhc",
    default=False,
    action="store_true",
    help="Random values instead for MHC binding prediction")

parser.add_argument("--iedb-mhc",
    default=False,
    action="store_true",
    help="Use IEDB's web API for MHC binding")

parser.add_argument("--netmhc-cons",
    default=False,
    action="store_true",
    help="Use local NetMHCcons binding predictor")

parser.add_argument("--output-epitopes-file",
    help="Output CSV file for dataframe containing scored epitopes",
    required=False)

parser.add_argument("--print-epitopes",
    help="Print dataframe with epitope scores",
    default=False,
    action="store_true")

parser.add_argument("--print-peptides",
    default=False, 
    help="Print dataframe with vaccine peptide scores",
    action="store_true")

parser.add_argument("--output-report-file",
    default="report.html",
    help="Path to HTML report containing scored vaccine peptides")

def print_mutation_report(
        input_filename,
        variant_report,
        raw_genomic_mutation_df,
        transcripts_df):
    print 
    print "MUTATION REPORT FOR", input_filename 
    print 
    last_mutation = None 
    for (mut_description, transcript_id), msg in variant_report.iteritems():
        if mut_description != last_mutation:
            print mut_description
            last_mutation = mut_description
        print "--", transcript_id, ":", msg 

    logging.info("---")
    logging.info("FILE LOADING SUMMARY FOR %s", input_filename)
    logging.info("---")
    logging.info("# original mutations: %d", len(raw_genomic_mutation_df))
    logging.info(
        "# mutations with annotations: %d", 
        len(transcripts_df.groupby(['chr', 'pos', 'ref', 'alt'])))
    logging.info("# transcripts: %d", len(transcripts_df))

def group_epitopes(scored_epitopes):
    """
    Given a DataFrame with fields:
      - TranscriptId
      - SourceSequence
      - MutationStart
      - MutationEnd
      - GeneMutationInfo
      - PeptideMutationInfo 
      - Gene
      - GeneInfo
      - Epitope
      - EpitopeStart
      - EpitopeEnd
      - MHC_IC50
      - MHC_PercentileRank
    
    Group epitopes under their originating transcript and
    make nested lists of dictionaries to contain the binding scores
    for MHC alleles. 
    
    Return a list of dictionaries for each mutated transcript with fields:
      - TranscriptId
      - SourceSequence
      - MutationStart
      - MutationEnd
      - GeneMutationInfo
      - PeptideMutationInfo 
      - Gene
      - GeneInfo
      - Epitopes : list of dictionaries
      
    Each entry of the 'Epitopes' list contains the following fields:
      - 'Epitope'
      - 'EpitopeStart'
      - 'EpitopeEnd'
      - 'MHC_AlleleScores' : list of allele-specific entries 
    
    Each entry of 'MHC_AlleleScores' the following fields:
      - 'Allele' 
      - 'MHC_PercentileRank'
      - 'MHC_IC50' 
    """
    peptides = []
    for (transcript_id, seq), transcript_group in \
            scored_epitopes.groupby(["TranscriptId", "SourceSequence"]):
        peptide_entry = {}
        peptide_entry["Peptide"] = seq
        peptide_entry['TranscriptId'] = transcript_id
        head = transcript_group.to_records()[0]
        peptide_entry["MutationStart"] = head.MutationStart
        peptide_entry["MutationEnd"] = head.MutationEnd
        peptide_entry["GeneMutationInfo"] = head.GeneMutationInfo
        peptide_entry["PeptideMutationInfo"] = head.PeptideMutationInfo
        peptide_entry["GeneInfo"] = head.GeneInfo
        peptide_entry['Gene'] = head.Gene
        peptide_entry['Description'] = "%s (%s) : %s" % (
                head.Gene, head.TranscriptId, head.PeptideMutationInfo
        ) 
        peptide_entry['Epitopes'] = []
        for (epitope, epitope_start, epitope_end), epitope_group in \
                transcript_group.groupby(
                    ['Epitope', 'EpitopeStart', 'EpitopeEnd']):
            epitope_entry = {
                'Epitope' : epitope, 
                'EpitopeStart' : epitope_start, 
                'EpitopeEnd' : epitope_end, 
                'MHC_Allele_Scores' : []
            }
            seen_alleles = set([])
            for epitope_allele_row in epitope_group.to_records():
                allele = epitope_allele_row['Allele']
                assert allele not in seen_alleles, \
                    "Repeated entry %s" % epitope_allele_row
                seen_alleles.add(allele)
                percentile_rank = epitope_allele_row['MHC_PercentileRank']
                ic50 = epitope_allele_row['MHC_IC50']
                allele_entry = {
                    'Allele': allele, 
                    'MHC_PercentileRank' : percentile_rank,
                    'MHC_IC50' : ic50,
                }
                epitope_entry['MHC_Allele_Scores'].append(allele_entry)
            peptide_entry['Epitopes'].append(epitope_entry)
        peptides.append(peptide_entry)
    return peptides

def render_report(input_names, peptides, alleles):
    template_lookup = TemplateLookup(
        directories=['.', 'viz'], 
        default_filters=['literal']
    )
    template = Template(
        filename = 'viz/index.html.template', 
        lookup = template_lookup
    )

    html = template.render(
        peptides = peptides, 
        vcf_filename = input_names, 
        hla_alleles = alleles,
    )
    return html 
    


if __name__ == '__main__':
    args = parser.parse_args()

    init_logging(args.quiet)

    peptide_length = int(args.peptide_length)

    # get rid of gene descriptions if they're in the dataframe
    if args.hla_file:
        alleles = [normalize_hla_allele_name(l) for l in open(args.hla_file)]
    elif args.hla:
        alleles = [normalize_hla_allele_name(l) for l in args.hla.split(",")]
    else:
        alleles = [normalize_hla_allele_name(DEFAULT_ALLELE)]

    # stack up the dataframes and later concatenate in case we
    # want both commandline strings (for weird mutations like translocations)
    # and files
    mutated_region_dfs = []

    if args.string:
        df = load_comma_string(args.string)
        mutated_region_dfs.append(df)

    # loop over all the input files and
    # load each one into a dataframe

    for input_filename in args.input_file:
        transcripts_df, raw_genomic_mutation_df, variant_report = \
            load_file(input_filename, max_peptide_length = peptide_length)
        mutated_region_dfs.append(transcripts_df)

        # print each genetic mutation applied to each possible transcript
        # and either why it failed or what protein mutation resulted
        if not args.quiet:
            print_mutation_report(
                input_filename, 
                variant_report, 
                raw_genomic_mutation_df, 
                transcripts_df)
        
    if len(mutated_region_dfs) == 0:
        parser.print_help()
        print "\nERROR: Must supply at least --string or --input-file"
        sys.exit()

    mutated_regions = pd.concat(mutated_region_dfs)

    if args.random_mhc:
        scored_epitopes = \
            mhc_random.generate_scored_epitopes(mutated_regions, alleles)
    elif args.iedb_mhc:
        mhc = IEDBMHCBinding(name = 'mhc', alleles=alleles)
        scored_epitopes = mhc.predict(mutated_regions)
    elif args.netmhc_cons:
        predictor = ConsensusBindingPredictor(alleles)
        scored_epitopes = predictor.predict(mutated_regions)
    else:
        predictor = PanBindingPredictor(alleles)
        scored_epitopes = predictor.predict(mutated_regions)

    imm = ImmunogenicityPredictor(alleles = alleles)
    scored_epitopes = imm.predict(scored_epitopes)


    if 'MHC_PercentileRank' in scored_epitopes:
        scored_epitopes = scored_epitopes.sort(['MHC_PercentileRank'])

    if args.output_epitopes_file:
        scored_epitopes.to_csv(args.output_epitopes_file, index=False)
        
    if args.print_epitopes:
        print scored_epitopes.to_string()

    peptides = group_epitopes(scored_epitopes)
   

    if args.print_peptides:
        for pep in peptides:
            print pep 
    
    
    input_names = ";".join(args.input_file)
    if args.string:
        input_names += ";" + args.string
    
    html = render_report(input_names, peptides, alleles)
    with open(args.output_report_file, 'w') as f:
        f.write(html)
