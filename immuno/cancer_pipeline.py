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

import argparse
import sys

import pandas as pd
from Bio import SeqIO
import numpy as np

from common import peptide_substrings
from immunogenicity import ImmunogenicityRFModel
from mhc import IEDBMHCBinding, normalize_hla_allele_name
from load_file import load_file
from strings import load_comma_string
from vaccine_peptides import build_peptides_dataframe
from mako.template import Template
from mako.lookup import TemplateLookup

DEFAULT_ALLELE = 'HLA-A*02:01'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # must supply either an input file or an amino acid string
    input_group = parser.add_argument_group()
    input_group.add_argument("--input", action="append", default=[],
        help="input file name (i.e. FASTA, MAF, VCF)")
    input_group.add_argument("--string", default = None,
        help="amino acid string")
    parser.add_argument("--peptide-length",
        default=31,
        type = int,
        help="length of vaccine peptides (may contain multiple epitopes)")
    parser.add_argument("--min-peptide-padding", 
        default = 5, 
        type = int, 
        help = "minimum number of wildtype residues before or after a mutation")
    parser.add_argument("--hla-file",
        help="file with one HLA allele per line")
    parser.add_argument(
        "--hla",
        help="comma separated list of allele (default HLA-A*02:01)")
    parser.add_argument(
        "--epitopes-output",
        help="output file for dataframe containing scored epitopes",
        required=False)
    parser.add_argument(
        "--peptides-output",
        help="output file for dataframe containing scored vaccine peptides",
        required=False)
    parser.add_argument(
        "--print-epitopes",
        help="print dataframe with epitope scores",
        default=False,
        action="store_true")
    parser.add_argument(
        "--print-peptides",
        default = False, 
        help="print dataframe with vaccine peptide scores",
        action="store_true")
    parser.add_argument(
        "--html-report",
        default = "report.html",
        help = "Path to HTML report containing scored peptides and epitopes")
    parser.add_argument("--skip-mhc",
        default=False,
        action="store_true",
        help="Don't predict MHC binding")

    args = parser.parse_args()


    peptide_length = int(args.peptide_length)


    # get rid of gene descriptions if they're in the dataframe
    if args.hla_file:
        alleles = [normalize_hla_allele_name(l) for l in open(args.hla_file)]
    elif args.hla:
        alleles = [normalize_hla_allele_name(l) for l in args.hla.split(",")]
    else:
        alleles = [DEFAULT_ALLELE]



    # stack up the dataframes and later concatenate in case we
    # want both commandline strings (for weird mutations like translocations)
    # and files
    mutated_region_dfs = []

    if args.string:
        df = load_comma_string(args.string)
        mutated_region_dfs.append(df)


    # loop over all the input files and
    # load each one into a dataframe

    for input_filename in args.input:
        df = load_file(input_filename, peptide_length)
        assert df is not None
        mutated_region_dfs.append(df)

    if len(mutated_region_dfs) == 0:
        parser.print_help()
        print "\nERROR: Must supply at least --string or --input"
        sys.exit()
    mutated_regions = pd.concat(mutated_region_dfs)


    if args.skip_mhc:
        records = []
        # if wer'e not running the MHC prediction then we have to manually
        # extract 9mer substrings
        for _, row in mutated_regions.iterrows():
            seq = row.SourceSequence
            epitope_length = 9
            for i in xrange(len(seq) - epitope_length):
                record = dict(row)
                record['Epitope'] = seq[i:i+epitope_length]
                record['EpitopeStart'] = i
                record['EpitopeEnd'] = i + epitope_length
                records.append(record)
        scored_epitopes = pd.DataFrame.from_records(records)
    else:
        mhc = IEDBMHCBinding(name = 'mhc', alleles=alleles)
        scored_epitopes = mhc.apply(mutated_regions)

        assert 'percentile_rank' in scored_epitopes, scored_epitopes.head()
        mhc_percentile = scored_epitopes['percentile_rank']
        mhc_score = (100.0 - mhc_percentile) / 100.0
        mhc_binding_category = mhc_percentile <= 1.0
        scored_epitopes['mhc_score'] = mhc_score
        scored_epitopes['mhc_binding_category'] = mhc_binding_category

    if 'mhc_score' in scored_epitopes:
        scored_epitopes = scored_epitopes.sort(['mhc_score'])

    if args.epitopes_output:
        scored_epitopes.to_csv(args.epitopes_output, index=False)
    if args.print_epitopes:
        print scored_epitopes.to_string()

    peptides = build_peptides_dataframe(scored_epitopes,
        peptide_length = peptide_length, 
        min_peptide_padding = args.min_peptide_padding)
    if args.peptides_output:
        peptides.to_csv(args.peptides_output, index=False)
    
    if args.print_peptides:
        print peptides.to_string()

    input_names = ";".join(args.input)
    if args.string:
        input_names += ";" + args.string
    template_lookup = TemplateLookup(directories=['.', 'viz'], default_filters=['literal'])
    template = Template(filename = 'viz/index.html.template', lookup = template_lookup)

    html = template.render(peptides = peptides) #(input_names, alleles, scored_epitopes, scored_peptides)
    
    with open(args.html_report, 'w') as f:
        f.write(html)

