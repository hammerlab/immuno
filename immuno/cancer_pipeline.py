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

from __future__ import print_function
import argparse
import sys


import pandas as pd
from Bio import SeqIO
import numpy as np

from common import peptide_substrings
from immunogenicity import ImmunogenicityRFModel
from binding import IEDBMHCBinding
from load_file import load_file
from strings import load_comma_string
from vaccine_peptides import build_peptides_dataframe
from report import build_html_report


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
    parser.add_argument("--allele-file",
        help="file with one allele per line")
    parser.add_argument(
        "--alleles",
        help="comma separated list of allele (default HLA-A*02:01)")
    parser.add_argument(
        "--output",
        help="output file for dataframe containing scored epitopes",
        required=False)
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
        print("\nERROR: Must supply at least --string or --input")
        sys.exit()
    mutated_regions = pd.concat(mutated_region_dfs)

    # get rid of gene descriptions if they're in the dataframe
    if args.allele_file:
        alleles = [l.strip() for l in open(args.allele_file)]
    elif args.alleles:
        alleles = [l.strip() for l in args.alleles.split(",")]
    else:
        alleles = [DEFAULT_ALLELE]



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
        # strong binders are considered percentile_rank < 2.0
        mhc_percentile = scored_epitopes['percentile_rank']
        mhc_score = (100.0 - mhc_percentile) / 100.0
        scored_epitopes['mhc_score'] = mhc_score

    immunogenicity = ImmunogenicityRFModel(name = 'immunogenicity')
    scored_epitopes = immunogenicity.apply(scored_epitopes)
    imm_score = scored_epitopes['immunogenicity']

    # TODO: make the imm score based on percentile in normal human proteins
    scored_epitopes['imm_score'] = imm_score


    if args.skip_mhc:
        combined_score = imm_score
    else:
        combined_score = (mhc_score + imm_score) / 2.0

    scored_epitopes['combined_score'] = combined_score
    scored_epitopes = scored_epitopes.sort(columns=('combined_score',))
    if args.output:
        scored_epitopes.to_csv(args.output, index=False)
    else:
        print(scored_epitopes.to_string())

    scored_peptides = build_peptides_dataframe(scored_epitopes,
        peptide_length = peptide_length)

    html = build_html_report(scored_epitopes, scored_peptides)
    with open(args.html_report, 'w') as f:
        f.write(html)

