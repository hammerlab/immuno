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

from load_maf import peptides_from_maf
from load_vcf import peptides_from_vcf
from load_snpeff import peptides_from_snpeff
from load_fasta import peptides_from_fasta

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
    parser.add_argument("--peptide_length",
        default=31,
        type = int,
        help="length of vaccine peptides (may contain multiple epitopes)")
    parser.add_argument("--allele_file",
        help="file with one allele per line")
    parser.add_argument("--alleles",
        help="comma separated list of allele (default HLA-A*02:01)")
    parser.add_argument("--output",
        help="output file for dataframes", required=False)
    parser.add_argument("--no-mhc",
        default=False,
        action="store_true",
        help="Don't predict MHC binding")

    args = parser.parse_args()

    # stack up the dataframes and later concatenate in case we
    # want both commandline strings (for weird mutations like translocations)
    # and files
    mutated_region_dfs = []

    if args.string:
        # allow multiple strings to be specified in comma-separated list
        starts = []
        stops = []
        full_peptides = []
        for string in args.string.split(","):
            full_peptide = string.upper().strip()
            # allow the user to specify mutated region of the amino acid
            # string QLSQ_Y_QQ (the full peptide is QLSQYQQ and Y is mutated)
            parts = full_peptide.split("_")
            if len(parts) == 1:
                full_peptide = parts[0]
                start = 0
                stop = len(full_peptide)
            elif len(parts) == 2:
                full_peptide = parts[0] + parts[1]
                start = len(parts[0])
                stop = len(full_peptide)
            else:
                assert len(parts) == 3, \
                    "Can't parse peptide string %s" % full_peptide
                full_peptide = parts[0] + parts[1] + parts[2]
                start = len(parts[0])
                stop = start + len(parts[1])
            full_peptides.append(full_peptide)
            starts = starts.append(start)
            stops.append(stop)
        df = pd.DataFrame({
            'SourceSequence': full_peptides,
            'MutationStart' : starts,
            'MutationEnd' : stops,
            'info' : ['commandline'] * len(full_peptides),
        })
        mutated_region_dfs.append(df)


    # loop over all the input files and
    # load each one into a dataframe

    for input_filename in args.input:
        if input_filename.endswith("eff.vcf"):
            df = peptides_from_snpeff(input_filename)
        if input_filename.endswith(".vcf"):
            df = peptides_from_vcf(input_filename)
        elif input_filename.endswith(".maf"):
            df = peptides_from_maf(input_filename)
        elif input_filename.endswith(".fasta") \
                or input_filename.endswith(".fa"):
            df = peptides_from_fasta(input_filename)
        else:
            assert False, "Unrecognized file type %s" % input_filename
        mutated_region_dfs.append(df)


    if len(mutated_region_dfs) == 0:
        parser.print_help()
        print("\nERROR: Must supply at least --string or --input")
        sys.exit()
    mutated_regions = pd.concat(mutated_region_dfs)

    peptide_length = int(args.peptide_length)
    # peptides = peptide_substrings(full_peptide, peptide_length)

    # get rid of gene descriptions if they're in the dataframe
    if args.allele_file:
        alleles = [l.strip() for l in open(args.allele_file)]
    elif args.alleles:
        alleles = [l.strip() for l in args.alleles.split(",")]
    else:
        alleles = [DEFAULT_ALLELE]

    mhc = IEDBMHCBinding(name = 'mhc', alleles=alleles)
    mhc_data = mhc.apply(mutated_regions)

    immunogenicity = ImmunogenicityRFModel(name = 'immunogenicity')
    scored_data = immunogenicity.apply(mhc_data)


    # IEDB returns a noisy column with 'method' descriptions, drop it
    if 'method' in scored_data.columns:
        scored_data = scored_data.drop('method', axis = 1)

    # strong binders are considered percentile_rank < 2.0
    mhc_percentile = scored_data['percentile_rank']
    # rescale [0,2] -> [0,1]
    mhc_score = mhc_percentile / 2
    # treat anything at percentile above 2 as just as bad as 2
    mhc_score[mhc_score > 1] = 1
    # lower percentiles are stronger binders
    mhc_score = 1.0 - mhc_score
    scored_data['mhc_score'] = mhc_score

    # rescale immune score so anything less than 0.5 is 0
    imm_score = scored_data['immunogenicity']
    imm_score -= 0.5
    imm_score *= 2
    imm_score[imm_score < 0] = 0
    scored_data['imm_score'] = imm_score

    scored_data['combined_score']= (mhc_score + imm_score) / 2.0

    scored_data = scored_data.sort(columns=('combined_score',))
    if args.output:
        scored_data.to_csv(args.output, index=False)
    else:
        print(scored_data.to_string())

    html = build_html_report(scored_data)
    with open('results.html', 'w') as f:
        f.write(html)

