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

import pandas as pd
from Bio import SeqIO

from pipeline import ImmunoPipeline
from immunogenicity import ImmunogenicityRFModel
from binding import IEDBMHCBinding
from maf_to_epitopes import get_eptiopes_from_maf
from epitope_generation import generate_epitopes_from_snpeff, generate_epitopes_from_vcf

def get_epitopes_from_fasta(fasta_files):
    epitope_dataframes = []
    for fasta_file in fasta_files:
        epitope_data = SeqIO.parse(fasta_file, 'fasta')
        peptides = list(set([e.seq for e in epitope_data]))
        df = pd.DataFrame({'peptide': peptides})
        epitope_dataframes.append(df)
    return pd.concat(epitope_dataframes)

def add_scoring(pipeline, alleles):
    mhc = IEDBMHCBinding(name = 'mhc', alleles=alleles)
    pipeline.add_scorer(mhc)

    immunogenicity = ImmunogenicityRFModel(name = 'immunogenicity')
    pipeline.add_scorer(immunogenicity)

    return pipeline

DEFAULT_ALLELE = 'HLA-A*02:01'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # must supply either an input file or an amino acid string
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", action="append", default=[],
        help="input file name (i.e. FASTA, MAF, VCF)")
    input_group.add_argument("--string", action="append", default=[],
        help="amino acid string")
    parser.add_argument("--allele_file",
        help="file with one allele per line")
    parser.add_argument("--alleles",
        help="comma separated list of allele")
    parser.add_argument("--output",
        help="output file for dataframes", required=False)

    args = parser.parse_args()

    if args.string:
        epitope_data = pd.DataFrame.from_dict({'Epitope' : args.string})
    elif len(args.input) > 0:
        input_filename = args.input[0]
        if input_filename.endswith("eff.vcf"):
            epitope_data = generate_epitopes_from_snpeff(input_filename)
        if input_filename.endswith(".vcf"):
            epitope_data = generate_epitopes_from_vcf(input_filename)
        elif input_filename.endswith(".maf"):
            epitope_data = get_eptiopes_from_maf(args.input)
        elif input_filename.endswith(".fasta") or input_filename.endswith(".fa"):
            epitope_data = get_epitopes_from_fasta(args.input)
        else:
            assert False, "Unrecognized file type %s" % input_filename
    else:
        assert False, \
            "Either amino acid string or input file required"

    if args.allele_file:
        alleles = [l.strip() for l in open(args.allele_file)]
    elif args.alleles:
        alleles = [l.strip() for l in args.alleles.split(",")]
    else:
        alleles = [DEFAULT_ALLELE]

    pipeline = ImmunoPipeline()
    add_scoring(pipeline, alleles)
    scored_data = pipeline.score(epitope_data)


    if args.output:
        scored_data.to_csv(args.output, index=False)
    else:
        print(scored_data)
