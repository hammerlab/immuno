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

from common import peptide_substrings
from pipeline import ImmunoPipeline
from immunogenicity import ImmunogenicityRFModel
from binding import IEDBMHCBinding

from load_maf import peptides_from_maf
from load_vcf import peptides_from_vcf
from load_snpeff import peptides_from_snpeff

def peptides_from_fasta(fasta_files, peptide_length):
    epitope_dataframes = []
    peptides = []
    source_seqs = []
    filenames = []
    for fasta_file in fasta_files:
        fasta_data = SeqIO.parse(fasta_file, 'fasta')
        seqs = list(set([e.seq for e in fasta_data]))
        for seq in seqs:
            curr_peptides = peptide_substrings(seq, peptide_length)
            peptides.extend(curr_peptides)
            source_seqs.extend([seq] * len(curr_peptides))
            filenames.extend([fasta_file] * len(curr_peptides))
    assert len(peptides) == len(source_seqs) == len(filenames)
    return pd.DataFrame({
        'Peptide': peptides,
        'SourceSequence': source_seqs,
        'Filename': filenames,
    })


DEFAULT_ALLELE = 'HLA-A*02:01'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # must supply either an input file or an amino acid string
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", action="append", default=[],
        help="input file name (i.e. FASTA, MAF, VCF)")
    input_group.add_argument("--string", default = None,
        help="amino acid string")
    parser.add_argument("--peptide_length", default=31,
        help="length of vaccine peptides (may contain multiple epitopes)")
    parser.add_argument("--allele_file",
        help="file with one allele per line")
    parser.add_argument("--alleles",
        help="comma separated list of allele (default HLA-A*02:01)")
    parser.add_argument("--output",
        help="output file for dataframes", required=False)
    parser.add_argument("--immunogenicity",
        default=True, help="Predict immunogenicity score")
    parser.add_argument("--mhc", default=True,
        help="Predict MHC binding")

    args = parser.parse_args()
    if args.string:
        full_peptide = args.string.upper().strip()
        n = len(full_peptide)
        peptide_length = min(int(args.peptide_length), n)
        peptides = peptide_substrings(full_peptide, peptide_length)
        epitope_data = pd.DataFrame({
            'Peptide': peptides,
            'SourceSequence': [full_peptide] * len(peptides)
        })

    # TODO: allow for multiple input files by concatenating
    # their dataframes
    elif len(args.input) > 0:
        input_filename = args.input[0]
        peptide_length = int(args.peptide_length)
        if input_filename.endswith("eff.vcf"):
            epitope_data = peptides_from_snpeff(
                input_filename, peptide_length)
        if input_filename.endswith(".vcf"):
            epitope_data = peptides_from_vcf(
                input_filename, peptide_length)
        elif input_filename.endswith(".maf"):
            epitope_data = peptides_from_maf(args.input, peptide_length)
        elif input_filename.endswith(".fasta") \
                or input_filename.endswith(".fa"):
            epitope_data = peptides_from_fasta(
                args.input, peptide_length)
        else:
            assert False, "Unrecognized file type %s" % input_filename
    else:
        assert False, \
            "Either amino acid string or input file required"

    # get rid of gene descriptions if they're in the dataframe
    if 'description_gene' in epitope_data.columns:
        epitope_data = epitope_data.drop('description_gene', axis = 1)

    if args.allele_file:
        alleles = [l.strip() for l in open(args.allele_file)]
    elif args.alleles:
        alleles = [l.strip() for l in args.alleles.split(",")]
    else:
        alleles = [DEFAULT_ALLELE]

    pipeline = ImmunoPipeline()

    if args.mhc not in (False, 0, "False", "0"):
        mhc = IEDBMHCBinding(name = 'mhc', alleles=alleles)
        pipeline.add_scorer(mhc)

    if args.immunogenicity not in (False, 0, "False", "0"):
        immunogenicity = ImmunogenicityRFModel(name = 'immunogenicity')
        pipeline.add_scorer(immunogenicity)

    scored_data = pipeline.score(epitope_data)

    # some of the MHC scores come back as all NaN so drop them
    scored_data = scored_data.dropna(axis=1, how='all')

    # IEDB returns a noisy column with 'method' descriptions, drop it
    if 'method' in scored_data.columns:
        scored_data = scored_data.drop('method', axis = 1)

    # TODO: combine based on commandline args mhc, immunogenicity, etc..
    scored_data['combined_score']= \
        scored_data['percentile_rank'] / 100.0 + \
        scored_data['immunogenicity']
    scored_data = scored_data.sort(columns=('combined_score',))
    if args.output:
        scored_data.to_csv(args.output, index=False)
    else:
        print(scored_data.to_string())
