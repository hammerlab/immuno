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
import datetime

import pandas as pd
from Bio import SeqIO
import numpy as np

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

def build_html_report(scored_data):
    scored_data = scored_data.sort(columns=('combined_score',), ascending=False)
    table = scored_data.to_html(
        index=False,
        na_rep="-",
        columns = [
            'SourceSequence', 'info', 'ref', 'alt', 'pos',
            'Epitope', 'EpitopeStart', 'EpitopeEnd',
            'percentile_rank', 'ann_rank', 'immunogenicity',
            'combined_score'
        ])

    # take each source sequence and shade its amino acid letters
    # based on the average score of each epitope containing that letter
    seq_divs = []
    seq_scores = []

    for seq in scored_data.SourceSequence.unique():
        scores = np.zeros(len(seq), dtype=float)
        imm_scores = np.zeros(len(seq), dtype=float)
        mhc_scores = np.zeros(len(seq), dtype=float)
        score_counts = np.ones(len(seq), dtype=int)
        rowslice = scored_data[scored_data.SourceSequence == seq]
        gene_info = None
        for _, row in rowslice.iterrows():
            gene_info = row['info']
            start = row['EpitopeStart'] - 1
            stop = row['EpitopeEnd']
            scores[start:stop] += row['combined_score']
            imm_scores[start:stop] += row['immunogenicity']
            mhc_scores[start:stop] += row['percentile_rank'] / 100.0
            score_counts[start:stop] += 1

        # default background for all letters of the sequence is gray
        # but make it more red as the score gets higher
        letters = []
        colors = []
        imm_colors = []
        mhc_colors = []
        scores /= score_counts
        imm_scores /= score_counts
        mhc_scores /= score_counts
        for i in xrange(len(seq)):
            letter = seq[i]
            score = scores[i]
            letter_td = "<td>%s</td>" % letter
            letters.append(letter_td)

            imm = imm_scores[i]
            mhc = mhc_scores[i]
            maxval = 256
            mhc_intensity = int(mhc*maxval)
            mhc_rgb = "rgb(%d, %d, %d)" % (mhc_intensity, mhc_intensity/2, 0)
            imm_intensity =  int(imm*maxval)
            imm_rgb = "rgb(%d, %d, %d)" % (0, imm_intensity/2, imm_intensity)


            color_cell = \
            """
            <td style="background-color: %s;">&nbsp;</td>
            """
            mhc_color_cell = color_cell %  mhc_rgb
            imm_color_cell = color_cell % imm_rgb
            imm_colors.append(imm_color_cell)
            mhc_colors.append(mhc_color_cell)

        median_score = np.median(scores)
        letters_cols = '\n\t'.join(letters)
        mhc_color_cols = '\n\t'.join(mhc_colors)
        imm_color_cols = '\n\t'.join(imm_colors)
        colored_letters_table = \
            """
            <center>
            <table border="1">
            <tr>
            <td style='background-color: rgb(190,190,190);'>Sequence</td>
            %s
            </tr>
            <tr>

            <td style='background-color: rgb(190,190,190);'>MHC Binding</td>
            %s
            </tr>
            <tr>

            <td style='background-color: rgb(190,190,190);'>Immunogenicity</td>
            %s
            </tr>
            </table>
            </center>
            """ % (letters_cols, mhc_color_cols, imm_color_cols)

        div = """
            <div style="border-bottom: 1px solid gray; margin-bottom: 1em;">
            <h3>Median Epitope Score = %0.4f (%s)</h3>
            %s
            <br>
            </div>
            """ % (median_score, gene_info, colored_letters_table)
        seq_divs.append(div)
        seq_scores.append(median_score)

    seq_order = reversed(np.argsort(seq_scores))
    seq_divs_html = "\n".join(seq_divs[i] for i in seq_order)

    page = """
        <html>
        <style>
            body { padding: 1em; }
            table, td, th
            {
                border:1px solid gray;
                text-align:center;
                padding: 0em;
            }
            td {
                height: 2em;
                width:1.5em;
                background-color:
                rgb(220,220,220);
            }
            th { background-color: rgb(90, 190, 240); }
        </style>
        <head><title>Immune Pipeline Results (%s)</title></head>
        <body>
        <h2>Mutation Regions</h2>
        %s
        <hr>
        <h2>Sorted Scores Results</h2>
        %s
        </body>
        </html>
    """ % (datetime.date.today(), seq_divs_html, table)
    with open('results.html', 'w') as f:
        f.write(page)


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

    # TODO: combine based on commandline args mhc, immunogenicity, etc..
    scored_data['combined_score']= \
        (
            scored_data['percentile_rank'] / 100.0 + \
            scored_data['immunogenicity']
        ) / 2.0
    scored_data = scored_data.sort(columns=('combined_score',))
    if args.output:
        scored_data.to_csv(args.output, index=False)
    else:
        print(scored_data.to_string())
    build_html_report(scored_data)

