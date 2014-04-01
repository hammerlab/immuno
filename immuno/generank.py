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
Prints ranks for immune-related genes in the given sample.
"""


import csv
import argparse

import pandas as pd
import numpy as np
import entrez
import rnaseq_atlas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample",
        help="Sample RNAseq file, one column for  gene IDs and another of RSEM values",
        required=True)
    parser.add_argument(
        "--sample-uses-hugo-ids",
        required=False,
        default=False,
        action="store_true",
        help="Don't convert sample IDs from Entrez, they're already Hugo names")
    parser.add_argument(
        "--normal",
        default = None,
        required=False,
        help="Reference values for normal tissues mapping Entrez to RSEM")
    parser.add_argument(
        "--tissue",
        help="Which column to select from the reference")
    parser.add_argument(
        "--percentile",
        help="Which rank percentile do we consider 'low'?",
        type = float,
        default = 0.25)

    args = parser.parse_args()

    if not args.sample_uses_hugo_ids:
        entrez_df = pd.read_csv(args.sample, sep='\t', names=('Entrez', 'RSEM'))

        # mapping from Entrez IDs to Hugo and RefSeq
        hugo_mapping = entrez.entrez_hugo_dataframe()
        # add Hugo column
        merged_df = entrez_df.merge(hugo_mapping, on='Entrez')
        # keep only the Hugo IDs and RSEM
        hugo_df = merged_df[['Hugo', 'RSEM']]
    else:
        hugo_df = pd.read_csv(args.sample, sep='\t', names=('Hugo', 'RSEM'))
    hugo_df = hugo_df.groupby("Hugo").mean().reset_index()

    group_rank_method = 'min'

    if args.normal is None:
        normal_df = rnaseq_atlas.hugo_to_rank(group_rank_method)
        tissue_cols = rnaseq_atlas.TISSUE_COLUMNS
    else:
        normal_df = pd.read_csv(args.normal)
        normal_df = hugo_mapping.merge(normal_df, on = "Entrez")
        normal_df = normal_df.drop("Entrez", axis=1)
        normal_df = normal_df.drop("RefSeq", axis=1)
        tissue_cols = normal_df.columns[1:]
        normal_df[tissue_cols] = \
            normal_df[tissue_cols].rank(method = group_rank_method) - 1
        normal_df[tissue_cols] /= len(normal_df)

    # result has a Hugo ID column, RSEM values,
    # and all the tissue-specific ranks from normal_df
    combined = hugo_df.merge(normal_df).set_index("Hugo")
    rsem_ranks = combined.pop('RSEM').rank(method = group_rank_method) - 1

    rsem_ranks /= len(rsem_ranks)
    n = len(combined)
    counts = pd.DataFrame({
        "Hugo": combined.index,
        "Up" : np.zeros(n, dtype=float),
        "Down" : np.zeros(n, dtype=float),
        "PrctDiff" :  np.zeros(n, dtype=float),
    }).set_index("Hugo")

    p = args.percentile
    if p > 1.0:
        p /= 100.0
    for tissue in tissue_cols:

        print tissue

        tissue_ranks = combined[tissue]

        normal_low = tissue_ranks < p
        normal_high = tissue_ranks > (1.0 - p)

        cancer_low = rsem_ranks < p
        cancer_high = rsem_ranks > (1.0 - p)

        decrease = normal_high & cancer_low
        increase = normal_low & cancer_high
        counts.Up += increase
        counts.Down += decrease
        counts.PrctDiff += (rsem_ranks - tissue_ranks)

    counts['CountDiff'] = counts.Up - counts.Down
    counts.PrctDiff /= len(tissue_cols)
    k = 30
    print "Top %d over-expressed genes" % k
    print counts.sort("PrctDiff", ascending=False)[:k]

    def tally(category, filename):
        with open(filename, 'r') as f:
            genes = f.read().splitlines()
        print
        print category
        score = 0
        for gene in genes:
            row = counts.ix[counts.index == gene]
            diff = float(row.PrctDiff)
            print "  - %s : %s" % (gene, diff)
            score += diff
        print "Mean: %s" % (score / len(genes))
        print
        print "---"
        print

    tally("HLA Type I", 'hla_type1_genes.txt')
    tally("HLA Type II", 'hla_type2_genes.txt')
    tally("Immune Genes", "immune_genes.txt")
