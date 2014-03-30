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

from copy import deepcopy

import pandas as pd
import numpy as np

import epitopes.mutate as mutate

from common import peptide_substrings
from ensembl import annotation
from ensembl.transcript_variant import peptide_from_transcript_variant


def _shorten_chromosome_name(chr):
    if chr.startswith('chr'):
        chr = chr[-1]
        if chr =='M':
            return 'MT'
        return chr
    return chr

def vcf_to_dataframe(vcf_filename):
    """
    Transforms a VCF file to a Pandas Dataframe

    Parameters
    ----------
    vcf_filename : Path to VCF file

    Returns Dataframe with fields
            - 'chr'
            - 'pos'
            - 'id'
            - 'ref'
            - 'alt'
            - 'qual'
            - 'filter'
            - 'info'
    """
    with open(vcf_filename) as fd:
        lines_to_skip = 0
        while next(fd).startswith('#'):
            lines_to_skip += 1

    # first 8 columns of a VCF file are required to be:
    #   chr    : chromosome (i.e., '20', 'chr20', 'MT')
    #   pos    : where's the variant on the chromsome?
    #   id     : dbSnp identifier of variant, if available (i.e. 'rs11449')
    #   ref    : reference letter(s)
    #   alt    : alternate letters(s) of variant
    #   qual   : phred-scaled quality score
    #   filter : "PASS" if variant passes all filters
    #   info   : optional info, comma separated key=value list
    header = ['chr', 'pos', 'id', 'ref', 'alt','qual', 'filter', 'info']

    df = pd.read_csv(
            vcf_filename,
            sep='\t',
            skiprows=lines_to_skip,
            header = None,
            names=header,
            usecols=header,
            dtype={'pos' : np.int32, 'chr':str})

    df['chr'] = df.chr.map(_shorten_chromosome_name)
    return df

def peptides_from_vcf(input_file, length=31, log_filename = 'run.log'):
    vcf_df = vcf_to_dataframe(input_file)
    transcripts_df = annotation.annotate_transcripts(vcf_df)

    def peptides_from_annotation(group):
        row = group.irow(0)

        transcript_id = row['stable_id_transcript']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        rows = []
        if transcript_id:
            full_peptide = \
                peptide_from_transcript_variant(
                    transcript_id, pos, ref, alt, min_padding = length)

        if full_peptide:
            peptides = peptide_substrings(full_peptide, length)
            for peptide in peptides:
                row = deepcopy(row)
                row['Peptide'] = peptide
                rows.append(row)
        new_df = pd.DataFrame.from_records(rows)
        return new_df
    cols = ['chr','pos', 'ref', 'alt']
    variants = transcripts_df.groupby(cols, group_keys=False)

    peptides = variants.apply(peptides_from_annotation)
    transcripts_df = transcripts_df.merge(peptides)
    if log_filename:
        transcripts_df.to_csv(log_filename, index=False)
    return transcripts_df
