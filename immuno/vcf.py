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
import logging

import pandas as pd
import numpy as np

from ensembl import annotation
from ensembl.transcript_variant import peptide_from_transcript_variant


def _shorten_chromosome_name(chr):
    chr = chr.replace('chr', '')
    if chr =='M':
        return 'MT'
    else:
        return chr

def parse_vcf(vcf_filename):
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

def load_vcf(
        input_filename,
        peptide_length=31,
        drop_low_quality = True,
        log_filename = 'vcf_csv.log'):
    """
    Parameters
    --------

    input_filename : str
        Path to VCF file

    peptide_length : int
        How long will the vaccine peptides be? Used to determine
        required padding.

    drop_low_quality : bool, optional
        Keep variants whose 'QUAL' columns doesn't say 'PASS' or '.'

    log_filename : str, optional
        Where should we log the parsed VCF dataframe?

    Returns a dataframe with columns:
        - chr : chomosome
        - pos : position in the chromosome
        - ref : reference DNA
        - alt : alternate DNA
        - info : gene name and entrez gene ID
        - stable_id_transcript : Ensembl transcript ID
        - SourceSequence : region of protein around mutation
        - MutationStart : first amino acid modified
        - MutationEnd : last mutated amino acid
        - MutationInfo : annotation i.e. V600E

    """
    vcf_df = parse_vcf(input_filename)

    # drop variants marked as low quality
    if drop_low_quality:
        qual = vcf_df['qual']
        mask = (qual == 'PASS') | (qual == '.')
        vcf_df = vcf_df[mask]

    logging.info("Loaded VCF %s with %d entries", input_filename, len(vcf_df))

    transcripts_df = annotation.annotate_vcf_transcripts(vcf_df)

    logging.info("Annotated VCF has %d entries", len(transcripts_df))

    new_rows = []
    group_cols = ['chr','pos', 'ref', 'alt', 'stable_id_transcript']
    for (_, pos, ref, alt, transcript_id), group in \
            transcripts_df.groupby(group_cols):
        row = group.irow(0)
        padding = peptide_length - 1 
        if transcript_id:
            logging.info("Getting peptide from transcript ID %s", transcript_id)
            seq, start, stop, annot = \
                peptide_from_transcript_variant(
                    transcript_id, pos, ref, alt,
                    padding = padding)
            assert isinstance(start, int), (start, type(start))
            assert isinstance(stop, int), (stop, type(stop))
        if seq:
            if '*' in seq:
                logging.warning(
                    "Found stop codon in peptide %s from transcript_id %s",
                    region.seq,
                    transcript_id)
            elif len(seq) <= padding:
                logging.info(
                    "Truncated peptide too short for transcript %s gene position %s %s > %s ", 
                    transcript_id, pos, ref, alt)
            else:
                row = deepcopy(row)
                row['SourceSequence'] = seq
                row['MutationStart'] = start
                row['MutationEnd'] = stop
                row['MutationInfo'] = annot
                new_rows.append(row)
    peptides = pd.DataFrame.from_records(new_rows)
    # peptides = variants.apply(peptides_from_annotation)
    transcripts_df = transcripts_df.merge(peptides)
    logging.info("Generated %d peptides from %s",
        len(transcripts_df), input_filename)

    # drop verbose or uninteresting columns from VCF
    for dumb_field in ('description_gene', 'filter', 'qual', 'id', 'name'):
        if dumb_field in transcripts_df.columns:
            transcripts_df = transcripts_df.drop(dumb_field, axis = 1)
    if log_filename:
        transcripts_df.to_csv(log_filename, index=False)
    return transcripts_df
