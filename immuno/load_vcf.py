from copy import deepcopy

import pandas as pd
import numpy as np

import epitopes.mutate as mutate

from common import peptide_substrings
import ensembl_annotation

_ensembl = EnsemblReferenceData()

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
    header = ['chr', 'pos', 'id', 'ref', 'alt','qual', 'filter', 'info']
    df = pd.read_csv(
            vcf_filename,
            sep='\t',
            skiprows=lines_to_skip,
            names=header,
            usecols=header,
            dtype={'pos' : np.int32})
    df['chr'] = df.chr.map(_shorten_chromosome_name)
    return df

def peptides_from_vcf(input_file, length=31):
    vcf_df = _vcf_to_dataframe(input_file)
    transcripts_df = ensembl_annotation.annotate_transcripts(vcf_df)

    def peptides_from_annotation(group):
        row = group.irow(0)
        transcript_id = row['stable_id_transcript']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        rows = []
        if transcript_id:
            full_peptide = \
                peptide_from_transcript(
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
    transcripts_df.to_csv('run.log', index=False)
    return transcripts_df
