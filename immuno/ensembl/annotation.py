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

import logging

import pandas as pd
from annotation_data import EnsemblAnnotationData

data = EnsemblAnnotationData()



def get_exons_from_transcript(transcript_id):
    """
    Filter exons down to those with this transcript_id

    Parameters
    ----------
    transcript_id : transcript id, of the from EST#####

    returns Pandas dataframe containing only those exons
    """
    exon_data = data.exon_data
    exons = exon_data[exon_data['stable_id_transcript'] == transcript_id]
    assert len(exons) > 0, \
        "Couldn't find exons for transcript %s" % transcript_id
    fields = ['stable_id_exon', 'seq_region_start_exon', 'seq_region_end_exon']
    return exons[fields]


def get_idx_from_interval(pos, intervals):
    idx = 0
    for (start, end) in intervals:
        if pos > end:
            idx += (end - start) + 1
        elif pos <= end and pos >= start:
            return idx + (pos - start)
        else:
            return None


def get_transcript_index_from_pos(pos, transcript_id):
    """
    Gets the index into to the transcript from genomic position
    The transcript is composed of spliced exons that have genomic start and
    stop positions.  This function searches for the exon that matches this
    genomic positions and returns the index into the transcript

    Parameters
    ----------
    position : int, genomic position in the contig
    transcript_id : transcript id, of the from EST#####

    """
    exons = get_exons_from_transcript(transcript_id)
    exons = exons.sort(columns=['seq_region_start_exon', 'seq_region_end_exon'])
    starts = exons['seq_region_start_exon']
    stops = exons['seq_region_end_exon']
    intervals = zip(starts, stops)
    result = get_idx_from_interval(pos, intervals)
    if result is None:
        logging.warning("Couldn't find position %d in transcript %s",
            pos, transcript_id)
    return result


def annotate(
        vcf_df,
        annotation_df,
        predicate,
        left_col = 'chr',
        right_col = 'name'):
    crossed = vcf_df.merge(
        annotation_df, left_on=left_col, right_on=right_col, how='left')
    annotated = crossed[crossed.apply(predicate, axis=1)]
    return annotated.drop_duplicates()


def _transcript_matches(transcript_row):
    position = transcript_row['pos']
    contig = transcript_row['chr']
    return transcript_row['seq_region_start_transcript'] < position \
            and transcript_row['seq_region_end_transcript'] > position \
            and transcript_row['name'] == contig

def annotate_vcf_transcripts(vcf_df):
    """
    Get list of transcript id from position

    Parameters
    ----------
    vcf : Pandas dataframe with chr, pos, ref, alt columns

    Return df with extra columns:
        'name', 'stable_id_gene', 'description_gene',
        'seq_region_start_gene', 'seq_region_end_gene',
        'stable_id_transcript', 'seq_region_start_transcript',
        'seq_region_end_transcript'

    """
    if 'gene_stable_id' in vcf_df.columns:
        annotated = annotate(vcf_df,
            data.transcript_data,
            _transcript_matches,
            left_col=['chr', 'gene_stable_id'],
            right=['name', 'gene_stable_id'])
    else:
        annotated = annotate(vcf_df, data.transcript_data, _transcript_matches)
    return annotated



def _gene_matches(gene_row):
    position = gene_row['pos']
    contig = gene_row['chr']
    return gene_row['seq_region_start_gene'] <= position \
            and gene_row['seq_region_end_gene'] > position \
            and gene_row['name'] == contig


def annotate_vcf_genes(vcf_df):
    """
    Get list of gene id from position

    Parameters
    ----------
    vcf_df : Pandas dataframe
        Expected to have columns 'chr', 'pos', 'ref', 'alt'

    Return df with extra columns:
        'name', 'stable_id_gene', 'description_gene',
        'seq_region_start_gene', 'seq_region_end_gene'

    """
    return annotate(vcf_df, data.gene_data, _gene_matches)

