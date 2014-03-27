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

import pandas as pd

from ensembl_io import download_transcript_metadata

class EnsemblData(object):
    def __init__(self):
        pass

    @property
    def transcript_metadata_path(self):
        if not hasattr(self, '_transcript_metadata_path'):
            self._transcript_metadata_path = download_transcript_metadata()
        return self._transcript_metadata_path

    @property
    def exon_data(self):
        """
        Dataframe containing exon data
        """
        if not hasattr(self, '_exon_data'):
            path = self.transcript_metadata_path
            self._exon_data = pd.read_csv(path, sep='\t')
        return self._exon_data

    @property
    def transcript_data(self):
        """
        Subset columns for transcript data only
        """
        if not hasattr(self, '_transcript_data'):
            transcript_cols = [
                'name', 'stable_id_gene', 'description_gene',
                'seq_region_start_gene', 'seq_region_end_gene',
                'stable_id_transcript', 'seq_region_start_transcript',
                'seq_region_end_transcript'
            ]
            self._transcript_data = self.exon_data[transcript_cols]
        return self._transcript_data

    @property
    def gene_data(self):
        """
        Subset columns for gene data only
        """
        if not hasattr(self, '_gene_data'):
            gene_cols = [
                'name', 'stable_id_gene', 'description_gene',
                'seq_region_start_gene', 'seq_region_end_gene'
            ]
            self._gene_data = self.transcript_data[gene_cols].drop_duplicates()

def _transcript_matches(transcript_row):
    position = transcript_row['pos']
    contig = transcript_row['chr']
    return transcript_row['seq_region_start_transcript'] < position \
            and transcript_row['seq_region_end_transcript'] > position \
            and transcript_row['name'] == contig:

def _gene_matches(gene_row):
    position = gene_row['pos']
    contig = gene_row['chr']
    return gene_row['seq_region_start_gene'] <= position \
            and gene_row['seq_region_end_gene'] > position \
            and gene_row['name'] == contig:

def annotate_transcripts(vcf_df):
    """
    Get list of transcript id from position

    Parameters
    ----------
    vcf : Pandas dataframe with chr, pos, ref, alt columns

    Return df with gene and transcript ids

    """
    if 'gene_stable_id' in vcf_df.columns:
        annotated = annotate(vcf_df,
            TRANSCRIPT_DATA, _transcript_matches, left_col=['chr', 'gene_stable_id'], right=['name', 'gene_stable_id'])
    else:
        annotated = annotate(vcf_df, TRANSCRIPT_DATA, _transcript_matches)

    return annotated

def annotate_genes(vcf_df):
    """
    Get list of gene id from position

    Parameters
    ----------
    vcf : Pandas dataframe with chr, pos, ref, alt columns

    Return df with gene ids

    """
    return annotate(vcf_df, GENE_DATA, _gene_matches)

def annotate(
        vcf_df,
        annotation_df,
        predicate,
        left_col = 'chr',
        right_col = 'name'):
    crossed = vcf_df.merge(
        annotation_df, left_on=left_col, right_on=right_col, how='left')
    annotated = crossed[crossed.apply(predicate, axis=1)]
    return annotated

def get_exons_from_transcript(transcript_id):
    """
    Filter exons down to those with this transcript_id

    Parameters
    ----------
    transcript_id : transcript id, of the from EST#####

    returns Pandas dataframe containing only those exons
    """
    exons = EXON_DATA[EXON_DATA['stable_id_transcript'] == transcript_id]
    fields = ['stable_id_exon', 'seq_region_start_exon', 'seq_region_end_exon']
    return exons[fields]

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

    returns
    """
    exons = get_exons_from_transcript(transcript_id)
    exons = exons.sort(columns=['seq_region_start_exon', 'seq_region_end_exon'])
    return get_idx_from_interval(
        pos, zip(exons['seq_region_start_exon'], exons['seq_region_end_exon']))

def get_idx_from_interval(pos, intervals):
    idx = 0
    for (start, end) in intervals:
        if pos > end:
            idx += (end - start) + 1
        elif pos <= end and pos >= start:
            return idx + (pos - start)
        else:
            ## error some
            return None

