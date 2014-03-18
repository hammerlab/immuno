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


TRANSCRIPT_META_DATA_FILE = download_transcript_metadata()
EXON_DATA = pd.read_csv(TRANSCRIPT_META_DATA_FILE, sep='\t')

# Subset columns for transcript data only
TRANSCRIPT_DATA = EXON_DATA[['name', 'stable_id_gene', 'description_gene', 'seq_region_start_gene', 'seq_region_end_gene', 'stable_id_transcript', 
                'seq_region_start_transcript', 'seq_region_end_transcript']].drop_duplicates()

# Subset columns for gene data only
GENE_DATA = TRANSCRIPT_DATA[['name', 'stable_id_gene', 'description_gene', 'seq_region_start_gene', 'seq_region_end_gene']].drop_duplicates()

def _transcript_matches(transcript_row):
    position = transcript_row['pos']
    contig = transcript_row['chr']
    if transcript_row['seq_region_start_transcript'] < position and transcript_row['seq_region_end_transcript'] > position \
        and transcript_row['name'] == contig:
        return True
    else:
        return False

def _gene_matches(gene_row):
    position = gene_row['pos']
    contig = gene_row['chr']
    if gene_row['seq_region_start_gene'] <= position and gene_row['seq_region_end_gene'] > position \
        and gene_row['name'] == contig:
        return True
    else:
        return False

def annotate_transcripts(vcf_df):
    """ get list of transcript id from position

        Parameters
        ----------
        vcf : Pandas dataframe with chr, pos, ref, alt columns

        Return df with gene and transcript ids

    """
    if 'gene_stable_id' in vcf_df.columns:
        annotated = annotate(vcf_df, TRANSCRIPT_DATA, _transcript_matches, left_col=['chr', 'gene_stable_id'], right=['name', 'gene_stable_id'])
    else:
        annotated = annotate(vcf_df, TRANSCRIPT_DATA, _transcript_matches)

    return annotated

def annotate_genes(vcf_df):
    """ get list of gene id from position

    Parameters
    ----------
    vcf : Pandas dataframe with chr, pos, ref, alt columns

    Return df with gene ids

    """
    return annotate(vcf_df, GENE_DATA, _gene_matches)

def annotate(vcf_df, annotation_df, predicate, left_col = 'chr', right_col = 'name'):
    crossed = vcf_df.merge(annotation_df, left_on=left_col, right_on=right_col, how='left')
    annotated = crossed[crossed.apply(predicate, axis=1)]

    return annotated

def get_exons_from_transcript(transcript_id):
    """ filter exons down to those with this transcript_id

    Parameters
    ----------
    transcript_id : transcript id, of the from EST#####

    returns Pandas dataframe containing only those exons
    """
    exons = EXON_DATA[EXON_DATA['stable_id_transcript'] == transcript_id]
    return exons[['stable_id_exon', 'seq_region_start_exon', 'seq_region_end_exon']]

def get_transcript_index_from_pos(pos, transcript_id):
    """ gets the index into to the transcript from genomic position
    The transcript is composed of spliced exons that have genomic start and
    stop positions.  This function searches for the exon that matches this
    genomic positions and returns the index into the transcript

    Parameters
    ----------
    position : int, genomic position in the contig
    transcript_id : transcript id, of the from EST#####

    returns 
    """
    exons = get_exons_from_transcript(transcript_id).sort(columns=['seq_region_start_exon', 'seq_region_end_exon'])
    transcript_idx = 0
    for (idx, row) in exons.iterrows():
        if pos > row['seq_region_end_exon']:
            transcript_idx += row['seq_region_end_exon'] - row['seq_region_start_exon']
        elif pos < row['seq_region_end_exon'] and pos >= row['seq_region_start_exon']:
            return transcript_idx + (pos - row['seq_region_start_exon'])
        else:
            ## error some
            return None
