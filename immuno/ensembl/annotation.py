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
from Bio.Seq import Seq

from annotation_data import EnsemblAnnotationData

data = EnsemblAnnotationData()

def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

def get_exons_from_transcript(transcript_id):
    """
    Filter exons down to those with this transcript_id

    Parameters
    ----------
    transcript_id : transcript id, of the from EST#####

    returns Pandas dataframe containing only those exons
    """
    assert transcript_id in data.transcript_exons_dict, (
            "Unknown transcript %s in a dict of size %s, "
            "with initial elements %s") % (
                    transcript_id,
                    len(data.transcript_exons_dict),
                    data.transcript_exons_dict.items()[:5])
    exons = data.transcript_exons_dict[transcript_id]
    assert len(exons) > 0, \
        "Couldn't find exons for transcript %s" % transcript_id
    fields = [
        'exon_id',
        'seq_start',
        'start_exon_id',
        'seq_end',
        'end_exon_id',
        'stable_id_exon',
        'seq_region_start_exon',
        'seq_region_end_exon'
    ]
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

def get_strand(transcript_id):
    """
    Gets the strand of the gene

    Parameters
    ----------
    transcript_id :
        Transcript id, of the from EST#####

    Return strand : int, +1 for forward strand else -1
    """
    try:
        exons = data.transcript_exons_dict[transcript_id]
    except KeyError:
        exons = []

    if len(exons) == 0:
        logging.warn("Transcript %s has no sequence information", transcript_id)
        return 1
    strand = exons['seq_region_strand_gene'].iget(0)
    return strand

def is_forward_strand(transcript_id):
    return get_strand(transcript_id) > 0

def is_incomplete_cds(transcript_id):
    """
    Compute 5 prime incomplete CDS - checks the start_phase of the first exon

    transcript_id :
        Transcript id, of the from EST#####
    """
    start_phase = get_cds_start_phase(transcript_id)
    if start_phase:
        return start_phase > 0
    else:
        return False

def get_start_exon(transcript_id):
    start_exons = data.start_exons_dataframe
    mask = start_exons['stable_id_transcript'] == transcript_id
    start_exon = start_exons[mask]
    if start_exon.empty:
        return None
    else:
        return start_exon.to_dict(outtype = 'records')[0]

def get_cds_start_phase(transcript_id):
    """
    Compute CDS start_phase - checks the phase of the first exon

    transcript_id :
        Transcript id, of the from EST#####
    """
    start_exon = get_start_exon(transcript_id)
    if start_exon is not None:
        return max(start_exon['phase'], 0)  # -1 means no phase (equiv. to zero)
    else:
        return None

def get_transcript_index_from_pos(
        pos,
        transcript_id,
        skip_untranslated_region= True):
    """
    Gets the index into to the transcript from genomic position
    The transcript is composed of spliced exons that have genomic start and
    stop positions.  This function searches for the exon that matches this
    genomic positions and returns the index into the transcript

    Parameters
    ----------
    position : int
        Genomic position in the contig

    transcript_id :
        Transcript id, of the from EST#####

    skip_untranslated_region : bool, optional
        If True (default), then give position in the CDS (coding sequence),
        otherwise give position in the longer full cDNA sequence.
    """
    exons = get_exons_from_transcript(transcript_id)
    exons = exons.sort(columns=['seq_region_start_exon', 'seq_region_end_exon'])
    exons['exon_length'] = \
        exons['seq_region_end_exon'] - exons['seq_region_start_exon'] + 1
    starts = exons['seq_region_start_exon']
    stops = exons['seq_region_end_exon']
    intervals = zip(starts, stops)

    transcript_length = exons['exon_length'].sum()
    transcript_idx = get_idx_from_interval(pos, intervals)

    if transcript_idx is None:
        logging.warning("Couldn't find position %d in transcript %s",
            pos, transcript_id)
    else:
        # Reverse array index if on reverse strand
        forward = is_forward_strand(transcript_id)
        transcript_idx = transcript_idx if forward else \
            transcript_length - transcript_idx - 1
        if skip_untranslated_region:
            # Adjust for translations (CDS) start region
            prefix_utr_length = get_five_prime_utr_length(exons, forward)
            if transcript_idx < prefix_utr_length:
                logging.warn(
                    "UTR mutation at cDNA position %d, transcript %s",
                    transcript_idx, transcript_id)
                return None
            else:
                transcript_idx -= prefix_utr_length

        # Adjust for CDS start phase if first exon is out of phase
        transcript_phase = get_cds_start_phase(transcript_id)
        transcript_idx += transcript_phase
        if transcript_phase > 0:
            logging.warn("Transcript %s is incomplete", transcript_id)

        # TODO: check that index is within the mRNA transcript
        # need to get the length of the coding region from the transcript_id
        #suffix_utr_length = get_three_prime_utr_length(exons, forward)
        #assert transcript_idx <= transcript_length + suffix_utr_length

    return transcript_idx

def get_five_prime_utr_length(exons_df, forward = True):
    """
    Gets the length of the 5' UTR from a set of sorted exons
    from a specifc transcript

    Parameters
    ----------
    exons : Pandas dataframe with 'exon_id', 'seq_region_end_exon',
            'seq_region_start_exon'
            Also, 'start_exon_id' marks the exon that starts translation and
            'seq_start' is the offset into the first translated exon

    forward : bool, default = True, is forward strand or not

    Return utr_length : int
    """
    exons_df = exons_df.sort(
        columns=['seq_region_start_exon', 'seq_region_end_exon'],
        ascending=[forward, forward])
    utr_length = 0
    for idx, row in exons_df.iterrows():
        if row['exon_id'] == row['start_exon_id']:
            utr_length += row['seq_start'] - 1
            return utr_length
        else:
            utr_length += \
                row['seq_region_end_exon'] - row['seq_region_start_exon'] + 1
    return None

def get_three_prime_utr_length(exons_df, forward = True):
    """
    Gets the length of the 3' UTR from a set of sorted exons from
    a specifc transcript

    Parameters
    ----------
    exons : Pandas dataframe with columns:
              'exon_id', 'seq_region_end_exon', 'seq_region_start_exon'
            Also, 'end_exon_id' marks the exon that starts translation and
            'seq_end' is the offset into the last translated exon

    forward : bool, default = True, is forward strand or not

    Return utr_length : int
    """
    reverse = not forward
    exons_df = exons_df.sort(
        columns=['seq_region_start_exon', 'seq_region_end_exon'],
        ascending=[reverse, reverse])
    utr_length = 0
    for idx, row in exons_df.iterrows():
        exon_length = \
            row['seq_region_end_exon'] - row['seq_region_start_exon'] + 1
        if row['exon_id'] == row['end_exon_id']:
            utr_length += exon_length - row['seq_end']
            return utr_length
        else:
            utr_length += exon_length
    return None

def annotate_vcf_transcripts(vcf_df):
    """
    Expand each variant in a DataFrame into multiple entries for
    all transcript_ids that could contain the mutated position
    Parameters
    ----------
    vcf_df : Pandas DataFrame with chr, pos, ref, alt columns

    Return DataFrame with extra columns:
        - 'name'
        - 'stable_id_gene'
        - 'description_gene'
        - 'seq_region_start_gene'
        - 'seq_region_end_gene'
        - 'stable_id_transcript'
        - 'seq_region_start_transcript'
        - 'seq_region_end_transcript'
    """

    genome_annotation_data_df = data.transcripts_dataframe

    # combine each variant entry from `vcf_df` with all possible gene
    # annotations  on the same chromosome from `genome_annotation_data_df`
    all_possible_transcripts = vcf_df.merge(
        genome_annotation_data_df, left_on='chr', right_on='name', how='left')
    variant_position = all_possible_transcripts['pos']
    transcript_start = all_possible_transcripts['seq_region_start_transcript']
    transcript_stop = all_possible_transcripts['seq_region_end_transcript']
    mask = (variant_position > transcript_start) & \
        (variant_position < transcript_stop)
    annotated = all_possible_transcripts[mask]
    return annotated.drop_duplicates()
