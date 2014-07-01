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
import numpy as np

def build_peptides_dataframe(
        epitopes_df,
        peptide_length,
        min_peptide_length = None,
        min_peptide_padding = 0):

    """
    Given a dataframe with a 'combined_score' for each short epitope
    create another dataframe of all longer vaccine peptides ranked by the
    median score of epitopes involving mutated amino acid residues minus
    epitopes which are unmodified by a mutation.

    Parameters
    ----------

    epitopes_df : pandas.DataFrame
        Must contain columns:
            - 'SourceSequence': window around the mutation in a protein
            - 'MutationStart': first mutated residue in the SourceSequence
            - 'MutationEnd' : last mutated reside in the SourceSequence
            - 'Epitope': short subset of the SourceSequence
            - 'EpitopeStart' : start position of epitope
            - 'EpitopeEnd' : end position of epitope
            - 'combined_score': score for the epitope
            - 'chr': chromosome of original DNA variant
            - 'pos': position  in the chromosome
            - 'ref': reference nucleotide(s) at `pos`
            - 'alt': alternate nucleotides found at `pos`
            - 'info': gene name and/or ID
            - 'stable_id_transcript': which transcript of the gene

    peptide_length : int
        How long should the subsets of SourceSequence be?

    min_peptide_length : int, optional
        If a SourceSequence is shorter than peptide_length, should we use it?
        Omitting min_peptide_length sets it equal to peptide_length.

    min_peptide_padding : int, optional 
        How many wildtype residues should be included before or after the mutation start? 

    Returns a new dataframe with columns:
        - 'Peptide':  amino acid sequence composed of multiple epitopes
        - 'SourceSequence' : longer amino acid sequence from which
            peptide was extracted
        - 'MutatedScore' : median score of epitopes with mutations
        - 'SelfScore' : median score of epitopes without mutations
        - 'Score' : log2(MutatedScore / SelfScore) if SelfScore > 0, else -inf
        - 'MutationStart' : first mutated reside in the peptide
        - 'MutationEnd' : last mutated residue in the peptide
        - 'chr': chromosome of original DNA variant
        - 'pos': position  in the chromosome
        - 'ref': reference nucleotide(s) at `pos`
        - 'alt': alternate nucleotides found at `pos`
        - 'info': gene name and/or ID
        - 'stable_id_transcript': which transcript of the gene
    """
    if min_peptide_length is None:
        min_peptide_length = peptide_length


    group_cols = [
        'SourceSequence',
        'MutationStart',
        'MutationEnd',
        'MutationInfo',
        'info',
        'stable_id_transcript'
    ]

    for required_field in group_cols:
        if required_field not in epitopes_df:
            epitopes_df[required_field] = '-'

    # if we loaded from VCF, these extra fields will be available:
    optional_cols = [
        'chr', 'pos', 'ref', 'alt',
    ]

    records = []
    # don't process a protein sequence / longer mutation region twice
    seen_sequences = set([])
    # don't process the same vaccine peptide twice
    seen_peptides = set([])
    for (seq, mut_start, mut_end, mut_info, info, transcript_id), group in\
            epitopes_df.groupby(group_cols):
        
        seq_cache_key = (seq, mut_start, mut_end)
        if seq_cache_key in seen_sequences:
            logging.info(
                "Already seen sequence %s (mut_start = %d, mut_end = %d), skipping", 
                seq,
                mut_start, 
                mut_end)
            continue
        else:
            logging.info(
                "Extracting vaccine peptides from sequence %s (mut_start = %d, mut_end=%d)",
                seq, 
                mut_start, 
                mut_end)
            seen_sequences.add(seq_cache_key)

        first_row = group.irow(0)

        # common properties that all the peptides we extract from
        # SourceSequence share
        base_record = {
            'SourceSequence' : seq,
            'MutationInfo' : mut_info,
            'TranscriptId' : transcript_id,
            'GeneInfo' : info,
        }
        for optional_field in optional_cols:
            if optional_field in first_row:
                base_record[optional_field] = first_row[optional_field]

        n = len(seq)
        if n >= peptide_length:
            window_size = peptide_length
        elif n >= min_peptide_length:
            window_size = min(n, min_peptide_length)
        else:
            logging.info("Skipping source sequence %s (length %d) from %s, shorter than %d",
                seq, len(seq), info, min_peptide_length)
            continue

        for peptide_start in xrange(min_peptide_padding, n  - window_size - min_peptide_padding + 1):

            peptide_end = peptide_start + window_size
            peptide = seq[peptide_start : peptide_end]
            
            # where is the mutation relative to this peptide?
            peptide_mut_start = mut_start - peptide_start
            peptide_mut_end = mut_end - peptide_start

            # if mutation isn't in the peptide, skip it
            if peptide_mut_start >= window_size or peptide_mut_end < 0:
                logging.info(
                    "Skipping self peptide %s from %s [%d:%d]",
                    peptide,
                    info,
                    peptide_start,
                    peptide_end)
                continue

            # if mutation goes outside the window, truncate it
            peptide_mut_start = max(peptide_mut_start, 0)
            peptide_mut_end = min(peptide_mut_end, window_size)

            # mutation start in the peptide should be between [0, len)
            # and end should be between [0, len]
            assert peptide_mut_start >= 0, "Expected non-negative start but got %s" % peptide_mut_start
            assert peptide_mut_start < window_size, "Expected peptide start <= %s but got %s" % (window_size, peptide_mut_start)
            assert peptide_mut_end >= 0, "Expected non-negative end but got %s" % peptide_mut_end
            assert peptide_mut_end <= window_size, "Expected peptide end <= %s but got %s" % (window_size, peptide_mut_end)
            
            peptide_cache_key = (peptide, peptide_mut_start, peptide_mut_end)
            if peptide_cache_key in seen_peptides:
                logging.info(
                    "Already added peptide %s mutation start = %d, mutation end = %d (from %s)", 
                    peptide, 
                    peptide_mut_start, 
                    peptide_mut_end, 
                    info)
                continue
            else:
                seen_peptides.add(peptide_cache_key)
            # copy fields of this record so we can add it to the data frame
            row = dict(base_record)
            row['Peptide'] = peptide
            row['PeptideStart'] = peptide_start
            row['PeptideEnd'] = peptide_end 
            row['PeptideMutationStart'] = peptide_mut_start
            row['PeptideMutationEnd'] = peptide_mut_end




            # To clarify the nomenclature (source seq vs. peptide vs. epitope)
            # look at this example.
            #
            # Let's say we have a small protein with this 48 residue sequence:
            #
            #   TAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQR
            #
            # and a somatic variant results in the a single amino acid change
            # of the first 'A' for a 'Q', then our mutated protein will be:
            #
            #   TQADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQR
            #
            # now either this full protein or some subset of it will become
            # the 'SourceSequence' depending on how large of a vaccine peptide
            # we're aiming to generate. If, for example, we only want 15-mer
            # vaccine peptides, then the SourceSequence will be the union
            # of all 15-mer windows containing the modified residue. In this
            # case, since the change was near the start of the protein, there
            # are only two such windows, yielding a SourceSequence of length 16:
            #
            #  TQADMAAQTTKHKWEA
            #
            # From this source sequence, we generate shorter substrings we
            # expect to bind to MHC molecules (typically 8-11 residues).
            # Assuming for simplicity that epitopes will all be 9-mers,
            # the set of epitopes will be:
            #
            #  TQADMAAQT
            #  QADMAAQTT
            #  ADMAAQTTK
            #  DMAAQTTKH
            #  MAAQTTKHK
            #  AAQTTKHKW
            #  AQTTKHKWE
            #
            # We look at the scores for each of these epitopes to generate
            # scores for the 15-mer peptides which contains these epitopes.
            # So, the first peptide "TQADMAAQTTKHKWE" contains all but the
            # the last epitope, and the second peptide "QADMAAQTTKHKWEA"
            # contains all but the first epitope.
            #

           
            epitope_mask = \
                (epitopes_df.SourceSequence == seq) & \
                (epitopes_df.EpitopeStart >= peptide_start) & \
                (epitopes_df.EpitopeEnd <= peptide_end) 
            epitope_subset = epitopes_df[epitope_mask]
            row['epitopes'] = epitope_subset

            records.append(row)
    assert len(records) > 0, "No vaccine peptides"
    return pd.DataFrame.from_records(records)





