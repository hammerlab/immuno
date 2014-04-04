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


def build_peptides_dataframe(
        epitopes_df,
        peptide_length,
        min_peptide_length = None):
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

    Returns a new dataframe with columns:
        - 'Peptide': longer amino acid sequence composed of multiple epitopes
        - 'MutatedScore' : median score of epitopes with mutations
        - 'SelfScore' : median score of epitopes without mutations
        - 'Score' : 'MutatedScore' - 'SelfScore'
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
        'SourceSequence', 'MutationStart', 'MutationEnd',
        'chr', 'pos', 'ref', 'alt', 'info', 'stable_id_transcript']
    for (seq, mut_start, mut_stop, chr, pos, ref, alt, info, t_id), rows \
    in epitopes_df.groupby(cols):
        peptides = []
        mut_scores = []
        self_scores = []
        mut_starts = []
        mut_ends = []

        n = len(seq)
        if n >= peptide_length:
            window_size = peptide_length
        elif n >= min_peptide_length:
            window_size = min(n, min_peptide_length)
        else:
            logging.info("Skipping source sequence %s from %s, shorter than %d",
                seq, min_peptide_length)
            continue
        for i in xrange(n + 1 - window_size):
            peptide_start = i
            peptide_end = i + window_size
            peptide = seq[peptide_start : peptide_end]

            # where is the mutation relative to this peptide?
            peptide_mut_start = min(mut_start - peptide_start, 0)
            peptide_mut_end = max(mut_end - peptide_start, window_size)

            # if mutation isn't in the peptide, skip it
            if peptide_mut_start >= window_size or peptide_mut_end <= 0:
                logging.info(
                    "Skipping self peptide %s from %s(transcript %s)[%d:%d]",
                    peptide,
                    info,
                    t_id,
                    peptide_start,
                    peptide_end)
                continue

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

            # which epitopes are from this sequence which overlap the mutations?
            mutated_epitope_mask = \
                (df.SourceSequence == seq) & \
                (df.EpitopeStart < mut_end) & \
                (df.EpitopeEnd > mut_start)

            mutated_epitopes = df[mutated_epitope_mask]
            mutated_epitope_score = mutated_epitopes['combined_score'].median()

            # self epitopes are those that overlap with the peptide but not
            # with the mutated epitopes
            self_epitope_mask = \
                (df.SourceSequence == seq) & \
                (df.EpitopeStart < peptide_end) & \
                (df.EpitopeEnd > peptide_start) & \
                ~mutated_epitope_mask

            self_epitopes = df[self_epitope_mask]
            self_epitope_score = self_epitopes['combined_score'].median()

            peptides.append(peptide)

            mut_scores.append(mutated_epitope_score)
            self_scores.append(self_epitope_score)
            mut_starts.append(peptide_mut_start)
            mut_ends.append(peptide_mut_start)






