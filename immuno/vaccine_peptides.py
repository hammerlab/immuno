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
from collections import namedtuple

import pandas as pd
import numpy as np

from epitope_scoring import simple_ic50_epitope_scorer
from peptide_binding_measure import IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
from immunogenicity import THYMIC_DELETION_FIELD_NAME


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

def is_mutant_epitope(epitope, mutation_start, mutation_end):
    """
    An epitope is considered mutant if it overlaps the mutated region
    of the source sequence and isn't similar to the thymically
    presented self epitopes.

    Parameters
    ----------

    epitope : dict

    mutation_start : int

    mutation_end : int
    """
    start = epitope['EpitopeStart']
    end = epitope['EpitopeEnd']
    overlaps = (start < mutation_end) and (end > mutation_start)
    if THYMIC_DELETION_FIELD_NAME in epitope:
        return not epitope[THYMIC_DELETION_FIELD_NAME] and overlaps
    else:
        return overlaps

def optimize_vaccine_peptide(
        seq,
        epitopes,
        mutation_start,
        mutation_end,
        epitope_scorer,
        result_length=31,
        padding=5):
    """
    Parameters
    ----------

    seq : str

    epitopes : list
        List of epitope records, each containing a nested list of per-allele
        binding predictions

    mutation_start : int
        Where in the given sequence is the first mutated residue?

    mutation_end : int
        Where in the given sequence is the last mutated residue?

    epitope_scorer : EpitopeScorer

    result_length : int
        How big of a substring are we looking to pull out as a vaccine peptide?

    padding : int
    """
    n = len(seq)

    if n <= result_length:
        # if source sequence is too short, just return whatever we have
        start = 0
        n_candidates = 1
    elif n <= result_length + 2 * padding:
        # if the source sequence is too short for the full amount of requested
        # padding, then center as best as we can
        actual_combined_padding = n - result_length
        start = actual_combined_padding / 2
        # if there are two equally good ways to center the insufficiently
        # padded sequence, try them both
        n_candidates = 1 if actual_combined_padding % 2 == 1 else 2
    else:
        start = padding
        n_candidates = n - result_length - 2 * padding

    # in case the mutation is at the beginning or end of the peptide,
    # make sure we cover it
    if mutation_start < start:
        difference = start - mutation_start
        start = mutation_start
        n_candidates += difference

    if mutation_start > start + result_length + n_candidates:
        difference = mutation_start - (start + result_length + n_candidates)
        n_candidates += difference

    # we're going to lexically sort each peptide by four criteria:
    #   - average score of its mutated epitopes
    #   - negative average score of it's wildtype epitopes
    #   - number of mutated residues covered
    #   - distance from the edge of the spurce sequence
    candidate_peptides = []

    for i in xrange(start, start+n_candidates):

        peptide_seq = seq[i:i+result_length]
        peptide_seq_len = len(peptide_seq)
        end = start + peptide_seq_len

        number_mutant_residues = \
            min(mutation_end, end) - max(mutation_start, start)
        peptide_mutation_start = mutation_start - i
        peptide_mutation_end = mutation_end - i

        half_len = peptide_seq_len / 2
        mutation_distance_from_edge = min(
            peptide_mutation_start,
            peptide_seq_len - peptide_mutation_start)

        mutant_score = 0.0
        wildtype_score = 0.0
        for epitope in epitopes:
            score = epitope_scorer.epitope_score(epitope)
            if is_mutant_epitope(epitope, mutation_start, mutation_end):
                mutant_score += score
            else:
                wildtype_score += score

        vaccine_peptide_record = {
            'VaccinePeptide' : peptide_seq,
            'VaccinePeptideMutationStart' : peptide_mutation_start,
            'VaccinePeptideMutationEnd' : peptide_mutation_end,
            'MutantEpitopeScore' : mutant_score,
            'VaccinePeptideStart' : i,
            'WildtypeEpitopeScore' : wildtype_score,
            'NumMutantResidues' : number_mutant_residues,
            'MutationDistanceFromEdge' : mutation_distance_from_edge,
        }
        candidate_peptides.append(vaccine_peptide_record)

    def score_tuple(record):
        """
        Create tuple of scores so that candidates get sorted lexicographically
        by multiple criteria. Make sure to make the wildtype epitope
        score negative (since we want fewer wildtype epitopes)
        """
        return (
            record['MutantEpitopeScore'],
            -record['WildtypeEpitopeScore'],
            record['NumMutantResidues'],
            record['MutationDistanceFromEdge'],
        )
    candidate_peptides.sort(key=score_tuple, reverse=True)
    best = candidate_peptides[0]
    return best

def select_vaccine_peptides(
        source_peptides,
        epitope_scorer=simple_ic50_epitope_scorer,
        vaccine_peptide_length=31,
        padding=5):
    """
    Given a set of longer peptides and their associated predicted epitopes,
    find the best vaccine peptide overlapping each mutation.

    Parameters
    ----------

    source_peptides : list
        List of peptide source sequence records

    epitope_scorer : EpitopeScorer instance

    vaccine_peptide_length : int

    padding : int
        Maximum distance from edges of vaccine peptide where mutation can start.
    """

    results = []
    for peptide_record in source_peptides:
        seq = peptide_record['SourceSequence']
        assert len(seq) > 0, "Invalid empty peptide"
        mutation_start = peptide_record['MutationStart']
        mutation_end = peptide_record["MutationEnd"]
        epitopes = peptide_record['Epitopes']

        vaccine_peptide_record = optimize_vaccine_peptide(
            seq,
            epitopes,
            mutation_start,
            mutation_end,
            epitope_scorer,
            result_length=vaccine_peptide_length,
            padding=padding,
        )
        start_idx = vaccine_peptide_record['VaccinePeptideStart']

        assert start_idx >= 0
        # must overlap the mutation to some degree
        assert start_idx <= mutation_end

        # augment the vaccine peptide with all info that was attached to its
        # source sequence
        for k,v in peptide_record.iteritems():
            if k not in vaccine_peptide_record:
                vaccine_peptide_record[k] = v

        results.append(vaccine_peptide_record)

    # Make sure that sort is in descending order by vaccine peptide score.
    # When comparing across genes/mutations we only care about the mutant
    # epitope score.
    results.sort(key=lambda record: record['MutantEpitopeScore'], reverse=True)
    return results
