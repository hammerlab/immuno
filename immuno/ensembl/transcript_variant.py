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

"""
Apply a mutation to a transcript and return a window of amino acids around
the mutated residue.
"""

import logging

from epitopes.mutate import mutate_protein_from_transcript, mutate

from transcript_data import EnsemblReferenceData
import annotation

_ensembl = EnsemblReferenceData()

def peptide_from_protein_transcript_variant(transcript_id, pos, ref, alt):
    """
    Given an ensembl transcript ID, mutate amino acid `ref` to `alt` at
    position `pos`.
    """
    transcript = _ensembl.get_protein(transcript_id)
    if transcript:
        try:
            return str(mutate(transcript, pos, ref, alt))
        except:
            logging.warning(
                "Failed to mutate transcript %s (ref %s, alt %s at pos %s)",
                ref,
                alt,
                pos)
            return None
    return None

def peptide_from_transcript_variant(
        transcript_id, pos, ref, alt,
        max_length = 500,
        min_padding=31):
    transcript = _ensembl.get_cdna(transcript_id)
    if not transcript:
        logging.warning("Couldn't find transcript for ID %s", transcript_id)
        return None, -1, -1
    idx = annotation.get_transcript_index_from_pos(pos, transcript_id)
    if idx is None:
        logging.warning(
            "Couldn't translate gene position %s into transcript index for %s",
            pos,
            transcript_id)
        return None, -1, -1
    try:
        mutated, start, stop = mutate_protein_from_transcript(
            transcript, idx, ref, alt,
            max_length = max_length,
            min_padding = min_padding,
            with_mutation_coordinates=True)
        return str(mutated), start, stop
    except AssertionError, error:
        logging.warning(
            "Failed to mutate %s (ref %s, alt %s at position %s)",
            transcript_id,
            ref,
            alt,
            pos)
        return None, -1, -1
