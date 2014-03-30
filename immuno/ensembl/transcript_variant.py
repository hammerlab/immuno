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
            return str(mutate(transcript.seq, pos, ref, alt))
        except:
            raise
            return None
    return None

def peptide_from_transcript_variant(
        transcript_id, pos, ref, alt,
        max_length = 500,
        min_padding=31):
    transcript = _ensembl.get_cdna(transcript_id)
    if transcript:

        idx = annotation.get_transcript_index_from_pos(pos, transcript_id)

        if idx is not None:
            try:
                mutated = mutate_protein_from_transcript(
                    transcript.seq, idx, ref, alt,
                    max_length = max_length,
                    min_padding = min_padding)
                return str(mutated)
            except AssertionError, error:
                return None
    return None
