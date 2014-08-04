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

from epitopes.mutate import mutate_protein_from_transcript, mutate, gene_mutation_description

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
                "Failed to mutate transcript %s (%s)", 
                transcript_id,
                gene_mutation_description(pos, ref, alt)
            )
            return None
    return None

def peptide_from_transcript_variant(
        transcript_id, pos, ref, alt,
        padding = None,
        max_length = None):
     
    
    forward = annotation.is_forward_strand(transcript_id)
    ref = ref if forward else annotation.reverse_complement(ref)
    alt = alt if forward else annotation.reverse_complement(alt)
    transcript = _ensembl.get_cds(transcript_id)
    
    def error_result(msg, *args):
        logging.warning(msg, *args)
        return None, -1, -1, msg % args 

    if not transcript:
        return error_result("Couldn't find transcript for ID %s", transcript_id)

    idx = annotation.get_transcript_index_from_pos(
        pos, 
        transcript_id,
        skip_untranslated_region = True)
    if idx is None:
        return error_result(
            "Couldn't translate gene position %s into transcript index for %s",
            pos,
            transcript_id)
    elif idx >= len(transcript):
        return error_result(
            "Can't get position %d in coding sequence of length %d for transcript %s (%s)",
            idx, 
            len(transcript),
            transcript_id,
            gene_mutation_description(pos, ref, alt))
    try:      
        idx = idx if forward else idx - len(ref) + 1  
        region = mutate_protein_from_transcript(
            transcript,
            idx,
            ref,
            alt,
            padding = padding)
        start = region.mutation_start
        stop = start + region.n_inserted
        if max_length and len(region.seq) > max_length:
            seq = region.seq[:max_length]
            stop = min(stop, max_length)
        else:
            seq = region.seq
        return seq, start, stop, region.annot
    except:
        raise
   