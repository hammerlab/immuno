# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



from ensembl_transcript_data import EnsemblReferenceData
import ensembl_annotation
from snpeff_effect import SnpEffEffect


"""
def peptide_from_protein_transcript(transcript_id, pos, ref, alt):
    transcript = _ensembl.get_protein(transcript_id)
    if transcript:
        try:
            return str(mutate.mutate(transcript.seq, pos, ref, alt))
        except:
            return None
    return None

def parse_aa_sift(self, sift_string):
    # ENSP00000427553:T109R
    protein_transcript_id, variation = sift_string.split(":")
    aa_orig = variation[0]
    aa_variant = variation[-1]
    aa_change_pos = int(variation[1:-1])
    return (protein_transcript_id, aa_orig, aa_variant, aa_change_pos)
"""