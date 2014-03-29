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

from immuno import transcript_variant


def test_peptide_from_protein_transcript():
    peptide = \
        transcript_variant.peptide_from_protein_transcript_variant(
            'ENSP00000427553', 109, 'Q', 'R')
    assert peptide is not None


def test_peptide_from_transcript():
    """
    Example:
    'ENST00000405570'
    pos: 41275636, ref : G, alt : A

    should return
    MAQNAVRLHYGLPVVVKLLHPPSHWPLIKATIGLIRNLALCPANHAPLREQGAIPRLVQLLVR
    """
    transcript_id = 'ENST00000405570'
    peptide = transcript_variant.peptide_from_transcript_variant(
        transcript_id, 41275636, ref='G', alt='A')
    assert peptide is not None
    n = len(peptide)
    print(str(peptide))
    assert(peptide[n/2] == 'I')

if __name__ == '__main__':
  from dsltools import testing_helpers
  testing_helpers.run_local_tests()


