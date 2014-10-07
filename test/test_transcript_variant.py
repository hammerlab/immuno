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

from immuno.ensembl import transcript_variant

protein_id = "ENSP00000427553"
protein_sequence = \
"""MRLPGAPALPDADFLVHLHFLVQTSWFICNFLVRIPWASALPDAPALLVILEKTFPEHATCRGCWVSGYLCWTAPGNCICSANVGFLKIENTYRQIHHTHMHRHTHTHTQTNPSHTHAQTHTHRVNDIGSQVELFVCLYLMQLLIHLSLELLFSFTYTVCLQILCINSFLSCRLLDDFPQLTLRTFEQTDTLKE"""

def test_peptide_from_protein_transcript():
    peptide = \
        transcript_variant.peptide_from_protein_transcript_variant(
            protein_id, 109, 'Q', 'R')
    assert peptide is not None
    assert peptide[109] == 'R', peptide

def test_peptide_from_transcript_variant_CTNNB1():
    """
    test_peptide_from_transcript:

    Apply Cosmic mutation COSM27279
    transcript = 'ENST00000405570'
    pos: 41265571,
    ref : A, alt : T
    amino acids = Q -> H  @ pos 4 (mutation = Q4H)
    """
    transcript_id = 'ENST00000405570'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            41265571,
            ref='A',
            alt='T',
            padding = None,
            max_length = None)
    assert peptide is not None
    n = len(peptide)
    assert n == 781, (n, peptide)
    assert(peptide[3] == 'H')

def test_peptide_from_transcript_variant_CASP9():
    transcript_id = 'ENST00000333868'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            15820483,
            ref='C',
            alt='G',
            padding = None,
            max_length = None)
    assert peptide is not None
    assert len(peptide) == 416
    # TODO: actually look up what this variant ought to be

def test_peptide_from_transcript_variant_CASP9_ref():
    transcript_id = 'ENST00000424908'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            15820483,
            ref='C',
            alt='G',
            padding = None,
            max_length = None)
    assert peptide is not None
    assert len(peptide) == 198
    assert peptide[0] == 'X'

def test_peptide_from_transcript_variant_VWA3A():
    transcript_id = 'ENST00000389397'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            22128096,
            ref='G',
            alt='A',
            padding = None,
            max_length = None)
    # expected mutation in the UTR to yield None
    assert peptide is None

def test_peptide_from_transcript_variant_RET():
    transcript_id = 'ENST00000340058'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            43617416,
            ref='T',
            alt='C',
            padding = None,
            max_length = None)
    assert peptide is not None

def test_peptide_from_transcript_variant_PAR():

    transcript_id = 'ENST00000371279'
    peptide, start, stop, annot = \
        transcript_variant.peptide_from_transcript_variant(
            transcript_id,
            55224569,
            ref='T',
            alt='G',
            padding = None,
            max_length = None)
    assert peptide is not None


if __name__ == '__main__':
  from dsltools import testing_helpers
  testing_helpers.run_local_tests()