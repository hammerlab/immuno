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

from immuno import epitope_generation
    

def test_generate_peptide_from_protein_transcript():
    epitope = epitope_generation.generate_peptide_from_protein_transcript('ENSP00000427553', 509, 'T', 'R')

    print(str(epitope))

def test_vcf2dataframe():
    vcf_file = 'example.vcf'
    df = epitope_generation.vcf2dataframe(vcf_file)

    assert(len(df) == 3)
    assert(len(df.columns) == 8)

def test_generate_peptide_from_transcript():
    """
    Example:
    'ENST00000405570'
    pos: 41275636, ref : G, alt : A

    should return MAQNAVRLHYGLPVVVKLLHPPSHWPLIKATIGLIRNLALCPANHAPLREQGAIPRLVQLLVR
    """
    transcript_id = 'ENST00000405570'
    peptide = epitope_generation.generate_peptide_from_transcript(transcript_id, 41275636, ref='G', alt='A')
    n = len(peptide)

    print(str(peptide))
    assert(peptide[n/2] == 'I')



