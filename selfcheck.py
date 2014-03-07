# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from pipeline import PipelineElement
from snpeff_effect import Ensembl


class SelfTypePeptide(PipelineElement):
  def __init__(self, name='selftypepeptide', protein_path='Homo_sapiens.GRCh37.74.pep.all.fa', index_lengths = [8,9,10,11]):
    self.name = name
    self._ensembl = Ensembl(protein_path=protein_path, index_peptides=True, index_lengths = index_lengths)

  def peptide_in_reference(self, peptide):
    self._ensembl.reference_contains_peptide(peptide)

  def _apply(self, data):
    return data.peptide.map(self.peptide_in_reference)




if __name__ == '__main__':
  selfcheck = SelfTypePeptide()
  selfcheck.peptide_in_reference("GRKAKGS")