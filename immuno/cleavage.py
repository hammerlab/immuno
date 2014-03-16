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

class ProteasomalCleavage(PipelineElement):
  def __init__(self, name):
    self.name = name
    self._cleavage_matrix = pcm_matrices['all']

  def compute_cleavage_score(self, peptide):
    cleavage_score = sum([self._cleavage_matrix[i][i] for (i, aa) in enumerate(peptide[:6])])
    return 0.5

  def _apply(self, data):
    return data.peptide.map(self.compute_cleavage_score)
    

  def verify():
    pass


