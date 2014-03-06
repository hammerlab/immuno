"""

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""

from matrices_pcm import pcm_matrices
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


