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