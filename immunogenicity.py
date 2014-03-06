"""

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""
from pipeline import PipelineElement
import iedb
from sklearn import ensemble

class ImmunogenicityRFModel(PipelineElement):
  def __init__(self, name, n_trees= 50, filename = 'tcell_compact.csv',
                 assay_group = 'cytotoxicity',
                 # 'drop' | 'keep' | 'positive' | 'negative'
                 noisy_labels = 'majority',
                 human = True,
                 hla_type1 = True,
                 exclude_hla_a2 = False,
                 only_hla_a2 = False,
                 max_ngram = 1, 
                 normalize_row = True, 
                 reduced_alphabet = None):
    self.name = name
    self._X, self._Y, self._feature_transformer = iedb.load_dataset(filename = filename,
                 assay_group=assay_group,
                 # 'drop' | 'keep' | 'positive' | 'negative'
                 noisy_labels = noisy_labels, 
                 human = human, 
                 hla_type1 = hla_type1,
                 exclude_hla_a2 = exclude_hla_a2, 
                 only_hla_a2 = only_hla_a2, 
                 max_ngram = max_ngram, 
                 normalize_row = normalize_row, 
                 reduced_alphabet = reduced_alphabet)
    self._model = ensemble.RandomForestClassifier(n_trees)
    self._model.fit(self._X, self._Y)

  def verify(self):
    pass
    #super(self)

  def _apply(self, data):
    transformed_data = self._feature_transformer.transform(data.peptide).todense()
    return [x[1] for x in self._model.predict_proba(transformed_data)]