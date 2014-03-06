"""

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""
from sklearn import ensemble

from pipeline import PipelineElement
from epitopes import iedb

class ImmunogenicityRFModel(PipelineElement):
  def __init__(self, name, n_trees= 50, 
                 assay_group = 'cytotoxicity',
                 # 'drop' | 'keep' | 'positive' | 'negative'
                 human = True,
                 hla_type = 1,
                 max_ngram = 1, 
                 reduced_alphabet = None):
    self.name = name
    self._X, self._Y, self._feature_transformer = iedb.load_tcell(
        assay_group=assay_group,
        noisy_labels = noisy_labels, 
        human = human, 
        mhc_class = mhc_class,
        max_ngram = max_ngram, 
        reduced_alphabet = reduced_alphabet)
    self._model = ensemble.RandomForestClassifier(n_trees)
    self._model.fit(self._X, self._Y)

  def verify(self):
    pass
    #super(self)

  def _apply(self, data):
    transformed_data = self._feature_transformer.transform(data.peptide).todense()
    return [x[1] for x in self._model.predict_proba(transformed_data)]
