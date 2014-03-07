# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. 

from sklearn import ensemble    

from pipeline import PipelineElement
from epitopes import iedb    

class ImmunogenicityRFModel(PipelineElement):
    def __init__(self, name, n_trees= 50, 
                assay_group = 'cytotoxicity',
                 # 'drop' | 'keep' | 'positive' | 'negative'
                human = True,
                hla_type = 1,
                ngram = 1, 
                reduced_alphabet = None):
        self.name = name
        self._X, self._Y, self._feature_transformer = iedb.load_tcell_dataset(
            assay_group=assay_group,
            human = human, 
            ngram = ngram, 
            reduced_alphabet = reduced_alphabet)
        self._model = ensemble.RandomForestClassifier(n_trees)
        self._model.fit(self._X, self._Y)

    def verify(self):
        pass

    def _apply(self, data):
        transformed_data = self._feature_transformer.transform(data.peptide).todense()
        return [x[1] for x in self._model.predict_proba(transformed_data)]
