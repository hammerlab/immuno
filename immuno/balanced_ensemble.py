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

from math import ceil

import numpy as np
import sklearn.ensemble
import sklearn.linear_model

from sklearn.base import ClassifierMixin

"""
Train multiple models for classification on an unbalanced dataset (one class
is more abundant than another. Implements the common "UnderBagging" approach
found to work well in Galar et al.'s "A Review on Ensembles for the Class Imbalance Problem: Bagging-, Boosting-, and Hybrid-Based Approaches", but instead of sampling with replacement (to attain diversity), just randomly subsample the minority class.
"""

def subsample_indices(n, k):
    idx = np.arange(n)
    np.random.shuffle(idx)
    return idx[:k]


class BalancedEnsembleClassifier(ClassifierMixin):
    def __init__(self, **kwargs):
        self.set_params(**kwargs)

    def get_params(self, deep=True):
        return {'n_estimators' : self.n_estimators}

    def set_params(self, **parameters):
        n_estimators = parameters.pop('n_estimators', 100)
        self.n_estimators = n_estimators
        assert len(parameters) == 0, \
            "Unexpected keywords %s" % (parameters.keys(),)


    def fit(self, X, Y):
        self.models = []
        n_total = len(Y)

        true_mask = Y > 0
        X_true = X[true_mask]
        Y_true = Y[true_mask]
        n_true = len(Y_true)

        false_mask = ~true_mask
        X_false = X[false_mask]
        Y_false = Y[false_mask]
        n_false = len(Y_false)

        n_min = min(n_true, n_false)

        # subsample even the minority class
        # to randomly drop outliers and
        # create more diverse population of models
        n_sub = int(0.9 * n_min)

        n_models = 2 * int(ceil(float(n_total) / n_sub)) + 1


        for i in xrange(n_models):

            # subsample the true and false datasets
            # and then concatenate both datasets into
            # a smaller X,Y with 2*n_sub rows

            true_idx = subsample_indices(n_true, n_sub)
            X_true_sub = X_true[true_idx]
            Y_true_sub = Y_true[true_idx]

            false_idx = subsample_indices(n_false, n_sub)
            X_false_sub  = X_false[false_idx]
            Y_false_sub = Y_false[false_idx]

            X_sub = np.vstack([X_true_sub, X_false_sub])
            Y_sub = np.concatenate([Y_true_sub, Y_false_sub])
            #clf = sklearn.linear_model.LogisticRegression()
            clf = sklearn.ensemble.RandomForestClassifier(
                n_estimators = self.n_estimators)
            clf.fit(X_sub, Y_sub)
            self.models.append(clf)

    def _predict_counts(self, X):
        Y = np.zeros(len(X), dtype=int)
        for model in self.models:
            Y += model.predict(X)
        return Y

    def _add_probs(self, X):
        Y = np.zeros(len(X), dtype=float)
        for model in self.models:
            Y += model.predict_proba(X)[:, 1]
        return Y

    def _mean_probs(self, X):
        probs = self._add_probs(X)
        probs /= len(self.models)
        return probs

    def decision_function(self, X):
        return self._mean_probs(X)

    def predict_proba(self, X):
        Y_prob = self.decision_function(X)
        return np.vstack([1.0 - Y_prob, Y_prob]).T

    def predict(self, X):
        Y_prob = self.decision_function(X)
        return Y_prob > 0.5

