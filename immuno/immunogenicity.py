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

import cPickle
import logging
from os.path import exists

from sklearn import ensemble
import numpy as np
import pandas as pd

from epitopes import \
    (cri_tumor_antigens, iedb, features, reduced_alphabet, reference)

from pipeline import PipelineElement
from balanced_ensemble import BalancedEnsembleClassifier

def train(
        assay_group = 'cytotoxicity',
        mhc_class = 1,
        max_ngram = 2,
        alphabet = 'sdm12'):
    params_str = "assay='%s', mhc=%s, ngram=%s, alphabet=%s" % \
        (assay_group, mhc_class, max_ngram, alphabet)
    logging.info("Training immunogenicity classifier (%s)" % params_str)
    alphabet_dict = getattr(reduced_alphabet, alphabet)
    X, Y, vectorizer = iedb.load_tcell_ngrams(
        assay_group = assay_group,
        human = True,
        mhc_class = mhc_class,
        max_ngram = max_ngram,
        reduced_alphabet = alphabet_dict,
        min_count = None,
        return_transformer = True)
    ensemble = BalancedEnsembleClassifier()
    ensemble.fit(X, Y)
    return vectorizer, ensemble


def train_cached(
        assay_group = 'cytotoxicity',
        mhc_class = 1,
        max_ngram = 2,
        alphabet = 'sdm12',
        suffix = '.pickle'):
    """
    Trains and caches both a feature transforming vectorizer
    which turns peptide strings into fixed-length vectors
    and a classifier which maps vectors to immunogenicity labels
    """

    model_base_filename = \
        "immunogenicity_classifier_assay_%s_mhc_%s_ngram_%s_alphabet_%s" % \
        (assay_group, mhc_class, max_ngram, alphabet)
    model_filename = model_base_filename + suffix

    vectorizer_base_filename = \
        "immunogenicity_vectorizer_assay_%s_mhc_%s_ngram_%s_alphabet_%s" % \
        (assay_group, mhc_class, max_ngram, alphabet)
    vectorizer_filename = vectorizer_base_filename + suffix

    if exists(model_filename) and exists(vectorizer_filename):
        logging.debug("Loading cached immunogenicity classifier and vectorizer")
        with open(vectorizer_filename, 'r') as vectorizer_file:
            vectorizer = cPickle.load(vectorizer_file)

        with open(model_filename, 'r') as model_file:
            clf = cPickle.load(model_file)
        return vectorizer, clf

    vectorizer, clf = train()

    with open(vectorizer_filename, 'w') as vectorizer_file:
        cPickle.dump(vectorizer, vectorizer_file, cPickle.HIGHEST_PROTOCOL)

    with open(model_filename, 'w') as model_file:
        cPickle.dump(clf, model_file, cPickle.HIGHEST_PROTOCOL)

    return vectorizer, clf

class ImmunogenicityRFModel(PipelineElement):

    def __init__(
            self, name = "immunogenicity",
            classifier = None, vectorizer = None):
        self.name = name
        if classifier is None and vectorizer is None:
            vectorizer, classifier = train_cached()
        self.vectorizer = vectorizer
        self.classifier = classifier


    def verify(self):
        pass

    def _apply(self, df):
        if self.vectorizer:
            X = self.vectorizer.transform(df.Epitope)
        else:
            X = df.Epitope
        return self.classifier.decision_function(X)
