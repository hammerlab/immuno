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
from os.path import exists, join 

from sklearn import ensemble
import numpy as np
import pandas as pd
import scipy.stats 

from epitopes import \
    (cri_tumor_antigens, iedb, features, reduced_alphabet, reference)

from balanced_ensemble import BalancedEnsembleClassifier
from common import MODELS_DIR

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

    reference_peptides = reference.load_peptide_set(peptide_length = 9, nrows = 1000)
    # filter out stop codons and unknown residues
    reference_peptides = [p for p in reference_peptides if '*' not in p and 'X' not in p]
    reference_peptide_vectors = vectorizer.transform(reference_peptides)
    reference_scores = ensemble.decision_function(reference_peptide_vectors)

    return vectorizer, ensemble, reference_scores 


def make_path(base_filename, suffix):
    filename = base_filename + suffix
    path = join(MODELS_DIR, filename)
    return path

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
    model_path = make_path(model_base_filename, suffix)

    vectorizer_base_filename = \
        "immunogenicity_vectorizer_assay_%s_mhc_%s_ngram_%s_alphabet_%s" % \
        (assay_group, mhc_class, max_ngram, alphabet)
    vectorizer_path = make_path(vectorizer_base_filename, suffix)

    reference_scores_base_filename = \
        "immunogenicity_reference_scores_assay_%s_mhc_%s_ngram_%s_alphabet_%s" % \
        (assay_group, mhc_class, max_ngram, alphabet)

    reference_scores_path = make_path(reference_scores_base_filename, suffix)

    if all(map(exists, [model_path, vectorizer_path, reference_scores_path])):
        logging.debug("Loading cached immunogenicity classifier and vectorizer")
        with open(vectorizer_path, 'r') as vectorizer_file:
            vectorizer = cPickle.load(vectorizer_file)

        with open(model_path, 'r') as model_file:
            clf = cPickle.load(model_file)

        with open(reference_scores_path, 'r') as scores_file:
            reference_scores = cPickle.load(scores_file)
        return vectorizer, clf, reference_scores 

    vectorizer, clf, reference_scores = train()

    with open(vectorizer_path, 'w') as vectorizer_file:
        cPickle.dump(vectorizer, vectorizer_file, cPickle.HIGHEST_PROTOCOL)

    with open(model_path, 'w') as model_file:
        cPickle.dump(clf, model_file, cPickle.HIGHEST_PROTOCOL)


    with open(reference_scores_path, 'w') as scores_file:
        cPickle.dump(reference_scores, scores_file, cPickle.HIGHEST_PROTOCOL)

    return vectorizer, clf, reference_scores 

class ImmunogenicityRFModel(PipelineElement):

    def __init__(
            self, 
            name = "immunogenicity",
            classifier = None, 
            vectorizer = None):
        self.name = name
        if classifier is None and vectorizer is None:
            vectorizer, classifier, reference_scores = train_cached()
        self.vectorizer = vectorizer
        self.classifier = classifier
        self.reference_scores = reference_scores 

    def apply(self, df):
        if self.vectorizer:
            X = self.vectorizer.transform(df.Epitope)
        else:
            X = df.Epitope
        probs = self.classifier.decision_function(X)
        percentiles = []
        for prob in probs:
            prctile = scipy.stats.percentileofscore(self.reference_scores, prob)
            percentiles.append(prctile)
        df[self.name] = np.array(percentiles) / 100.0
        return df