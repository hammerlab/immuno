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

import urllib2
import urllib
from StringIO import StringIO
import logging
import re

import pandas as pd

from mhc_base_predictor import MHCBasePredictor
from mhc_common import normalize_hla_allele_name, seq_to_str, convert_str
from peptide_binding_measure import (
        IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
)

"""
A note about prediction methods, copied from the IEDB website:

The prediction method list box allows choosing from a number of MHC class I
binding prediction methods:
- Artificial neural network (ANN),
- Stabilized matrix method (SMM),
- SMM with a Peptide:MHC Binding Energy Covariance matrix (SMMPMBEC),
- Scoring Matrices  from Combinatorial Peptide Libraries (Comblib_Sidney2008),
- Consensus,
- NetMHCpan.

IEDB recommended is the default prediction method selection.
Based on availability of predictors and previously observed predictive
performance, this selection tries to use the best possible method for a given
MHC molecule. Currently for peptide:MHC-I binding prediction, for a given MHC
molecule, IEDB Recommended uses the Consensus method consisting of ANN, SMM,
and CombLib if any corresponding predictor is available for the molecule.
Otherwise, NetMHCpan is used. This choice was motivated by the expected
predictive performance of the methods in decreasing order:
    Consensus > ANN > SMM > NetMHCpan > CombLib.
"""

VALID_IEDB_METHODS = [
    'recommended',
    'consensus',
    'netmhcpan',
    'ann',
    'smmpmbec',
    'smm',
    'comblib_sidney2008'
]


def _parse_iedb_response(response):
    """
    Take the binding predictions returned by IEDB's web API
    and parse them into a DataFrame
    """
    lines = response.split("\n")

    # manually parsing since Pandas is insane
    header_names = lines[0].split("\t")

    d = {}
    for col_name in header_names:
        d[col_name] = []

    for line in lines[1:]:
        line = line.strip()
        if len(line) > 0:
            fields = line.split('\t')
            for i, header_name in enumerate(header_names):
                value = convert_str(fields[i] if len(fields) > i else None)
                d[header_name].append(value)
    return pd.DataFrame(d)


def _query_iedb(request_values, url):
    """
    Call into IEDB's web API for MHC binding prediction using request dictionary
    with fields:
        - "method"
        - "length"
        - "sequence_text"
        - "allele"

    Parse the response into a DataFrame.
    """
    data = urllib.urlencode(request_values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req).read()
    return _parse_iedb_response(response)


class IEDB_MHC_Binding_Predictor(MHCBasePredictor):

  def __init__(
        self,
        hla_alleles,
        epitope_lengths,
        method,
        url):
    MHCBasePredictor.__init__(
        self,
        hla_alleles=hla_alleles,
        epitope_lengths=epitope_lengths)

    if method not in VALID_IEDB_METHODS:
        raise ValueError(
            "Invalid IEDB MHC binding prediction method: %s" % (method,))

    self.method = method

    if not isinstance(url, str):
        raise TypeError("Expected URL to be string, not %s : %s" % (
            url, type(url)))
    self.url = url


  def _get_iedb_request_params(self, sequence, allele):

    params = {
        "method" : seq_to_str(self.method),
        "length" : seq_to_str(self.epitope_lengths),
        "sequence_text" : sequence,
        # have to repeat allele for each length
        "allele" : ",".join([allele] * len(self.epitope_lengths)),
    }
    return params


  def predict(self, data):
    """
    Given a dataframe with long amino acid sequences in the
    'SourceSequence' field, return an augmented dataframe
    with shorter k-mers in the 'Epitope' column and several
    columns of MHC binding predictions with names such as 'percentile_rank'
    """
    # take each mutated sequence in the dataframe
    # and general MHC binding scores for all k-mer substrings
    responses = {}
    for i, peptide in enumerate(data.SourceSequence):
        for allele in self.alleles:
            key = (peptide, allele)
            if key not in responses:
                request = self._get_iedb_request_params(peptide, allele)
                logging.info(
                    "Calling IEDB (%s) with request %s",
                    self._url,
                    request)
                response_df = _query_iedb(request, self._url)
                response_df.rename(
                    columns={
                        'peptide': 'Epitope',
                        'length' : 'EpitopeLength',
                        'start' : 'EpitopeStart',
                        'end' : 'EpitopeEnd',
                        'allele' : 'Allele',
                    },
                    inplace=True)
                response_df['EpitopeStart'] -= 1
                response_df['EpitopeEnd'] -= 1
                responses[key] = response_df
            else:
                logging.info(
                    "Already made predictions for peptide %s with allele %s",
                    peptide,
                    allele)

    # concatenating the responses makes a MultiIndex with two columns
    # - SourceSequence
    # - index of epitope from that sequence's IEDB call
    #
    # ...when we reset the index, we turn these into two columns
    # named 'level_0', and 'level_1'. We want to rename the former
    # and delete the latter.
    responses = pd.concat(responses).reset_index()
    responses['SourceSequence'] = responses['level_0']
    del responses['level_0']
    del responses['level_1']

    # IEDB has inclusive end positions, change to exclusive
    responses['EpitopeEnd'] += 1

    assert 'ann_rank' in responses, responses.head()
    responses[PERCENTILE_RANK_FIELD_NAME] = responses['ann_rank']

    assert 'ann_ic50' in responses, responses.head()
    responses[IC50_FIELD_NAME] = responses['ann_ic50']

    # instead of just building up a new dataframe I'm expliciting
    # dropping fields here to document what other information is available
    drop_fields = (
        'seq_num',
        'method',
        'ann_ic50',
        'ann_rank',
        'consensus_percentile_rank',
        'smm_ic50',
        'smm_rank',
        'comblib_sidney2008_score',
        'comblib_sidney2008_rank'
    )
    for field in drop_fields:
        if field in responses:
            responses = responses.drop(field, axis = 1)

    result = data.merge(responses, on='SourceSequence')

    # some of the MHC scores come back as all NaN so drop them
    result = result.dropna(axis=1, how='all')

    return result

class IEDB_MHC1(IEDB_MHC_Binding_Predictor):
    def __init__(self,
        alleles,
        epitope_lengths=[9],
        method='recommended',
        url='http://tools-api.iedb.org/tools_api/mhci/'):

        IEDB_MHC_Binding_Predictor.__init__(
            self,
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            method=method,
            url=url)

class IEDB_MHC2(IEDB_MHC_Binding_Predictor):
    def __init__(self,
            alleles,
            method='recommended',
            url='http://tools-api.iedb.org/tools_api/mhcii/'):

      IEDB_MHC_Binding_Predictor.__init__(
        self,
        alleles=alleles,
        # only epitope lengths of 15 currently supported by IEDB's web API
        epitope_lengths=[15],
        method=method,
        url=url)

