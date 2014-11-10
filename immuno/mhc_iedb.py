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

from mhc_common import normalize_hla_allele_name, seq_to_str, convert_str
from peptide_binding_measure import (
        IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
)


class IEDBMHCBinding(object):

  def __init__(
        self,
        alleles=[],
        name="IEDB-MHC-Binding",
        method=['recommended'],
        lengths = [9],
        url='http://tools.iedb.org/tools_api/mhci/'):
    self.name = name
    self._method = method
    self._lengths = lengths
    self._url = url
    self._alleles = alleles

  def _get_iedb_request_params(self, sequence, allele=None):
    if not allele:
        allele = seq_to_str(self._alleles)
    params = {
        "method" : seq_to_str(self._method),
        "length" : seq_to_str(self._lengths),
        "sequence_text" : sequence,
        "allele" : allele,
    }
    return params

  def query_iedb(self, sequence, gene_info, allele=None):
    request_values = self._get_iedb_request_params(sequence, allele=allele)
    logging.info("Calling IEDB with {} {}, {}".format(
        gene_info, sequence, allele))
    try:
        data = urllib.urlencode(request_values)
        req = urllib2.Request(self._url, data)
        response = urllib2.urlopen(req).read()
        lines = response.split("\n")

        # manually parsing since pandas is insane
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
    except KeyboardInterrupt:
        raise
    except:
        raise
        logging.error(
            "Connection error: Failed on sequence {}".format(sequence))
        return pd.DataFrame()


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
        for allele in self._alleles:
            key = (peptide, allele)
            if key not in responses:
                gene_info = data['GeneInfo'][i]
                response = self.query_iedb(peptide, gene_info, allele)
                response.rename(
                    columns={
                        'peptide': 'Epitope',
                        'length' : 'EpitopeLength',
                        'start' : 'EpitopeStart',
                        'end' : 'EpitopeEnd',
                        'allele' : 'Allele',
                    },
                    inplace=True)
                response['EpitopeStart'] -= 1
                response['EpitopeEnd'] -= 1
                responses[key] = response
            else:
                logging.info(
                    "Skipping binding for peptide %s / allele %s, already queried",
                    peptide, allele)

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

class IEDBMHC1Binding(IEDBMHCBinding):
    def __init__(self,
            name = 'IEDB-MHC1-Binding',
            url='http://tools.iedb.org/tools_api/mhci/',
            alleles=[]):

        super(IEDBMHC1Binding, self).__init__(
            name = name,
            url = url,
            alleles = alleles)


class IEDBMHC2Binding(IEDBMHCBinding):
    def __init__(self,
            name = 'IEDB-MHC2-Binding',
            url='http://tools.iedb.org/tools_api/mhcii/',
            alleles=[]):

      super(IEDBMHC2Binding, self).__init__(
        name = name,
        url = url,
        alleles = alleles)

    def _get_iedb_request_params(self, sequence):
      params = {
        "method" : seq_to_str(self._method),
        "sequence_text" : sequence,
        "allele" : seq_to_str(self._alleles),
      }
      return params
