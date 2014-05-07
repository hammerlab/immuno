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

from os.path import join 
import urllib2
import urllib
from StringIO import StringIO
import logging
import re

import pandas as pd

from pipeline import PipelineElement


def normalize_hla_allele_name(hla):
    """
    HLA allele names can look like:
        - HLA-A*03:02
        - HLA-A02:03
        - HLA-A:02:03
        - HLA-A2
        - A2 
        - A*03:02
        - A02:02
        - A:02:03
    ...should all be normalized to:
        HLA-A*03:02:03
    """
    hla = hla.strip().upper()
    match = re.match('(HLA\-)?([A-Z])(\*|:)?([0-9][0-9]?):?([0-9][0-9]?)$', hla)
    assert match, "Malformed HLA type %s" % hla 
    (_, gene, _, family, protein) = match.groups()
    if len(family) == 1:
        family = "0" + family 
    if len(protein) == 1:
        protein = "0" + protein 
    return "HLA-%s*%s:%s" % (gene, family, protein )

class IEDBMHCBinding(PipelineElement):

  def __init__(
        self, 
        alleles=[],
        name="IEDB-MHC-Binding",
        method='recommended',
        lengths = [9],
        url='http://tools.iedb.org/tools_api/mhci/'):
    self.name = name
    self._method = method
    self._lengths = lengths
    self._url = url
    self._alleles = alleles

  def _get_iedb_request_params(self, sequence):
    params = {
        "method" : self._method,
        "length" : ",".join(str(l) for l in self._lengths),
        "sequence_text" : sequence,
        "allele" : ",".join(self._alleles),
    }
    return params

  def query_iedb(self, sequence, gene_info):
    request_values = self._get_iedb_request_params(sequence)
    logging.info("Calling iedb with {} {}, {}".format(
        gene_info, sequence, self._alleles))
    try:
      data = urllib.urlencode(request_values)
      req = urllib2.Request(self._url, data)
      response = urllib2.urlopen(req).read()
      return pd.read_csv(StringIO(response), sep='\t', na_values=['-'])
    except KeyboardInterrupt:
        raise
    except:
        logging.error(
            "Connection error: Failed on sequence {}".format(sequence))
        return pd.DataFrame()


  def apply(self, data):
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
        if peptide not in responses:
            response = self.query_iedb(peptide, data['info'][i])
            responses[peptide] = response 
        else:
            logging.info(
                "Skipping binding for peptide %s, already queried",
                peptide)
    responses = pd.concat(responses).reset_index(0)
    responses.rename(
        columns={
            'level_0':'SourceSequence',
            'peptide': 'Epitope',
            'length' : 'EpitopeLength',
            'start' : 'EpitopeStart',
            'end' : 'EpitopeEnd',
        },
        inplace=True)

    result = data.merge(responses, on='SourceSequence')

    # some of the MHC scores come back as all NaN so drop them
    result = result.dropna(axis=1, how='all')

    drop_fields = ('seq_num', 'method')
    for field in drop_fields:
        if field in result:
            result = result.drop(field, axis = 1)
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
        "method" : self._method,
        "sequence_text" : sequence,
        "allele" : self._alleles,
      }
      return params
