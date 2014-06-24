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

from pipeline import PipelineElement

def seq_to_str(obj):
    """
    Given a sequence convert it to a comma separated string. 
    If, however, the argument is a single object, return its string representation.
    """
    if isinstance(obj, (unicode, str)):
        return obj
    elif isinstance(obj, (list, tuple)):
        return  ",".join([str(x) for x in obj])
    else:
        return str(obj)

def convert_str(obj):
    """
    Given a string, convert it to an int or float if possible.
    """
    if obj is None: 
        return obj 
    try:
        try:
            return int(obj)
        except:
            return float(obj)
    except:
        return str(obj)

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
        method=['consensus'],
        lengths = [9],
        url='http://tools.iedb.org/tools_api/mhci/'):
    self.name = name
    self._method = method
    self._lengths = lengths
    self._url = url
    self._alleles = alleles

  def _get_iedb_request_params(self, sequence):
    
    params = {
        "method" : seq_to_str(self._method),
        "length" : seq_to_str(self._lengths),
        "sequence_text" : sequence,
        "allele" : seq_to_str(self._alleles)
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
            response.rename(
                columns={
                    'peptide': 'Epitope',
                    'length' : 'EpitopeLength',
                    'start' : 'EpitopeStart',
                    'end' : 'EpitopeEnd',
                    'allele' : 'Allele', 
                },
                inplace=True)
            responses[peptide] = response 
        else:
            logging.info(
                "Skipping binding for peptide %s, already queried",
                peptide)

   
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
    responses['MHC_PercentileRank'] = responses['ann_rank']

    assert 'ann_ic50' in responses, responses.head()
    responses['MHC_IC50'] = responses['ann_ic50']
     
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
