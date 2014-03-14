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

from pipeline import PipelineElement
import urllib2, urllib
import pandas as pd
from StringIO import StringIO

class IEDBMHCBinding(PipelineElement):

  def __init__(self, alleles=[], name="IEDB-MHC-Binding", method='recommended', lengths = [9,10,11], url='http://tools.iedb.org/tools_api/mhci/'):
    self.name = name
    self._method = method
    self._lengths = lengths
    self._url = url
    self._alleles = ",".join(alleles)

  def _get_iedb_request_params(self, sequence):
    params = {
        "method" : self._method,
        "length" : ",".join(str(l) for l in self._lengths),
        "sequence_text" : sequence,
        "allele" : self._alleles,
    }
    return params

  def query_iedb(self, sequence):
    request_values = self._get_iedb_request_params(sequence)
    print "Calling iedb with", sequence, self._alleles
    try:
      data = urllib.urlencode(request_values)
      req = urllib2.Request(self._url, data)
      response = urllib2.urlopen(req).read()

      return pd.read_csv(StringIO(response), sep='\t', na_values=['-'])
    except:
      print "Connection error: Failed on sequence", sequence
      return pd.DataFrame()

  def apply(self,data):
    responses = {}
    for epitope in data.Epitope:
       responses[epitope] = self.query_iedb(epitope)
    responses = pd.concat(responses).reset_index(0)
    responses.rename(columns={'level_0':'Epitope'}, inplace=True)
    return data.merge(responses, on='Epitope')

class IEDBMHC1Binding(IEDBMHCBinding):
    def __init__(self, name = 'IEDB-MHC1-Binding', url='http://tools.iedb.org/tools_api/mhci/', alleles=[]):
      super(IEDBMHC1Binding, self).__init__(name = name, url = url, alleles = alleles)


class IEDBMHC2Binding(IEDBMHCBinding):
    def __init__(self, name = 'IEDB-MHC2-Binding', url='http://tools.iedb.org/tools_api/mhcii/', alleles=[]):
      super(IEDBMHC2Binding, self).__init__(name = name, url = url, alleles = alleles)

    def _get_iedb_request_params(self, sequence):
      params = {
        "method" : self._method,
        "sequence_text" : sequence,
        "allele" : self._alleles,
      }
      return params



if __name__ == '__main__':
  iedb = IEDBMHC1Binding(alleles=['HLA-C*12:03', 'HLA-C*12:02'])
  print iedb.query_iedb("APHHSGVYPVNVQLYEAWKKV")

  iedb = IEDBMHC2Binding(alleles=['HLA-DRB1*01:01','H2-IAb'])
  print iedb.query_iedb("APHHSGVYPVNVQLYEAWKKV")
