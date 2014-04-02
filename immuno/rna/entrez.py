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


import pandas as pd
import epitopes

from common import memoize

ENTREZ_HUGO_URL = \
"http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=gd_pub_eg_id&status=Approved&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&hgnc_dbtag=on&submit=submit"

ENTREZ_HUGO_FILENAME = "entrez_huge_gene_mapping.txt"

@memoize
def entrez_hugo_dataframe():
    """
    Download the mappings between Entrez and Hugo gene IDs
    and load them in a dataframe.
    """
    path = epitopes.download.fetch_data(
        ENTREZ_HUGO_FILENAME, ENTREZ_HUGO_URL, subdir='immuno')
    df = pd.read_csv(
            path,
            sep='\t',
            names=('Hugo', 'Entrez'),
            header=0)
    return df

@memoize
def entrez_hugo_mapping():
    """
    Returns a dictionary mapping HUGO gene names
    to Entrez ID numbers.
    """
    df = entrez_hugo_dataframe()
    return dict(zip(df.Entrez, df.Hugo))

@memoize
def hugo_entrez_mapping():
    """
    Returns a dictionary mapping Entrez gene ID numbers
    to HUGO gene names.
    """
    df = entrez_hugo_dataframe()
    return dict(zip(df.Hugo, df.Entrez))

