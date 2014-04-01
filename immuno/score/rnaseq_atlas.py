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


from epitopes.download import fetch_data
import pandas as pd

from common import memoize
from entrez import entrez_hugo_dataframe

ATLAS_URL = "http://medicalgenomics.org/rna_seq_atlas/download?download_revision1=1"
ATLAS_FILENAME = "RNASeq_Atlas.txt"

TISSUE_COLUMNS = [
    'adipose',
    'colon',
    'heart',
    'hypothalamus',
    'kidney',
    'liver',
    'lung',
    'ovary',
    'skeletalmuscle',
    'spleen',
    'testes'
]

@memoize
def load_dataframe():
    """
    Return a dataframe mapping Entrez gene IDs to
    tissue-specific RPKMs. Index columns are
    'entrez_gene_id', 'ensembl_gene_id', and 'hgnc_symbol'.

    Tissue columns are:
        - adipose
        - colon
        - heart
        - hypothalamus
        - kidney
        - liver
        - lung
        - ovary
        - skeletalmuscle
        - spleen
        - testes

    """
    path = fetch_data(ATLAS_FILENAME, ATLAS_URL)
    return pd.read_csv(path, sep='\t', header=0)

def hugo_to_rpkm():
    df = load_dataframe()
    df = df[['hgnc_symbol'] + TISSUE_COLUMNS]
    df = df.rename(columns = {'hgnc_symbol' : 'Hugo'})
    # multiple transcripts result in multiple entries
    # for the same Hugo gene ID
    # average together their (similar values)
    # and then restore the 'Hugo' column with a call to 'reset_index'
    return df.groupby("Hugo").mean().reset_index()

def entrez_to_rpkm():
    df = load_dataframe()
    df = df[['entrez_gene_id'] + TISSUE_COLUMNS]
    df = df.rename(columns = {'entrez_gene_id' : 'Entrez'})
    # multiple transcripts result in multiple entries
    # for the same Entrez gene ID
    # average together their (similar values)
    # and then restore the 'Entrez' column with a call to 'reset_index'
    return df.groupby("Entrez").mean().reset_index()

def hugo_to_rank(group_rank_method='average'):
    df = hugo_to_rpkm()
    values = df[TISSUE_COLUMNS]
    ranks = values.rank(method=group_rank_method)  / len(values)
    df[TISSUE_COLUMNS] = ranks
    return df

def entrez_to_rank(group_rank_method='average'):
    df = entrez_to_rpkm()
    values = df[TISSUE_COLUMNS]
    ranks = values.rank(method=group_rank_method) / len(values)
    df[TISSUE_COLUMNS] = ranks
    return df

