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
import numpy as np

import entrez
import rnaseq_atlas

def load_rsem(filename,
        translate_entrez_ids = True,
        group_transcripts = True):
    if translate_entrez_ids:
        entrez_df = pd.read_csv(filename, sep='\t', names=('Entrez', 'RSEM'))

        # mapping from Entrez IDs to Hugo
        hugo_mapping = entrez.entrez_hugo_dataframe()
        # add Hugo column
        merged_df = entrez_df.merge(hugo_mapping, on='Entrez')
        # keep only the Hugo IDs and RSEM
        hugo_df = merged_df[['Hugo', 'RSEM']]
    else:
        hugo_df = pd.read_csv(args.sample, sep='\t', names=('Hugo', 'RSEM'))
    if group_transcripts:
        hugo_df = hugo_df.groupby("Hugo").max().reset_index()
    assert len(hugo_df) > 0, "No data found in %s" % filename
    return hugo_df