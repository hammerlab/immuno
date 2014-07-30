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

import logging

import pandas as pd
import numpy as np

from ensembl import annotation, gene_names 
from ensembl.transcript_variant import peptide_from_transcript_variant

def load_vcf(input_filename, drop_low_quality = True):
    """
    Parameters
    --------

    input_filename : str
        Path to VCF file

    drop_low_quality : bool, optional
        Keep variants whose 'QUAL' columns doesn't say 'PASS' or '.'

    Returns Dataframe with fields
            - 'chr'
            - 'pos'
            - 'id'
            - 'ref'
            - 'alt'
            - 'qual'
            - 'filter'
            - 'info'
    """
    with open(input_filename) as fd:
        lines_to_skip = 0
        while next(fd).startswith('#'):
            lines_to_skip += 1

    # first 8 columns of a VCF file are required to be:
    #   chr    : chromosome (i.e., '20', 'chr20', 'MT')
    #   pos    : where's the variant on the chromsome?
    #   id     : dbSnp identifier of variant, if available (i.e. 'rs11449')
    #   ref    : reference letter(s)
    #   alt    : alternate letters(s) of variant
    #   qual   : phred-scaled quality score
    #   filter : "PASS" if variant passes all filters
    #   info   : optional info, comma separated key=value list
    header = ['chr', 'pos', 'id', 'ref', 'alt','qual', 'filter', 'info']

    df = pd.read_csv(
            input_filename,
            sep='\t',
            skiprows=lines_to_skip,
            header = None,
            names=header,
            usecols=header,
            dtype={'pos' : np.int32, 'chr':str})

    
    # drop variants marked as low quality
    if drop_low_quality:
        qual = df['qual']
        mask = (qual == 'PASS') | (qual == '.')
        df = df[mask]

    return df 
