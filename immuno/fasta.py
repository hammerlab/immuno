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
from Bio import SeqIO

from strings import parse_string

def load_fasta(fasta_filename, peptide_length):
    seqs = []
    ids = []
    names = []
    descriptions = []
    starts = []
    stops = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq, start, stop = parse_string(str(record.seq))
        seqs.append(record.seq)
        ids.append(record.id)
        names.append(record.name)
        descriptions.append(record.description)

    assert len(ids) == len(seqs) == len(names) == len(descriptions) ==\
        len(starts) == len(stops)
    return pd.DataFrame({
        'SourceSequence': seqs,
        'MutationStart' : starts,
        'MutationStop' : stops,
        'Id' : ids,
        'Name' : names,
        'Description' : descriptions,
        'info': [fasta_filename] * len(seqs),
    })
