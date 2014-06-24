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

import pandas as pd

def parse_string(string):
    """
    Normalize a string by stripping out whitespace and converting to uppercase.
    In case the string contains annotations for a mutated region (denoted by
    characters between underscores (i.e. "QYY_LS_YY"), parse out the start/stop
    of the mutated region. If there are no underscores, the start is 0 and the
    stop is the length of the string.
    """

    full_peptide = string.upper().strip()
    # allow the user to specify mutated region of the amino acid
    # string QLSQ_Y_QQ (the full peptide is QLSQYQQ and Y is mutated)
    parts = full_peptide.split("_")
    if len(parts) == 1:
        full_peptide = parts[0]
        start = 0
        stop = len(full_peptide)
    elif len(parts) == 2:
        full_peptide = parts[0] + parts[1]
        start = len(parts[0])
        stop = len(full_peptide)
    else:
        assert len(parts) == 3, \
            "Can't parse peptide string %s" % full_peptide
        full_peptide = parts[0] + parts[1] + parts[2]
        start = len(parts[0])
        stop = start + len(parts[1])
    return full_peptide, start, stop


def load_strings(strings, source = 'commandline'):
    # allow multiple strings to be specified in comma-separated list

    rows = []
    for string in strings:
        full_peptide, start, stop = parse_string(string)
        row = {}
        row['SourceSequence'] = full_peptide
        row['MutationStart'] = start
        row['MutationEnd'] = stop
        row['info'] = source
        row['MutationInfo'] = "-"
        row['stable_id_transcript'] = "-"

        rows.append(row)
    return pd.DataFrame.from_records(rows,
        columns = (
            'SourceSequence',
            'MutationStart',
            'MutationEnd',
            'info',
            'MutationInfo',
            'stable_id_transcript',
        ))

def load_comma_string(s, source = 'literal'):
    return load_strings(s.split(","), source = source)