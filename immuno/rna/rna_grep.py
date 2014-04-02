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

import argparse
import pandas as pd
import numpy as np
import entrez
import rnaseq_atlas
from load_rna_data import load_rsem


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filename",
        "-f",
        required=True,
        help="""File with two tab separated columns: Entrez gene IDs and RSEM expression values""")
    parser.add_argument(
        "pattern",
        help="Gene name, part of a gene name, or regular expression ")

    args = parser.parse_args()
    df = load_rsem(args.filename, translate_entrez_ids = True)
    mask = df.Hugo.str.contains(args.pattern)
    subset = df[mask]
    if len(subset) == 0:
        print "Pattern %s not found" % args.pattern
    else:
        print subset