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
import os
from os.path import splitext, abspath, join
import time

import Bio.Seq

def peptide_substrings(full_peptide, window_length):
    n = len(full_peptide)
    return [full_peptide[i:i+window_length]
            for i in xrange(n + 1 - window_length)]

def normalize_chromosome_name(c):
    """
    standardize chromosome names by getting rid of any "chr"
    prefixes such as "chr17" and always naming the mitochondrial DNA "M"
    """

    # in case a chromosome is being represented numerically
    # convert it to a string (i.e. 1 -> "1")
    if  isinstance(c, (long, int)):
        c = str(c)
    # remove "chr" prefix (i.e. "chr1" -> "1")
    c = c.lower().replace('chr', '').upper()
    # depending on the reference, the mitochondrial
    # DNA is represented as either "M" or "MT",
    # standardize on "M"
    return "M" if (c == "MT") else c

VALID_AMINO_ACIDS = set([
    'A', 'R', 'N',
    'D', 'C', 'E',
    'Q', 'G', 'H',
    'I', 'L', 'K',
    'M', 'F', 'P',
    'S', 'T', 'W',
    'Y', 'V'
])

def is_valid_peptide(pep):
	return isinstance(pep, (str, Bio.Seq.Seq)) and \
        len(pep) > 0 and \
        all(residue in VALID_AMINO_ACIDS for residue in pep)

def init_logging(quiet = False):
    log_level = logging.WARNING if quiet else logging.DEBUG

    logging.basicConfig(
        format=\
            "[%(levelname)s %(filename)s:%(lineno)d %(funcName)s] %(message)s",
        level=log_level
    )

def splitext_permissive(path, ignored_exts):
    """
    Runs splitext on the path. However, if ext is in ignored_exts, it re-runs
    splitext on the base: until finding an ext that is not in ignored_exts.

    Note that these exts include the dot: ".txt", ".bam", etc.
    """
    if "" in ignored_exts:
        raise ValueError("ignored_exts cannot contain the empty string")
    base, ext = splitext(path)
    if ext in ignored_exts:
        return splitext_permissive(base, ignored_exts)
    return (base, ext)


def find_paths(filename_string = "", directory_string = "", extensions = None):
    """
    Parse input comma separated list of files and comma separated list
    of directories, collecting all of their full file paths.
    """
    paths = []
    if filename_string:
        for filename in filename_string.split(","):
            paths.append(abspath(filename.strip()))
    if directory_string:
        for dirpath in directory_string.split(","):
            for filename in os.listdir(dirpath):
                path = join(dirpath, filename)
                paths.append(path)
    if extensions is not None:
        paths = [
            p for p in paths
            if any(p.endswith(ext) for ext in extensions)
        ]
    return paths

