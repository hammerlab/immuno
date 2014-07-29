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

import appdirs
import logging


def peptide_substrings(full_peptide, window_length):
    n = len(full_peptide)
    return [full_peptide[i:i+window_length]
            for i in xrange(n + 1 - window_length)]

def normalize_chromosome_name(c):
	"""
	standardize chromosome names by getting rid of any "chr" prefixes such as "chr17" and
	always naming the mitochondrial DNA "M"
	"""
	c = c.lower().replace('chr', '').upper()
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
	return all(residue in VALID_AMINO_ACIDS for residue in pep)

def init_logging(quiet = False):
    log_level = logging.WARNING if quiet else logging.DEBUG
       
    logging.basicConfig(
        format="[%(levelname)s %(filename)s:%(lineno)d %(funcName)s] %(message)s",
        level=log_level
    )