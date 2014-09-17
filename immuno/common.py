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
from os.path import splitext
from subprocess import Popen, CalledProcessError
import time 

import appdirs

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


def run_command(args):
    """
    Given a list whose first element is a command name, followed by arguments, 
    execute it and show timing info. 
    """
    cmd = args[0]
    start_time = time.time()
    with open(os.devnull, 'w') as devnull:
        process = Popen(args, stdout = devnull)

    ret_code = process.wait()

    if ret_code:
        logging.info(
            "%s finished with return code %s", 
            cmd,
            ret_code)
        raise CalledProcessError(ret_code, cmd)
    else:
        elapsed_time = time.time() - start_time
        logging.info("%s took %0.4f seconds", cmd, elapsed_time)


class CleanupFiles(object):
    """
    Context manager that deletes a set of files at the end of a block or
    if an exception gets raised 
    """
    def __init__(self, files = [], filenames = [], dictionaries = []):
        self.files = files 
        self.filenames = []
        self.dictionaries = []

    def __enter__(self):
        pass 

    def __exit__(self, type, value, traceback):
        files = self.files 
        
        # extend files with all values from all dictionaries we're tracking
        for d in self.dictionaries:
            files.extend(d.values())
        
        for f in files:
            logging.info("Cleaning up %s", f)
            try:
                f.close()
            except:
                pass 

            try:
                os.remove(f.name)
            except:
                pass 

        for name in self.filenames:
            logging.info("Cleaning up %s", name)
            try:
                os.remove(name)
            except:
                pass 