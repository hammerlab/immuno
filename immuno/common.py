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

class AsyncProcess(object):
    """
    A thin wrapper around Popen which starts a process asynchronously,
    suppresses stdout printing, and raises an exception if the return code 
    of wait() isn't 0
    """

    def __init__(self, args, suppress_stdout = True, suppress_stderr = False):
        assert len(args) > 0
        self.cmd = args[0]
        with open(os.devnull, 'w') as devnull:
            stderr = devnull if suppress_stderr else None
            stdout = devnull if suppress_stdout else None
            self.process = Popen(args, stdout = stdout, stderr = stderr)

    def wait(self):
        ret_code = self.process.wait()
        logging.info(
            "%s finished with return code %s", 
            self.cmd,
            ret_code)
        if ret_code:
            raise CalledProcessError(ret_code, self.cmd)
        return ret_code

def run_command(args, **kwargs):
    """
    Given a list whose first element is a command name, followed by arguments, 
    execute it and show timing info. 
    """
    assert len(args) > 0
    cmd = args[0]
    start_time = time.time()
    process = AsyncProcess(args, **kwargs)
    process.wait()
    elapsed_time = time.time() - start_time
    logging.info("%s took %0.4f seconds", cmd, elapsed_time)

def run_multiple_commands(multiple_args_lists, print_commands=True, **kwargs):
    assert len(multiple_args_lists) > 0
    assert all(len(args) > 0 for args in multiple_args_lists)
    start_time = time.time()
    command_names = [args[0] for args in multiple_args_lists]
    processes = []
    for args in multiple_args_lists:
        if print_commands:
            print " ".join(args)
        p = AsyncProcess(args, **kwargs)
        processes.append(p)

    for p in processes:
        p.wait()

    elapsed_time = time.time() - start_time
    logging.info("Ran %d commands (%s) in %0.4f seconds",
        len(multiple_args_lists),
        ",".join(command_names),
        elapsed_time
    )

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
                logging.warning("Failed to close %s", f) 

            try:
                os.remove(f.name)
            except:
                logging.warning("Failed to remove file %s", f) 

        for name in self.filenames:
            logging.info("Cleaning up %s", name)
            try:
                os.remove(name)
            except:
                logging.warning("Failed to remove filename %s", name) 