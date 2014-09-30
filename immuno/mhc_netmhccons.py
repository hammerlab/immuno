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
import subprocess
import tempfile
import time 

import numpy as np
import pandas as pd 
from epitopes.mutate import gene_mutation_description

from common import run_multiple_commands, CleanupFiles
from mhc_common import normalize_hla_allele_name
from mhc_formats import create_input_fasta_file, parse_xls_file


class ConsensusBindingPredictor(object):

    def __init__(
            self, 
            hla_alleles, 
            netmhc_command = "netMHCcons"):
        self.netmhc_command = netmhc_command
        
        try:
            subprocess.check_output([self.netmhc_command],
                stderr=subprocess.STDOUT)
        except:
            assert False, "Failed to run %s" % self.netmhc_command

        # normalize alleles and keep only unique names
        normalized_alleles = {
            normalize_hla_allele_name(allele.strip().upper()).replace("*", "")
            for allele in hla_alleles
        }

        self.alleles = []

        # try running "netMHCcons -a" with each allele name
        # and check if it gives you back a "wrong format" error 
        for allele in normalized_alleles:
            try:
                subprocess.check_output(
                    [self.netmhc_command, '-a', allele],
                    stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError, e:
                if "allele" in e.output and "wrong format" in e.output:
                    logging.warning(
                        "Allele %s not recognized by NetMHCcons", allele)
                    continue
            except:
                pass 
            logging.info("Normalize HLA allele %s", allele)
            self.alleles.append(allele)


    def predict(self, df, mutation_window_size = None):
        """
        Given a dataframe of mutated amino acid sequences, run each sequence 
        through NetMHCcons. 
        If mutation_window_size is not None then only make predictions for that
        number residues away from mutations. 

        Expects the input DataFrame to have the following fields: 
            - SourceSequence
            - MutationStart
            - MutationEnd
            - GeneInfo
            - Gene
            - GeneMutationInfo
            - PeptideMutationInfo
            - TranscriptId
        """

        input_filename, peptide_entries = create_input_fasta_file(
            df,
            mutation_window_size=mutation_window_size
        )

        output_files = {}
           
        commands_list = []
        for i, allele in enumerate(self.alleles):

            output_file = tempfile.NamedTemporaryFile(
                    "r+", 
                    prefix="netMHCcons_output_%d" % i,
                    delete=False)
            output_files[allele] = output_file
            command = [
                self.netmhc_command,  
                    "-xls",
                    "-xlsfile", output_file.name,
                    "-length", "9",
                    "-f", input_filename,
                    "-a", allele]
            commands_list.append(command)
            
        results = []
        # Cleanup either when finished or if an exception gets raised by 
        # deleting the input and output files
        with CleanupFiles(
                filenames = [input_filename], 
                dictionaries = [output_files]):
            run_multiple_commands(commands_list, print_commands = True)
            results.extend(
                    parse_xls_file(
                        output_file.read(),
                        peptide_entries, 
                        mutation_window_size=mutation_window_size))
        assert len(results) > 0, "No epitopes from netMHCcons"
        return pd.DataFrame.from_records(results)
