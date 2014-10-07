
import tempfile
import os
import logging
import time

import numpy as np
import pandas as pd

from cleanup_context import CleanupFiles
from process_helpers import run_command
from mhc_common import normalize_hla_allele_name
from mhc_formats import create_input_fasta_file, parse_xls_file


class PanBindingPredictor(object):

    def __init__(
            self,
            hla_alleles,
            netmhc_command = "netMHCpan"):
        self.netmhc_command = netmhc_command

        try:
            run_command([self.netmhc_command])
        except:
            assert False, "Failed to run %s" % self.netmhc_command

        try:
            valid_alleles_str = check_output([self.netmhc_command, "-listMHC"])
            assert len(valid_alleles_str) > 0, \
                "%s returned empty allele list" % self.self.netmhc_command
            valid_alleles = set([])
            for line in valid_alleles_str.split("\n"):
                if not line.startswith("#"):
                    valid_alleles.add(line)
        except:
            logging.warning("Failed to run %s -listMHC", self.netmhc_command)
            valid_alleles = None

        self.alleles = []
        for allele in hla_alleles:
            allele = normalize_hla_allele_name(allele.strip().upper())
            # for some reason netMHCpan drop the "*" in names
            # such as "HLA-A*03:01" becomes "HLA-A03:01"
            if valid_alleles and allele.replace("*", "") not in valid_alleles:
                print "Skipping %s (not available in NetMHCpan)" % allele
            else:
                self.alleles.append(allele)
        # don't run the MHC predictor twice for homozygous alleles,
        # only run it for unique alleles
        self.alleles = set(self.alleles)


    def predict(self, df, mutation_window_size = None):
        """
        Given a dataframe of mutated amino acid sequences, run each sequence
        through NetMHCpan.
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

        alleles_str = \
            ",".join(allele.replace("*", "") for allele in self.alleles)
        output_file =  tempfile.NamedTemporaryFile(
                "r+",
                prefix="netMHCpan_output",
                delete=False)
        command = [
            self.netmhc_command,
                "-xls",
                "-xlsfile", output_file.name,
                 "-l", "9",
                  "-f", input_filename,
                  "-a", alleles_str]
        print " ".join(command)

        with CleanupFiles(
                filenames = [input_filename],
                files = [output_file]):
            run_command(command)
            results = parse_xls_file(
                output_file.read(), peptide_entries,
                mutation_window_size=mutation_window_size
            )

        assert len(results) > 0, "No epitopes from netMHCpan"
        return pd.DataFrame.from_records(results)
