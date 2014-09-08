
from subprocess import Popen, CalledProcessError, check_output, PIPE
import tempfile
import os
import logging
import time 

import numpy as np
import pandas as pd 
from epitopes.mutate import gene_mutation_description

from mhc_common import normalize_hla_allele_name

def create_input_fasta_file(df, mutation_window_size = None):
    """
    Turn peptide entries from a dataframe into a FASTA file. If mutation_window_size is 
    an integer >0 then only use subsequence around mutated residues. 

    Return the name of closed file which has to be manually deleted, and a dictionary from 
    FASTA IDs to peptide records. 
    """
    input_file = tempfile.NamedTemporaryFile("w", prefix="peptide", delete=False)
    
    peptide_entries = {}
    records = df.to_records()
    n_records = len(records)
    # create input file for all peptide sequences and also
    # put the entries into a dictionary so we can read out the results later
    for i, mutation_entry in enumerate(records):
        seq =  mutation_entry['SourceSequence']
        if mutation_window_size:
            start = max(0, mutation_entry.MutationStart - mutation_window_size)
            stop = min(len(seq), mutation_entry.MutationEnd + mutation_window_size)
            seq = seq[start:stop]
        identifier = "%s_%s" % (i, mutation_entry['Gene'][:5])
        peptide_entries[identifier] = mutation_entry

        input_file.write(">%s\n" % identifier)
        input_file.write(seq)
        # newline unless at end of file
        if n_records > i + 1:
            input_file.write("\n")
    input_file.close()  
    return input_file.name, peptide_entries

def bad_binding_score(x):
    return x < 0 or np.isnan(x) or np.isinf(x)

def build_output_rows(lines, peptide_entries, mutation_window_size = None):
    """
    First line of XLS file format has HLA alleles
    and second line has fields like:
        ['Pos', 'Peptide', 'ID', 
         '1-log50k', 'nM', 'Rank', 
         '1-log50k', 'nM', 'Rank', 
         '1-log50k', 'nM', 'Rank', 
         ...'Ave', 'NB']
    """
    lines = [line.split("\t") for line in lines if len(line) > 0]
    alleles = [x for x in lines[0] if len(x) > 0]
    # skip alleles and column headers
    lines = lines[2:]
    results = []
    for line in lines[2:]:
        pos = int(line[0])
        epitope = line[1]
        identifier = line[2]
        assert identifier in peptide_entries, "Bad identifier %s, epitopes = %s" % (identifier, epitopes.head())
        mutation_entry = peptide_entries[identifier]

        if mutation_window_size:
            # if we clipped parts of the amino acid sequence which don't overlap mutations
            # then we have to offset epitope positions by however much was removed from the 
            # beginning of the sequence
            original_start = max(0, mutation_entry.MutationStart - mutation_window_size)
            pos += original_start

        for i, allele in enumerate(alleles):
            
            # we start at an offset of 3 to skip the allele-invariant pos, epitope, identifier columns
            # each allele has three columns: log IC50, IC50, rank 
            log_ic50 = float(line[3+3*i])
            ic50 = float(line[3+3*i+1])
            rank = float(line[3+3*i+2])
            
            # if we have a bad IC50 score we might still get a salvageable 
            # log of the score. Strangely, this is necessary sometimes! 
            if bad_binding_score(ic50):
                ic50 = 50000 ** (-log_ic50 + 1)

            if bad_binding_score(ic50): 
                logging.warn("Invalid IC50 value %0.4f for %s w/ allele %s" % (ic50, epitope, allele))
                continue 
            elif bad_binding_score(rank) or rank > 100:
                logging.warn("Invalid percentile rank %s for %s w/ allele %s" % (rank, epitope, allele))
                continue 

            # keep track of original genetic variant that gave rise to this epitope
            new_row = {}
            # fields shared by all epitopes from this sequence 
            new_row['SourceSequence'] = mutation_entry.SourceSequence
            new_row['MutationStart'] = mutation_entry.MutationStart
            new_row['MutationEnd'] = mutation_entry.MutationEnd
            new_row['GeneInfo'] = mutation_entry.GeneInfo
            new_row['Gene'] = mutation_entry.Gene
            new_row["GeneMutationInfo"] = mutation_entry.GeneMutationInfo
            new_row['PeptideMutationInfo'] = mutation_entry.PeptideMutationInfo
            new_row['TranscriptId'] = mutation_entry.TranscriptId

            # fields specific to this epitope 
            new_row['Allele'] = normalize_hla_allele_name(allele)
            new_row['EpitopeStart'] = pos 
            new_row['EpitopeEnd'] = pos + len(epitope)
            new_row['Epitope'] = epitope 
            new_row['MHC_IC50'] = ic50
            new_row['MHC_PercentileRank'] = rank 
            results.append(new_row) 
    return results 

class PanBindingPredictor(object):

    def __init__(self, hla_alleles, netmhc_command = "netMHCpan"):
        self.netmhc_command = netmhc_command
        
        try:
            valid_alleles_str = check_output([self.netmhc_command, "-listMHC"])
            assert len(valid_alleles_str) > 0, "%s returned empty allele list" % self.self.netmhc_command
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
            # for some reason netMHCpan drop the "*" in names such as "HLA-A*03:01" becomes "HLA-A03:01"
            if  allele.replace("*", "") not in valid_alleles:
                print "Skipping %s (not available in NetMHCpan)" % allele
            else:
                self.alleles.append(allele)
        # don't run the MHC predictor twice for homozygous alleles,
        # only run it for unique alleles
        self.alleles = set(self.alleles)


    def predict(self, df, mutation_window_size = None):
        """
        Given a dataframe of mutated amino acid sequences, run each sequence through NetMHCpan. 
        If mutation_window_size is not None then only make predictions for that number residues
        away from mutations. 

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

        input_filename, peptide_entries = \
            create_input_fasta_file(df, mutation_window_size = mutation_window_size)

        output_files = {}
        processes = {}
        process_commands = {}

        def cleanup():
            """
            Cleanup either when finished or if an exception gets raised by 
            killing all the running netMHCpan processes and deleting the input 
            and output files
            """
            for p in processes.itervalues():
                try:
                    p.kill()
                except:
                    pass 
            
            for output_file in output_files.itervalues():
                try:
                    output_file.close()
                    os.remove(output_file.name)
                except:
                    pass 

            os.remove(input_filename)


        alleles_str = ",".join(allele.replace("*", "") for allele in self.alleles)
        output_file = \
            tempfile.NamedTemporaryFile(
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
        try: 
            start_time = time.time()
            with open(os.devnull, 'w') as devnull:
                process = Popen(command, stdout = devnull)
            
            ret_code = process.wait()
            
            if ret_code:
                logging.info("netMHCpan finished with return code %s", ret_code)
                raise CalledProcessError(ret_code, process_commands[allele])
            else:
                elapsed_time = time.time() - start_time
                logging.info("netMHCpan took %0.4f seconds", elapsed_time)
                lines = output_file.read().split("\n")
                results = build_output_rows(lines, peptide_entries, mutation_window_size = mutation_window_size)
        except:
            cleanup()
            raise 
        cleanup()
        assert len(results) > 0, "No epitopes from netMHCpan"
        return pd.DataFrame.from_records(results)
