
import subprocess 
import tempfile
import os
import logging

import numpy as np
import pandas as pd 
from epitopes.mutate import gene_mutation_description

from mhc_common import normalize_hla_allele_name

def create_input_fasta_file(df):
    """
    Turn peptide entries from a dataframe into a FASTA file. 
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

def build_output_rows(epitopes_df, peptide_entries, allele):
    results = []
    for identifier, group in epitopes_df.groupby("ID"):
        assert identifier in peptide_entries, "Bad identifier %s, epitopes = %s" % (identifier, epitopes.head())
        mutation_entry = peptide_entries[identifier]
        # columns: Pos    Peptide   ID  1-log50k          nM  Rank
        for epitope_row in group.to_records():
            pos = epitope_row['Pos']
            epitope = epitope_row['Peptide']
            ic50 = epitope_row['nM']
            rank = epitope_row['Rank']
            # if we have a bad IC50 score we might still get a salvageable 
            # log of the score. Strangely, this is necessary sometimes! 
            if bad_binding_score(ic50):
                log_ic50 = epitope_row['1-log50k']
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
            new_row['Allele'] = allele
            new_row['EpitopeStart'] = pos 
            new_row['EpitopeEnd'] = pos + len(epitope)
            new_row['Epitope'] = epitope 
            new_row['MHC_IC50'] = ic50
            new_row['MHC_PercentileRank'] = rank 
            results.append(new_row) 
    return results 

class PanBindingPredictor(object):
    def __init__(self, hla_alleles):
        valid_alleles_str = subprocess.check_output(["netMHCpan", "-listMHC"])
        valid_alleles = set([])
        for line in valid_alleles_str.split("\n"):
            if not line.startswith("#"):
                valid_alleles.add(line)
        self.alleles = []
        for allele in hla_alleles:
            allele = normalize_hla_allele_name(allele.strip().upper())
            # for some reason netMHCpan drop the "*" in names such as "HLA-A*03:01" becomes "HLA-A03:01"
            if  allele.replace("*", "") not in valid_alleles:
                print "Skipping %s (not available in NetMHCpan)" % allele
            else:
                self.alleles.append(allele)


    def predict(self, df):
        """
        Given a dataframe of mutated amino acid sequences, run each sequence through NetMHCpan. 

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


        input_filename, peptide_entries = create_input_fasta_file(df)

        results = []
        for allele in self.alleles:
            output_file = tempfile.NamedTemporaryFile("w", prefix="netMHCpan_output", delete=False)
            command = ["netMHCpan",  "-xls", "-xlsfile", output_file.name, "-l", "9", "-f", input_filename, "-a", allele.replace("*", "")]
            print "Calling netMHCpan for %s" % allele 
            print " ".join(command)
            subprocess.check_output(command)
            epitopes_df = pd.read_csv(output_file.name, sep='\t', skiprows = 1)
            output_file.close()
            os.remove(output_file.name)
            results.extend(build_output_rows(epitopes_df, peptide_entries, allele))
        os.remove(input_filename)
        assert len(results) > 0, "No epitopes from netMHCpan"
        return pd.DataFrame.from_records(results)
