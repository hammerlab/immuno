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


import tempfile

def create_input_fasta_file(df, mutation_window_size = None):
    """
    Turn peptide entries from a dataframe into a FASTA file. 
    If mutation_window_size is an integer >0 then only use subsequence
    around mutated residues. 

    Return the name of closed file which has to be manually deleted, 
    and a dictionary from FASTA IDs to peptide records. 
    """
    input_file = tempfile.NamedTemporaryFile(
        "w", prefix="peptide", delete=False)
    
    peptide_entries = {}
    records = df.to_records()
    n_records = len(records)
    # create input file for all peptide sequences and also
    # put the entries into a dictionary so we can read out the results later
    for i, mutation_entry in enumerate(records):
        seq =  mutation_entry['SourceSequence']
        if mutation_window_size:
            start = max(
                0, 
                mutation_entry.MutationStart - mutation_window_size)
            stop = min(
                len(seq), 
                mutation_entry.MutationEnd + mutation_window_size)
            seq = seq[start:stop]
        identifier = "%s_%s" % (i, mutation_entry['Gene'][:5])
        peptide_entries[identifier] = mutation_entry

        input_file.write(">%s\n" % identifier)
        input_file.write(seq)
        # newline unless at end of file
        if  i + 1 < n_records:
            input_file.write("\n")
    input_file.close()  
    return input_file.name, peptide_entries

def invalid_binding_score(x):
    return x < 0 or np.isnan(x) or np.isinf(x)

def parse_xls_file(contents, peptide_entries, mutation_window_size = None):
    """
    XLS is a wacky output format used by NetMHCpan and NetMHCcons
    for peptide binding predictions. 

    First line of XLS file format has HLA alleles
    and second line has fields like:
        ['Pos', 'Peptide', 'ID', 
         '1-log50k', 'nM', 'Rank', 
         '1-log50k', 'nM', 'Rank', 
         '1-log50k', 'nM', 'Rank', 
         ...'Ave', 'NB']
    """
    lines = [line.split("\t") 
             for line in contents.split("\n")
             if len(line) > 0]
    # top line of XLS file has alleles
    alleles = [x for x in lines[0] if len(x) > 0]
    # skip alleles and column headers
    lines = lines[2:]
    results = []
    for line in lines:
        pos = int(line[0])
        epitope = line[1]
        identifier = line[2]
        
        assert identifier in peptide_entries, \
            "Bad identifier %s, epitopes = %s" % (identifier, epitopes.head())
        mutation_entry = peptide_entries[identifier]

        if mutation_window_size:
            # if we clipped parts of the amino acid sequence which don't
            # overlap mutations then we have to offset epitope positions by
            # however much was removed from the beginning of the sequence    
            original_start = max(
                0, 
                mutation_entry.MutationStart - mutation_window_size
            )
            pos += original_start

        for i, allele in enumerate(alleles):
            
            # we start at an offset of 3 to skip the allele-invariant
            # pos, epitope, identifier columns
            # each allele has three columns: log IC50, IC50, rank 
            log_ic50 = float(line[3+3*i])
            ic50 = float(line[3+3*i+1])
            rank = float(line[3+3*i+2])
            
            # if we have a bad IC50 score we might still get a salvageable 
            # log of the score. Strangely, this is necessary sometimes! 
            if invalid_binding_score(ic50):
                ic50 = 50000 ** (-log_ic50 + 1)

            if invalid_binding_score(ic50): 
                logging.warn(
                    "Invalid IC50 value %0.4f for %s w/ allele %s",
                    ic50,
                    epitope,
                    allele)
                continue 
            elif invalid_binding_score(rank) or rank > 100:
                logging.warn(
                    "Invalid percentile rank %s for %s w/ allele %s", 
                    rank, epitope, allele)
                continue 

            # keep track of original genetic variant that 
            # gave rise to this epitope
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
