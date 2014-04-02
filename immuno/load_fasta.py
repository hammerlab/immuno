import pandas as pd
from Bio import SeqIO

def peptides_from_fasta(fasta_files, peptide_length):
    epitope_dataframes = []
    peptides = []
    source_seqs = []
    filenames = []
    for fasta_file in fasta_files:
        fasta_data = SeqIO.parse(fasta_file, 'fasta')
        seqs = list(set([e.seq for e in fasta_data]))
        for seq in seqs:
            curr_peptides = peptide_substrings(seq, peptide_length)
            peptides.extend(curr_peptides)
            source_seqs.extend([seq] * len(curr_peptides))
            filenames.extend([fasta_file] * len(curr_peptides))
    assert len(peptides) == len(source_seqs) == len(filenames)
    return pd.DataFrame({
        'Peptide': peptides,
        'SourceSequence': source_seqs,
        'Filename': filenames,
    })
