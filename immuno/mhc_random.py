import pandas as pd
import random
def generate_scored_epitopes(mutated_regions):
	records = []
    # if wer'e not running the MHC prediction then we have to manually
    # extract 9mer substrings
    for _, row in mutated_regions.iterrows():
        seq = row.SourceSequence
        epitope_length = 9
        for i in xrange(len(seq) - epitope_length + 1):
            record = {}
            record['Epitope'] = seq[i:i+epitope_length]
            record['EpitopeStart'] = i
            record['EpitopeEnd'] = i + epitope_length
            record['SourceSequence'] = seq
            record['MutationStart'] = row['MutationStart']
            record['MutationEnd'] = row['MutationEnd']
            record['MutationInfo'] = row['MutationInfo']
            record['GeneInfo'] = row['GeneInfo']
            record['TranscriptId'] = row['TranscriptId']
            record['Gene'] = row['Gene']
            record['MHC_PercentileRank'] = random.randint(0,99)
            record['MHC_IC50'] = random.rand() * 10000.0
            records.append(record)
    scored_epitopes = pd.DataFrame.from_records(records)
    return scored_epitopes
