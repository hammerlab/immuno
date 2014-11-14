
from ensembl.gene_names import transcript_id_to_transcript_name
from peptide_binding_measure import (
        IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
)

def group_epitopes_dataframe(scored_epitopes, use_transcript_name = True):
    """
    Given a DataFrame with fields:
        - chr
        - pos
        - ref
        - alt
        - TranscriptId
        - SourceSequence
        - MutationStart
        - MutationEnd
        - GeneMutationInfo
        - PeptideMutationInfo
        - Gene
        - GeneInfo
        - Epitope
        - EpitopeStart
        - EpitopeEnd
        - Allele
        - MHC_IC50
        - MHC_Percentile_Rank

    Group epitopes under their originating transcript and
    make nested lists of dictionaries to contain the binding scores
    for MHC alleles.

    Return a list of dictionaries for each mutated transcript with fields:
        - chr
        - pos
        - ref
        - alt
        - TranscriptId
        - SourceSequence
        - MutationStart
        - MutationEnd
        - GeneMutationInfo
        - PeptideMutationInfo
        - Gene
        - GeneInfo
        - Epitopes : list of dictionaries

    Each entry of the 'Epitopes' list contains the following fields:
        - 'Epitope'
        - 'EpitopeStart'
        - 'EpitopeEnd'
        - 'MHC_AlleleScores' : list of allele-specific entries

    Each entry of 'MHC_AlleleScores' the following fields:
        - 'Allele'
        - 'MHC_PercentileRank'
        - 'MHC_IC50'
    """

    peptides = []

    for (transcript_id, seq), transcript_group in \
            scored_epitopes.groupby(["TranscriptId", "SourceSequence"]):
        peptide_entry = {}
        peptide_entry['SourceSequence'] = seq
        peptide_entry["Peptide"] = seq

        if use_transcript_name:
            transcript_name =  transcript_id_to_transcript_name(transcript_id)
            peptide_entry['TranscriptId'] = transcript_name
        else:
            peptide_entry['TranscriptId'] = transcript_id

        head = transcript_group.to_records()[0]
        peptide_entry['chr'] = head.chr
        peptide_entry['pos'] = head.pos
        peptide_entry['ref'] = head.ref
        peptide_entry['alt'] = head.alt
        peptide_entry["MutationStart"] = head.MutationStart
        peptide_entry["MutationEnd"] = head.MutationEnd
        peptide_entry["GeneMutationInfo"] = head.GeneMutationInfo
        peptide_entry["PeptideMutationInfo"] = head.PeptideMutationInfo
        peptide_entry["GeneInfo"] = head.GeneInfo
        peptide_entry['Gene'] = head.Gene
        peptide_entry['Epitopes'] = []
        for (epitope, epitope_start, epitope_end), epitope_group in \
                transcript_group.groupby(
                    ['Epitope', 'EpitopeStart', 'EpitopeEnd']):
            epitope_entry = {
                'Epitope' : epitope,
                'EpitopeStart' : epitope_start,
                'EpitopeEnd' : epitope_end,
                'MHC_Allele_Scores' : []
            }
            seen_alleles = set([])
            for epitope_allele_row in epitope_group.to_records():
                allele = epitope_allele_row['Allele']
                if allele in seen_alleles:
                    logging.warn("Repeated entry %s", epitope_allele_row)
                    continue
                seen_alleles.add(allele)
                percentile_rank = epitope_allele_row[PERCENTILE_RANK_FIELD_NAME]
                ic50 = epitope_allele_row[IC50_FIELD_NAME]
                allele_entry = {
                    'Allele': allele,
                    PERCENTILE_RANK_FIELD_NAME : percentile_rank,
                    IC50_FIELD_NAME : ic50,
                }
                epitope_entry['MHC_Allele_Scores'].append(allele_entry)
            peptide_entry['Epitopes'].append(epitope_entry)
        peptides.append(peptide_entry)
    return peptides
