import pandas as pd
from immuno.group_epitopes import group_epitopes_dataframe
from immuno.peptide_binding_measure import (
        IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
)

"""
We're constructing a DataFrame with the following fields:
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
    - MHC_PercentileRank
"""
epitopes_df = pd.DataFrame()
epitopes_df['chr'] = ['1', '1', 'X', 'X',]
epitopes_df['pos'] = [10, 10, 2000, 2000]
epitopes_df['ref'] = ['A', 'A', '', '']
epitopes_df['alt'] = ['T', 'T', 'C', 'C']
epitopes_df['TranscriptId'] = [
    'ENST00000528762', 'ENST00000528762',
    'ENST00000544455', 'ENST00000544455'
]
epitopes_df['SourceSequence'] = [
    'ASIINFKELA', 'ASIINFKELA',
    'ASILLLVFYW', 'ASILLLVFYW',
]
epitopes_df['MutationStart'] = [3, 3, 5, 5]
epitopes_df['MutationEnd'] = [4, 4, 6, 6]
epitopes_df['GeneMutationInfo'] = ['A>T', 'A>T', 'insC', 'insC']
epitopes_df['PeptideMutationInfo'] = ['L>I', 'L>I', 'fs', 'fs']
epitopes_df['Gene'] = ['SMAD4', 'SMAD4', 'TP53', 'TP53']
epitopes_df['GeneInfo'] = [None, None, None, None]
epitopes_df['Epitope'] = ['SIINFKEL', 'SIINFKEL', 'SILLLVFY', 'SILLLVFY']
epitopes_df['EpitopeStart'] = [1, 1, 1, 1]
epitopes_df['EpitopeEnd'] = [10, 10, 10, 10]
epitopes_df['Allele']= [
    'HLA-A*02:01',
    'HLA-B*08:02',
    'HLA-A*02:01',
    'HLA-B*08:02'
]
epitopes_df[IC50_FIELD_NAME] = [0.9, 205.9, 5039.0, 112.9]
epitopes_df[PERCENTILE_RANK_FIELD_NAME] = [0.1, 9.2, 25.2, 3.4]

def test_group_epitopes_dataframe():
    grouped = group_epitopes_dataframe(epitopes_df)
    assert len(grouped) == 2
    assert isinstance(grouped, list)

    for elt in grouped:
        assert isinstance(elt, dict), \
            "Wrong type for %s : %s" % (elt, type(elt))
        assert elt['Gene'] in ("SMAD4", "TP53")
        assert "Epitopes" in elt, elt.keys()
        assert len(elt['Epitopes']) == 1

        epitopes = elt['Epitopes']
        assert isinstance(epitopes, list)

        epitope = epitopes[0]

        assert isinstance(epitope, dict), "Wrong type for epitope %s : %s" % (
            epitope, type(epitope))
        assert 'MHC_Allele_Scores' in epitope, epitope.keys()
        assert len(epitope['MHC_Allele_Scores']) == 2

if __name__ == '__main__':
    test_group_epitopes_dataframe()