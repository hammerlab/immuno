import numpy as np

from immuno.vaccine_peptides import select_vaccine_peptides
from immuno.peptide_binding_measure import (
    IC50_FIELD_NAME,
    PERCENTILE_RANK_FIELD_NAME,
)
from immuno.epitope_scoring import (
    simple_ic50_epitope_scorer,
    logistic_ic50_epitope_scorer
)

# the specific values in all these fields, aside from offsets into the
# source sequence, are nonsense. Don't expect them to make sense together.

candidates = [
    {
        # this field isn't normally here, just using it for testing
        'ExpectedRanking' : 0,
        'SourceSequence' : 'AAASIINFKELAAA',
        'MutationStart' : 5,
        'MutationEnd' : 6,
        'Gene' : 'BETTER',
        'GeneMutationInfo' : 'chr1 g.3944 C>T',
        'PeptideMutationInfo' : 'S399L',
        'TranscriptId' : 'BETTER',
        'Epitopes' : [
            {
                'Epitope': 'SIINKFEL',
                'EpitopeStart' : 3,
                'EpitopeEnd' : 12,
                'MHC_Allele_Scores' : [
                    {
                        'Allele' : 'HLA-A*02:01',
                        IC50_FIELD_NAME : 20.0,
                        PERCENTILE_RANK_FIELD_NAME : 0.3,
                    },
                    {
                        'Allele' : 'HLA-B*08:02',
                        IC50_FIELD_NAME : 2300.0,
                        PERCENTILE_RANK_FIELD_NAME : 50.0,
                    }
                ]

            }
        ]
    },
    {
        'ExpectedRanking' : 1,
        'SourceSequence' : 'AAAAATAAAAAAAAAAAAAAA',
        'MutationStart' : 5,
        'MutationEnd' : 6,
        'Gene' : 'WORSE',
        'GeneMutationInfo' : 'chrX g.2929494 A>T',
        'PeptideMutationInfo' : 'A201T',
        'TranscriptId' : 'WORSE-1',
        'Epitopes' : [
            {
                'Epitope': 'AAAATAAAA',
                'EpitopeStart' : 1,
                'EpitopeEnd' : 10,
                'MHC_Allele_Scores' : [
                    {
                        'Allele' : 'HLA-A*02:01',
                        IC50_FIELD_NAME : 600.0,
                        PERCENTILE_RANK_FIELD_NAME : 2.5,
                    },
                    {
                        'Allele' : 'HLA-B*08:02',
                        IC50_FIELD_NAME : 3300.0,
                        PERCENTILE_RANK_FIELD_NAME : 50.0,
                    }
                ]
            },
            {
                'Epitope': 'AAAAAAAAA',
                'EpitopeStart' : 6,
                'EpitopeEnd' : 16,
                'MHC_Allele_Scores' : [
                    {
                        'Allele' : 'HLA-A*02:01',
                        IC50_FIELD_NAME : 9600.0,
                        PERCENTILE_RANK_FIELD_NAME : 12.0,
                    },
                    {
                        'Allele' : 'HLA-B*08:02',
                        IC50_FIELD_NAME : 14300.0,
                        PERCENTILE_RANK_FIELD_NAME : 50.0,
                    }
                ]
            }
        ]
    }
]

VACCINE_PEPTIDE_LENGTH = 10

def check_vaccine_peptides(results):
    assert len(results) == len(candidates)
    last_score = np.inf
    for i, result in enumerate(results):
        assert len(result['VaccinePeptide']) == VACCINE_PEPTIDE_LENGTH, result
        score = result['MutantEpitopeScore']
        assert score <= last_score, \
            "Vaccine peptides not sorted in decreasing order"
        last_score = score
        assert i == result['ExpectedRanking'], \
            "Unexpected ranking %d for %s" % (i, result)

def test_select_vaccine_peptides_simple_scorer():
    results = select_vaccine_peptides(
        candidates,
        vaccine_peptide_length=VACCINE_PEPTIDE_LENGTH,
        epitope_scorer=simple_ic50_epitope_scorer)
    check_vaccine_peptides(results)


def test_select_vaccine_peptides_logistic_scorer():
    results = select_vaccine_peptides(
        candidates,
        vaccine_peptide_length=VACCINE_PEPTIDE_LENGTH,
        epitope_scorer=logistic_ic50_epitope_scorer)
    check_vaccine_peptides(results)

if __name__ == '__main__':
    test_select_vaccine_peptides_simple_scorer()
    test_select_vaccine_peptides_logistic_scorer()
