import numpy as np


from immuno.epitope_scoring import (
    simple_ic50_epitope_scorer,
    logistic_ic50_epitope_scorer
)
from immuno.peptide_binding_measure import (
    IC50_FIELD_NAME,
    PERCENTILE_RANK_FIELD_NAME,
)
from immuno.vaccine_peptides import (
    clamp,
    is_mutant_epitope,
    generate_candidate_vaccine_peptides,
    select_vaccine_peptide,
    select_vaccine_peptides
)
from immuno.immunogenicity import THYMIC_DELETION_FIELD_NAME


def test_clamp():
    assert clamp(-1, 10) == 0
    assert clamp(5, 10) == 5
    assert clamp(10, 10) == 10
    assert clamp(11, 10) == 10

def test_is_mutant_epitope():
    mutation_position = 12
    mutant_epitope = {
        'EpitopeStart' : 10,
        'EpitopeEnd' : 19,
        'Epitope' : 'SIINKFKEL',
        THYMIC_DELETION_FIELD_NAME : False,

    }
    assert is_mutant_epitope(
        mutant_epitope,
        mutation_start=mutation_position,
        mutation_end=mutation_position+1)

    mutant_edge_epitope  = {
        'EpitopeStart' : 12,
        'EpitopeEnd' : 21,
        'Epitope' : 'LLLLLLLLL',
        THYMIC_DELETION_FIELD_NAME : False,

    }
    assert is_mutant_epitope(
        mutant_edge_epitope,
        mutation_start=mutation_position,
        mutation_end=mutation_position+1)

    wildtype_epitope = {
        'EpitopeStart' : 1,
        'EpitopeEnd' : 10,
        'Epitope' : 'AAAAAAAAA',
        THYMIC_DELETION_FIELD_NAME : False,
    }
    assert not is_mutant_epitope(
        wildtype_epitope,
        mutation_start=mutation_position,
        mutation_end=mutation_position+1
    )

    wildtype_edge_epitope =  {
        'EpitopeStart' : 3,
        'EpitopeEnd' : 12,
        'Epitope' : 'AAAAAAAAA',
        THYMIC_DELETION_FIELD_NAME : False,
    }

    assert not is_mutant_epitope(
        wildtype_edge_epitope,
        mutation_start=mutation_position,
        mutation_end=mutation_position+1
    )

    thymically_deleted_epitope = {
        'EpitopeStart' : 9,
        'EpitopeEnd' : 18,
        'Epitope' : 'TSIINKFE',
        THYMIC_DELETION_FIELD_NAME : True,
    }

    assert not is_mutant_epitope(
        thymically_deleted_epitope,
        mutation_start=mutation_position,
        mutation_end=mutation_position+1
    )




VACCINE_PEPTIDE_LENGTH = 10
SOURCE_SEQUENCE = 'AAAAAAAAASIINKFEL'

def make_epitopes(source_sequence, epitope_length=9):
    """
    Construct epitopes with better binding strength toward the end
    """

    epitopes = []
    for i in xrange(len(source_sequence)-epitope_length+1):
        # IC50s will be:
        # 2000, 1000, 666, 500, 400, 333, 285, 250.0, 222
        ic50 = 2000.0 / (i+1)

        # percentiles will be:
        # 8, 4, 2.66, 2 1.6, 1.33, 1.14, 1, 0.88
        percentile = (len(source_sequence) - epitope_length) / float(i+1)
        epitope = {
            'Epitope': source_sequence[i:i+epitope_length],
            'EpitopeStart': i,
            'EpitopeEnd': i+epitope_length,
            'MHC_Allele_Scores' : [
                {
                    'Allele' : 'HLA-A*02:01',
                    IC50_FIELD_NAME : ic50,
                    PERCENTILE_RANK_FIELD_NAME : percentile,
                }
            ],
            THYMIC_DELETION_FIELD_NAME : False,
        }
        epitopes.append(epitope)
    return epitopes

def test_generate_candidate_vaccine_peptides():
    epitopes = make_epitopes(SOURCE_SEQUENCE, epitope_length=9)
    mutation_start = len(SOURCE_SEQUENCE) - 4
    mutation_end = len(SOURCE_SEQUENCE) - 3
    candidates = generate_candidate_vaccine_peptides(
        seq=SOURCE_SEQUENCE,
        epitopes=epitopes,
        mutation_start=mutation_start,
        mutation_end=mutation_end,
        epitope_scorer=logistic_ic50_epitope_scorer,
        result_length=VACCINE_PEPTIDE_LENGTH,
        padding=0)

    expected_n_candidates = len(SOURCE_SEQUENCE) - VACCINE_PEPTIDE_LENGTH + 1
    assert len(candidates) == expected_n_candidates, \
        "Expected %d candidate vaccine peptides but only got %d: %s" % (
            expected_n_candidates,
            len(candidates),
            candidates
        )

    # Check each candidate vaccine peptide
    # if it contains any mutant epitopes it should have a non-zero
    # mutant epitope score.
    # Also check that vaccine peptides overlapping the mutation are
    # marked with a non-zero number of mutant residues.
    for i, candidate in enumerate(candidates):
        if i <= mutation_start - VACCINE_PEPTIDE_LENGTH:
            assert candidate['NumMutantResidues'] == 0, candidate
        else:
            assert candidate['NumMutantResidues'] > 0, candidate

        mutant_binder = None
        for epitope in epitopes:
            epitope_start = epitope['EpitopeStart']
            epitope_end = epitope['EpitopeEnd']
            overlaps_mutation = (
                epitope_start < mutation_end
                and epitope_end > mutation_start
            )
            contained_in_vaccine_peptide = (
                epitope_start >= candidate['VaccinePeptideStart'] and
                epitope_end <= candidate['VaccinePeptideEnd']
            )
            if contained_in_vaccine_peptide and overlaps_mutation:
                ic50 = epitope['MHC_Allele_Scores'][0][IC50_FIELD_NAME]
                epitope_score = \
                    logistic_ic50_epitope_scorer.binding_value_score(ic50)
                if epitope_score > 0:
                    mutant_binder = epitope
                    break

        if mutant_binder:
            assert candidate['MutantEpitopeScore'] > 0, \
                """
                Candidate peptide #%d (%s)
                    -- from %d to %d
                    -- was given mutant epitope score 0.0
                    -- contains epitope: %s
                """ % (
                    i+1,
                    candidate['VaccinePeptide'],
                    candidate['VaccinePeptideStart'],
                    candidate['VaccinePeptideEnd'],
                    mutant_binder
                )
        else:
            assert candidate['MutantEpitopeScore'] == 0, candidate

def test_select_vaccine_peptide():
    epitopes = make_epitopes(SOURCE_SEQUENCE)
    unpadded_result = select_vaccine_peptide(
        SOURCE_SEQUENCE,
        epitopes,
        mutation_start=len(SOURCE_SEQUENCE)-3,
        mutation_end=len(SOURCE_SEQUENCE)-2,
        epitope_scorer=logistic_ic50_epitope_scorer,
        result_length=VACCINE_PEPTIDE_LENGTH,
        padding=0)

    unpadded_expected = SOURCE_SEQUENCE[-VACCINE_PEPTIDE_LENGTH:]

    assert unpadded_result['VaccinePeptide'] == unpadded_expected, \
        "Expected vaccine peptide %s but got %s " % (
            unpadded_expected,
            unpadded_result)
    assert unpadded_result['NumMutantResidues'] == 1

    padded_result = select_vaccine_peptide(
        SOURCE_SEQUENCE,
        epitopes,
        mutation_start=len(SOURCE_SEQUENCE)-5,
        mutation_end=len(SOURCE_SEQUENCE)-4,
        epitope_scorer=logistic_ic50_epitope_scorer,
        result_length=VACCINE_PEPTIDE_LENGTH,
        padding=1)

    padded_expected = SOURCE_SEQUENCE[-VACCINE_PEPTIDE_LENGTH-1:-1]

    assert padded_result['VaccinePeptide'] == padded_expected, \
        "Expected vaccine peptide %s but got %s" % (
            padded_expected,
            padded_result)
    assert padded_result['NumMutantResidues'] == 1


# the specific values in all these fields, aside from offsets into the
# source sequence, are nonsense. Don't expect them to make sense together.

candidate_source_sequences = [
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

def check_vaccine_peptides(results):
    assert len(results) == len(candidate_source_sequences)
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
        candidate_source_sequences,
        vaccine_peptide_length=VACCINE_PEPTIDE_LENGTH,
        epitope_scorer=simple_ic50_epitope_scorer)
    check_vaccine_peptides(results)


def test_select_vaccine_peptides_logistic_scorer():
    results = select_vaccine_peptides(
        candidate_source_sequences,
        vaccine_peptide_length=VACCINE_PEPTIDE_LENGTH,
        epitope_scorer=logistic_ic50_epitope_scorer)
    check_vaccine_peptides(results)

if __name__ == '__main__':
    test_clamp()
    test_is_mutant_epitope()
    test_generate_candidate_vaccine_peptides()
    test_select_vaccine_peptide()
    test_select_vaccine_peptides_simple_scorer()
    test_select_vaccine_peptides_logistic_scorer()
