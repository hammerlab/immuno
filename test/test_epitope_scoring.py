from immuno.epitope_scoring import (
    DecreasingLogisticFunction,
    simple_ic50_epitope_scorer,
    logistic_ic50_epitope_scorer
)
from immuno.peptide_binding_measure import (
    IC50_FIELD_NAME,
    PERCENTILE_RANK_FIELD_NAME
)

from nose.tools import raises

def test_logistic():
    fn = DecreasingLogisticFunction(midpoint=300.0, width=75.0)
    assert fn(0) > 0.95, fn(0)
    assert fn(50000) < 0.001, fn(50000)

def test_simple_ic50_epitope_scorer_allele_weight():
    # since default scorer doesn't use allele weights, expect 1.0
    assert simple_ic50_epitope_scorer.allele_weight("HLA-A*02:01") == 1.0


def test_simple_ic50_epitope_scorer_binding_value_score():
    # since default scorer doesn't use allele weights, expect accessing
    # the allele weight to raise an exception
    cutoff = simple_ic50_epitope_scorer.cutoff
    assert simple_ic50_epitope_scorer.binding_value_score(cutoff-1) == 1
    assert simple_ic50_epitope_scorer.binding_value_score(cutoff+1) == 0

def test_logistic_ic50_scorer_binding_value_score():
    cutoff = logistic_ic50_epitope_scorer.cutoff
    # just being stronger than the cutoff should still make for a very low
    # score
    assert simple_ic50_epitope_scorer.binding_value_score(cutoff-1) < 0.001

    # whereas very high affinity scores like 1nM should be scored highly
    assert simple_ic50_epitope_scorer.binding_value_score(1.0) > 0.9


binder_prediction_record = {
    IC50_FIELD_NAME : 39.0,
    PERCENTILE_RANK_FIELD_NAME : 1.2,
}

nonbinder_prediction_record = {
    IC50_FIELD_NAME : 5300.0,
    PERCENTILE_RANK_FIELD_NAME : 8.9,
}


def test_simple_ic50_scorer_binding_record_score():
    assert simple_ic50_epitope_scorer.binding_record_score(
        nonbinder_prediction_record) == 0

    assert simple_ic50_epitope_scorer.binding_record_score(
        binder_prediction_record) == 1

def test_logistic_ic50_scorer_binding_record_score():
    assert logistic_ic50_epitope_scorer.binding_record_score(
        nonbinder_prediction_record) < 0.001

    assert logistic_ic50_epitope_scorer.binding_record_score(
        binder_prediction_record) > 0.9

epitope = {
    'Epitope': 'SIINFKEL',
    'MHC_Allele_Scores' : [
        {
            'Allele' : 'HLA-A*02:01',
            IC50_FIELD_NAME : 32.0,
            PERCENTILE_RANK_FIELD_NAME : 9.8
        },
        {
            'Allele' : 'HLA-B*08:02',
            IC50_FIELD_NAME : 1900.0,
            PERCENTILE_RANK_FIELD_NAME : 39.2,
        }
    ]
}

def test_simple_ic50_scorer_epitope():
    score = simple_ic50_epitope_scorer.epitope_score(epitope)
    assert score == 0.5

def test_logistic_ic50_scorer_epitope():
    score = logistic_ic50_epitope_scorer.epitope_score(epitope)
    assert score > 0.25
    assert score < 0.5
