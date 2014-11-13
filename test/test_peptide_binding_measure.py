from nose.tools import assert_raises, istest

from immuno.peptide_binding_measure import (
    ic50_binding_measure,
    percentile_binding_measure,
)

ic50_cutoff = 500
percentile_cutoff = 2.0

def test_ic50_value_cutoff():
    assert ic50_binding_measure.value_is_binder(300, cutoff=ic50_cutoff)
    assert not ic50_binding_measure.value_is_binder(600, cutoff=ic50_cutoff)


def test_percentile_value_cutoff():
    assert percentile_binding_measure.value_is_binder(1.0,
        cutoff=percentile_cutoff)
    assert not percentile_binding_measure.value_is_binder(3.0,
        cutoff=percentile_cutoff)

def test_raises_on_negative_ic50():
    with assert_raises(Exception):
        ic50_binding_measure.value_is_binder(-1, cutoff=ic50_cutoff)


def test_raises_on_negative_percentile():
    with assert_raises(Exception):
        percentile_binding_measure.value_is_binder(-1, cutoff=percentile_cutoff)

def test_raises_on_percentile_greater_than_100():
    with assert_raises(Exception):
        percentile_binding_measure.value_is_binder(101, cutoff=percentile_cutoff)

binder_prediction_record = {
    'MHC_IC50' : 39.0,
    'MHC_Percentile_Rank': 1.2,
}

nonbinder_prediction_record = {
    'MHC_IC50' : 2300.0,
    'MHC_Percentile_Rank' : 4.9,
}

def test_ic50_record_cutoff():
    assert ic50_binding_measure.record_is_binder(
        binder_prediction_record, cutoff=ic50_cutoff)

    assert not ic50_binding_measure.record_is_binder(
        nonbinder_prediction_record, cutoff=ic50_cutoff)


def test_percentile_record_cutoff():
    assert percentile_binding_measure.record_is_binder(
        binder_prediction_record,
        cutoff=percentile_cutoff)

    assert not percentile_binding_measure.record_is_binder(
        nonbinder_prediction_record,
        cutoff=percentile_cutoff)
