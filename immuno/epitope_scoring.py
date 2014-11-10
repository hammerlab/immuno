import numpy as np

from peptide_binding_measure import ic50_binding_measure


class EpitopeScorer(object):

    def __init__(
            self,
            binding_measure=ic50_binding_measure,
            cutoff=500.0,
            transformation_function=None,
            allele_weights=None):
        """
        Parameters
        ----------

        binding_measure : PeptideBindingMeasure
            Object with method record_is_binder(prediction_record, cutoff)

        cutoff : float
            Value above which we ignore epitopes (could also be above,
            depending on `prediction_type`)

        transformation_function : fn : float -> float (optional)
            Given the predicted binding quantity, transform it into a value
            between 0.0 and 1.0

        allele_weights : dict (optional)
            Associates each MHC allele name with an importance weight
        """
        self.binding_measure = binding_measure
        self.cutoff = cutoff
        self.transformation_function = transformation_function
        self.allele_weights = allele_weights

    def allele_weight(self, allele_name):
        if self.allele_weights:
            assert allele_name in allele_weights, \
                "Allele not in weight dictionary: %s" % (allele_name,)
            return allele_weights[allele_name]
        else:
            return 1.0

    def binding_value_score(self, value, allele_name=None):
        """
        Transform a binding prediction for a particular allele into a
        rescaled score (non-binders get a score of 0.0, strongest
        predicted binders should have a score of 1.0)

        Parameters
        ----------

        value : float

        allele_name : string (optional)
        """
        # if value is too weak (e.g. IC50 below 500nM), return a score of 0.0
        if not self.binding_measure.value_is_binder(value, self.cutoff):
            return 0.0

        allele_weight = self.allele_weight(allele_name) if allele_name else 1.0

        if self.transformation_function:
            score = self.transformation_function(value)
        else:
            score = 1.0

        return allele_weight * score

    def binding_record_score(self, record, allele_name=None):
        value = self.binding_measure.extract_value(record)
        return self.binding_value_score(value, allele_name=allele_name)

    def mean_binding_record_scores(self, epitope_binding_prediction_records):
        """

        Parameters
        ----------

        epitope_binding_prediction_records : dictionary of binding predictions
            Keys are alleles, values are binding prediction records with
            fields such as 'MHC_IC50' and 'MHC_Percentile_Rank'

        Returns sum of normalized binding prediction scores.
        """
        total = 0.0
        if isinstance(epitope_binding_prediction_records, dict):
            items = epitope_binding_prediction_records.iteritems()
        else:
            assert isinstance(epitope_binding_prediction_records, list)
            items = [
                (item['Allele'], item)
                for item in epitope_binding_prediction_records
            ]

        for allele, record in items:
            total += self.binding_record_score(record, allele_name=allele)
        return total / len(epitope_binding_prediction_records)


    def epitope_score(self, epitope):
        binding_records = epitope['MHC_Allele_Scores']
        return self.mean_binding_record_scores(binding_records)


class DecreasingLogisticFunction(object):
    def __init__(self, midpoint, width):
        assert width > 0
        self.midpoint = midpoint
        self.width = width

    def __call__(self, value):
        normalized = (float(value) - self.midpoint) / self.width
        # simplification of 1.0 - logistic(x) = logistic(-x)
        return 1.0 / (1.0 + np.exp(normalized))

# add up all the epitopes with IC50 <= 500nM
simple_ic50_epitope_scorer = EpitopeScorer(
    binding_measure=ic50_binding_measure,
    cutoff=500.0)

logistic_fn = DecreasingLogisticFunction(midpoint=300.0, width=75.0)

logistic_ic50_epitope_scorer = EpitopeScorer(
    binding_measure=ic50_binding_measure,
    cutoff=2000.0,
    transformation_function=logistic_fn,
)
