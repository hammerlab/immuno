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

import logging
from os import environ, listdir
from os.path import exists, split, join

from mhc_common import compact_hla_allele_name
from peptide_binding_measure import IC50_FIELD_NAME


DEFAULT_PEPTIDE_DIR = environ.get(
    "IMMUNO_THYMIC_PEPTIDES",
    join(split(__file__)[0], "thymic_peptides"))

THYMIC_DELETION_FIELD_NAME = 'ThymicDeletion'

def _load_allele_mapping_dict(path):
    """
    Since some alleles have identical peptide sets as others, we compress
    the stored data by only retaining one allele from each equivalence class
    and using a mappings file to figure out which allele is retained.
    """
    result = {}
    with open(path, 'r') as f:
        for line in f.read().split("\n"):
            if len(line) > 0:
                k, v = line.split("\t")
            result[k] = v
    return result


class ImmunogenicityPredictor(object):

    """
    Predict whether some T-cell in a person's circulating repertoire could
    recognize a particular pattern. The subset of the 'self' proteome which
    binds to an individual's HLA alleles tells us which T-cells were removed
    by negative selection. T-cells inspect peptides more strongly along
    interior residues (positions 3-8), so we restrict our query only to those
    positions.
    """

    def __init__(
            self,
            alleles,
            data_path = DEFAULT_PEPTIDE_DIR,
            binding_threshold = 500,
            first_position = 3,
            last_position = 8):
        """
        Parameters
        --------

        alleles : list of strings

        data_path : str, optional

        first_position : int, optional
            Start position for extracting substring of
            query peptide (indexed starting from 1)

        last_position : int, optional
            Last position for extracting substring of
            query peptide (indexed starting from 1)
        """

        self.binding_threshold = binding_threshold
        self.first_position = first_position
        self.last_position = last_position
        self.alleles = {
            compact_hla_allele_name(allele) for allele in alleles
        }
        self.data_path = data_path

        assert exists(self.data_path), \
            "Directory with thymic peptides (%s) does not exist" % \
            self.data_path

        available_alleles = listdir(self.data_path)

        mappings_file_path = join(self.data_path, 'mappings')
        if exists(mappings_file_path):
            self.allele_mappings = \
                _load_allele_mapping_dict(mappings_file_path)
        else:
            self.allele_mappings = \
                dict(zip(available_alleles, available_alleles))

        self.peptide_sets = {}

        for allele in self.alleles:
            if allele not in self.allele_mappings:
                logging.warn(
                    "No MHC peptide set available for HLA allele %s", allele)
                continue
            else:
                logging.info(
                    "Loading thymic MHC peptide set for HLA allele %s", allele)

            filename = self.allele_mappings[allele]
            assert filename in available_alleles, \
                "No MHC peptide set available for HLA allele %s (file = %s)" % \
                    (allele,filename)

            with open(join(self.data_path, filename), 'r') as f:
                peptide_set = {l for l in f.read().split("\n") if len(l) > 0}
            self.peptide_sets[allele] = peptide_set

    def predict(self, peptides_df):
        """
        Determine whether 9-mer peptide is immunogenic by checking

        1) that the epitope binds strongly to a particular MHC allele
        2) the "core" of the peptide (positions 3-8) don't overlap with any
           other peptides in the self/thymic MHC ligand sets of that HLA allele

        Returns DataFrame with two extra columns:
            - ThymicDeletion: Was this epitope deleted during thymic selection
              (and thus can't be recognize by T-cells)?
            - Immunogenic: Is this epitope a sufficiently strong binder that
              wasn't deleted during thymic selection?
        """

        thymic_peptide_sets = self.peptide_sets.values()

        # assume a peptide is non-immunogenic unless not in thymic sets
        # We do this in case some alleles are missing, resulting in all
        # their associated ligands being considered non-immunogenic
        peptides_df[THYMIC_DELETION_FIELD_NAME] = True
        for i in xrange(len(peptides_df)):
            row = peptides_df.ix[i]
            peptide = row.Epitope
            allele = compact_hla_allele_name(row.Allele)
            if allele in self.peptide_sets:
                # positions in the epitope are indexed starting from 1 to
                # match immunology nomenclature
                substring = \
                    peptide[self.first_position - 1 : self.last_position]
                peptides_df[THYMIC_DELETION_FIELD_NAME].ix[i] = \
                    substring in self.peptide_sets[allele]

        peptides_df["Immunogenic"] = \
            ~peptides_df[THYMIC_DELETION_FIELD_NAME] & \
            (peptides_df[IC50_FIELD_NAME] <= self.binding_threshold)

        return peptides_df

