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

import cPickle
import logging
from os import environ 
from os.path import exists, split, join 

from mhc_common import compact_hla_allele_name
DEFAULT_PEPTIDE_DIR = environ.get("IMMUNO_THYMIC_PEPTIDES", join(split(__file__)[0], "thymic_peptides"))

class ImmunogenicityPredictor(object):

	def __init__(self, alleles, data_path = None):
		self.alleles = [compact_hla_allele_name(allele) for allele in alleles]
		if data_path is None:
			self.data_path = DEFAULT_PEPTIDE_DIR
		else:
			self.data_path = data_path 
		assert exists(self.data_path)


	def predict(self, peptide):
		"""
		Determine whether 9-mer peptide is immunogenic by checking for 
		membership of its 3rd through 8th residues in the "self"/thymic MHC 
		ligand sets for each HLA allele

