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

from copy import deepcopy
import logging

from vcf import load_vcf
from snpeff_effect import SnpEffEffect

from ensembl.transcript_variant import peptide_from_transcript_variant

def _parse_effects(info_field):
    info_fields = info_field.split(";")
    for field in info_fields:
        if field.startswith("EFF"):
            key, value = field.split("=")
            return [SnpEffEffect(effect) for effect in value.split(",")]
    logging.warning("Couldn't find EFF info in %s", info_field)
    return None

def load_snpeff(snpeff_annotated_file, window=7):
    data = load_vcf(snpeff_annotated_file)
    results = []
    for idx, variant in data.iterrows():
      variant = dict(variant)
      info = variant.pop('info')
      for effect in _parse_effects(info):
        variant = deepcopy(variant)
        variant['Transcript'] = effect.transcript_id
        peptide = peptide_from_transcript_variant(
            effect.transcript_id,
            variant['pos'],
            variant['ref'],
            variant['alt'])
        if peptide:
          variant['Peptide'] = peptide
          results.append(variant)

    return pd.DataFrame.from_records(results)

