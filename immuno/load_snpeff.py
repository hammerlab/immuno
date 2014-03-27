from copy import deepcopy

from load_vcf import vcf_to_dataframe
from snpeff_effect import SnpEffEffect


def _parse_effects(info_field):
    info_fields = info_field.split(";")
    for field in info_fields:
        if field.startswith("EFF"):
            key, value = field.split("=")
            return [SnpEffEffect(effect) for effect in value.split(",")]
    return None

def peptides_from_snpeff(snpeff_annotated_file, window=7):
    data = vcf_to_dataframe(snpeff_annotated_file)
    results = []
    for idx, variant in data.iterrows():
      variant = dict(variant)
      info = variant.pop('info')
      for effect in _parse_effects(info):
        variant = deepcopy(variant)
        variant['Transcript'] = effect.transcript_id
        peptide = peptide_from_transcript(
            effect.transcript_id,
            variant['pos'],
            variant['ref'],
            variant['alt'])
        if peptide:
          variant['Peptide'] = peptide
          results.append(variant)
    return pd.DataFrame.from_records(results)

