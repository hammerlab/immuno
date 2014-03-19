# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from copy import deepcopy

import pandas as pd
import numpy as np

import epitopes.mutate as mutate

from ensembl_transcript_data import EnsemblReferenceData
import ensembl_annotation
from snpeff_effect import SnpEffEffect


_ensembl = _ensembl = EnsemblReferenceData()

def liftover(chr):
    if chr.startswith('chr'):
        chr = chr[-1]
        if chr =='M':
            return 'MT'
        return chr
    return chr

def vcf2dataframe(vcffile):
    with open(vcffile) as fd:
        lines_to_skip = 0
        while next(fd).startswith('#'):
            lines_to_skip += 1
    header = ['chr', 'pos', 'id', 'ref', 'alt','qual', 'filter', 'info']
    df = pd.read_csv(vcffile, sep='\t', skiprows=lines_to_skip, names=header, usecols=header, dtype={'pos' : np.int32})
    df['chr'] = df.chr.map(liftover)
    return df

def generate_epitopes_from_vcf(input_file):
    vcf_df = vcf2dataframe(input_file)
    transcripts_df = ensembl_annotation.annotate_transcripts(vcf_df)
    def epitope_from_annotation(row):
        transcript_id = row['stable_id_transcript']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        if transcript_id:
            return generate_epitope_from_transcript(transcript_id, pos, ref, alt)
        else:
            return None

    epitopes = transcripts_df.apply(epitope_from_annotation, axis=1)
    transcripts_df['Epitope'] = pd.Series(epitopes, index=transcripts_df)
    return transcripts_df


def generate_epitopes_from_snpeff(snpeff_annotated_file, window=7):
    data = vcf2dataframe(snpeff_annotated_file)
    results = []
    for idx, variant in data.iterrows():
      variant = dict(variant)
      info = variant['info']
      del variant['info']
      for effect in parse_effects(info):
        variant = deepcopy(variant)
        variant['Transcript'] = effect.transcript_id   
        epitope = generate_epitope_from_transcript(effect.transcript_id, variant['pos'], variant['ref'], variant['alt'])
        if epitope:
          variant['Epitope'] = epitope
          results.append(variant)

    return pd.DataFrame.from_records(results)

def generate_epitopes_from_protein_transcript(transcript_id, pos, ref, variant):
    transcript = _ensembl.get_protein(transcript_id)
    if transcript:
        try:
            return mutate.mutate_protein_from_transcript(transcript.seq, pos, ref, variant)
        except:
            return None
    return None

def generate_epitope_from_transcript(transcript_id, pos, ref, variant):
    transcript = _ensembl.get_cdna(transcript_id)
    if transcript:
      idx = ensembl_annotation.get_transcript_index_from_pos(pos, transcript_id)
      if idx:
        try:
            return mutate.mutate_protein_from_transcript(transcript.seq, pos, ref, variant)
        except:
            return None
    return None

def parse_effects(info_field):
    info_fields = info_field.split(";")
    for field in info_fields:
        if field.startswith("EFF"):
            key, value = field.split("=")
            return [SnpEffEffect(effect) for effect in value.split(",")]

def parse_aa_sift(self, sift_string):
# ENSP00000427553:T109R
    protein_transcript_id, variation = sift_string.split(":")
    aa_orig = variation[0]
    aa_variant = variation[-1]
    aa_change_pos = int(variation[1:-1])
    return (protein_transcript_id, aa_orig, aa_variant, aa_change_pos)

# def generate_epitopes_from_annotations(self, annotated_file, window=7):
#     data = pd.read_csv(annotated_file)
#     results = []
#     for idx, variant in data.iterrows():
#         variant = dict(variant)
#         aapos_sift = str(variant['aapos_sift'])
#         if aapos_sift != '.' and aapos_sift != 'nan':
#             (protein_transcript_id, aa_orig, aa_variant, aa_change_pos) = self.parse_aa_sift(aapos_sift)
#             variant['Epitope']= generate_epitopes_from_protein_transcript(protein_transcript_id, aa_change_pos, aa_variant)
#             results.append(variant)
#         else:
#             transcripts = variant['ensembl_transcriptid'].split(";")
#             for transcript_id in transcripts:
#                 variant = deepcopy(variant)
#                 epitope = generate_epitope_from_transcript(transcript_id, variant['pos'], variant['ref'], variant['alt'], window=window)
#                 if epitope:
#                     variant['Epitope'] = epitope
#                     results.append(variant)
#     return pd.DataFrame.from_records(results)
