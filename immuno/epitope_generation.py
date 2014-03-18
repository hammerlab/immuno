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
# limitations under the License.     import argparse
from copy import deepcopy

import pandas as pd
import numpy as np

from ensembl import Ensembl
from snpeff_effect import SnpEffEffect



class Variant2Epitope(object):
    def __init__(self,
                 protein_path = 'Homo_sapiens.GRCh37.74.pep.all.fa',
                 cdna_file = 'Homo_sapiens.GRCh37.74.cdna.all.fa',
                 cds_path = 'Homo_sapiens.GRCh37.74.cds.all.fa'):
        if cdna_file:
          self._ensembl = Ensembl(cdna_path = cdna_file, cds_path=cds_path, protein_path=protein_path)

    @staticmethod
    def vcf2dataframe(vcffile):
        with open(vcffile) as fd:
            lines_to_skip = 0
            while next(fd).startswith('#'):
                lines_to_skip += 1
        header = ['chr', 'pos', 'id', 'ref', 'alt','qual', 'filter', 'info']
        df = pd.read_csv(vcffile, sep='\t', skiprows=lines_to_skip, names=header, usecols=header, dtype={'pos' : np.int32})
        return df

    def _insert_mutation(self, sequence, pos, ref, variant):
        mut_sequence = sequence.tomutable()
        try:
          mut_sequence[pos:pos+len(ref)] = variant
        except IndexError:
          print pos, sequence
        return mut_sequence


    def generate_epitopes_from_annotations(self, annotated_file, window=7):
        data = pd.read_csv(annotated_file)
        results = []
        for idx, variant in data.iterrows():
            variant = dict(variant)
            aapos_sift = str(variant['aapos_sift'])
            if aapos_sift != '.' and aapos_sift != 'nan':
                (protein_transcript_id, aa_orig, aa_variant, aa_change_pos) = self.parse_aa_sift(aapos_sift)
                variant['Epitope']= self.generate_epitopes_from_protein_transcript(protein_transcript_id, aa_change_pos, aa_variant)
                results.append(variant)
            else:
                transcripts = variant['ensembl_transcriptid'].split(";")
                for transcript_id in transcripts:
                    variant = deepcopy(variant)
                    epitope = self.generate_epitope_from_transcript(transcript_id, variant['pos'], variant['ref'], variant['alt'], window=window)
                    if epitope:
                        variant['Epitope'] = epitope
                        results.append(variant)
        return pd.DataFrame.from_records(results)

    def generate_epitopes_from_snpeff(self, snpeff_annotated_file, window=7):
        data = Variant2Epitope.vcf2dataframe(snpeff_annotated_file)
        results = []
        for idx, variant in data.iterrows():
          variant = dict(variant)
          info = variant['info']
          del variant['info']
          for effect in Variant2Epitope.parse_effects(info):
            variant = deepcopy(variant)
            variant['Transcript'] = effect.transcript_id
            epitope = self.generate_epitope_from_transcript(
                effect.transcript_id, variant['pos'],
                variant['ref'],
                variant['alt'],
                window=window)
            if epitope:
              variant['Epitope'] = epitope
              results.append(variant)

        return pd.DataFrame.from_records(results)

    @staticmethod
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


    def generate_epitopes_from_protein_transcript(
            self,
            transcript_id,
            pos,
            ref,
            variant,
            window=10):
        (transcript_start, transcript_end, transcript) = \
            self._ensembl.get_protein(transcript_id)
        if transcript:
            return self.get_flanking_epitope(
                self._insert_mutation(transcript.seq, pos, ref, variant),
                pos, window=window)
        return None

    def get_flanking_epitope(self, sequence, pos, window=10):
        flanked = sequence[pos-window:pos+window].toseq()
        return flanked

    def generate_epitope_from_transcript(
            self,
            transcript_id,
            pos,
            ref,
            variant,
            window=10):
        (transcript_start, transcript_end, transcript) = \
            self._ensembl.get_cdna(transcript_id)
        if transcript:
          idx = self._ensembl.get_transcript_index_from_pos(pos, transcript_id)
          if idx:
              mutated_sequence = self._insert_mutation(
                    transcript.seq, idx, ref, variant)
              polymer = mutated_sequence.toseq().translate().tomutable()
              epitope = self.get_flanking_epitope(polymer, idx/3, window=window)
              return str(epitope)
        return None
