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

import argparse

import pandas as pd
import vcf

from ensembl import Ensembl
from snpeff import SnpEffEffect

class Variant2Epitope():
  def __init__(self, 
               protein_path = 'Homo_sapiens.GRCh37.74.pep.all.fa', 
               cdna_file = 'Homo_sapiens.GRCh37.74.cdna.all.fa', 
               cds_path = 'Homo_sapiens.GRCh37.74.cds.all.fa'):
    if cdna_file:
      self._ensembl = Ensembl(cdna_path = cdna_file, cds_path=cds_path, protein_path=protein_path)

  def parse_snpeff(self, vcf_file):
    reader = vcf.Reader(open(vcf_file))
    effects = []
    for record in reader:
      effects += [SnpEffEffect(effect) for effect in record.INFO['EFF']]
    return effects

  def get_transcripts(self):
    transcriptIds = [effect.transcriptId for effect in self._effects] 
    return transcriptIds

  def _insert_mutation(self, sequence, pos, variant):
    mut_sequence = sequence.tomutable()
    try:
      mut_sequence[pos] = variant
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
      # else:
      #   transcripts = variant['ensembl_transcriptid'].split(";")
      #   for transcript_id in transcripts:
      #     variant['Epitope']=  self.generate_epitope_from_transcript(transcript_id, variant['pos'], variant['alt'])
      #     results.append(variant)

    return pd.DataFrame.from_records(results)



  def parse_aa_sift(self, sift_string):
    # ENSP00000427553:T109R
    print sift_string
    protein_transcript_id, variation = sift_string.split(":")
    aa_orig = variation[0]
    aa_variant = variation[-1]
    aa_change_pos = int(variation[1:-1])
    return (protein_transcript_id, aa_orig, aa_variant, aa_change_pos)


  def generate_epitopes_from_protein_transcript(self, transcript_id, pos, variant, window=7):
    (transcript_start, transcript_end, transcript) = self._ensembl.get_protein(transcript_id)
    if transcript:
      return self.get_flanking_epitope(self._insert_mutation(transcript.seq, pos, variant), pos, window=window)
    return None

  # def generate_epitopes(self, vcf_file, window = 7):
  #   reader = vcf.Reader(open(vcf_file))
  #   epitopes = []
  #   for record in reader:
  #     effects = [SnpEffEffect(effect) for effect in record.INFO['EFF']]
  #     for effect in effects:
  #       if self._ensembl:
  #         peptide = self.get_peptide(effect.transcriptId)
  #         if peptide:
  #           epitopes.append(peptide)
    
  #   return epitopes

  def get_flanking_epitope(self, sequence, pos, window=7):
    return sequence[pos-window:pos+window].toseq()

  def generate_epitope_from_transcript(self, transcript_id, pos, variant, window=7):
    (transcript_start, transcript_end, transcript) = self._ensembl.get_cdna(transcript_id)
    if transcript:
      idx = self._ensembl.get_transcript_index_from_pos(transcript_id)
      mutated_sequence = self._insert_mutation(transcript.seq, idx, variant)
      polymer = mutated_sequence.translate().tomutable()
      return self.get_flanking_epitope(polymer, idx/3)
    return None
    


def main(args):
  generator = Variant2Epitope(args.vcf)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("--vcf", help="vcf file to process")
  parser.add_argument("--output", help="output file for data")

  args = parser.parse_args()

  main(args)
