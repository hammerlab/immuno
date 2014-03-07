#
# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
 

import vcf
import Bio
from Bio import SeqIO
from snpeff_effect import SnpEffEffect, Ensembl
import argparse


class Variant2Epitope():
  def __init__(self, vcf_file=None, cdna_file='Homo_sapiens.GRCh37.74.cdna.all.fa', cds_path='Homo_sapiens.GRCh37.74.cds.all.fa'):
    if cdna_file:
      self._ensembl = Ensembl(cdna_path = cdna_file, cds_path=cds_path)
      #self._transcripts = SeqIO.to_dict(SeqIO.parse(ensembl_file, 'fasta'))
      #print self._transcripts.items()[:20]
    if vcf_file:
      self._epitopes = self.generate_epitopes(vcf_file)
    
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
    mut_sequence[pos] = variant
    return mut_sequence.toseq()

  def generate_epitopes(self, vcf_file, window = 7):
    reader = vcf.Reader(open(vcf_file))
    epitopes = []
    for record in reader:
      effects = [SnpEffEffect(effect) for effect in record.INFO['EFF']]
      for effect in effects:
        if self._ensembl:
        #  print effect.effect
        #  print effect.full_line
          peptide = self.get_peptide(effect.transcriptId)
          if peptide:
            epitopes.append(peptide)
    
    return epitopes

  def get_peptide(self, transcriptId):
    (transcript_start, transcript_end, transcript) = self._ensembl.get_cdna(transcriptId)
    if transcript:
      # pos = record.POS - transcript_start 
      # start = pos - window
      # # end = pos + window + 1
      # # mutated_sequence = self._insert_mutation(transcript.seq, 0, str(record.ALT[0]))
      polymer = transcript.seq.translate().tomutable()

      peptide = polymer
      return peptide
    return None
    


def main(args):
  generator = Variant2Epitope(args.vcf)
  #rint generator.get_transcripts()
  #data.to_csv(args.output, index=False)

if __name__ == '__main__':
  # parser = argparse.ArgumentParser()
  # parser.add_argument("--vcf", help="vcf file to process")
  # parser.add_argument("--output", help="output file for data")

  # args = parser.parse_args()

  # main(args)
  generator = Variant2Epitope()
  print generator.get_peptide('ENST00000407673')



