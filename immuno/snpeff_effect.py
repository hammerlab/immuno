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

import re

class SnpEffEffect():
  """

  EFF=EXON
  (MODIFIER|||||FTCD|processed_transcript|CODING|ENST00000498355|6|1),
  FRAME_SHIFT(HIGH||-|-221|495|FTCD|protein_coding|CODING|ENST00000355384|6|1)

  """
  EFFECT_PATTERN = re.compile(r'(\w+)\(.+\)')
  EFFECT_FIELDS = re.compile(r'\w+\((.+)\)')


  def __init__(self, line = None):
    # self.effect 0
    # self.impact 1
    # self.functionalClass
    # self.codonChange
    # self.aminoAcidChange
    # self.geneName
    # self.geneBiotype
    # self.coding
    # self.transcriptID
    # self.exonID

    if line is not None:
      self.parse(line)


  def parse(self, record):
    self.effect = SnpEffEffect.EFFECT_PATTERN.match(record).group(1)
    fields = SnpEffEffect.EFFECT_FIELDS.match(record).group(1)
    self.full_line = fields
    fields = fields.split("|")

    self.impact = fields[0]
    self.functional_class = fields[1]
    self.transcript_id = fields[8]

if __name__ == '__main__':
  SnpEffEffect("EXON(MODIFIER|||||FTCD|processed_transcript|CODING|ENST00000498355|6|1)")

