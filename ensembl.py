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

from Bio import SeqIO
import pandas as pd

class Ensembl():

  def __init__(self, dna_path = None, cds_path = None, protein_path = None, cdna_path = None, index_peptides = False, index_lengths=[8,9,10,11]):
    self.load_ensembl(dna_path, cds_path, protein_path, cdna_path)
    if index_peptides:
      self._indexed = True
      self._lengths = index_lengths
      self._index_peptides()

  def load_ensembl( self, dna_path = None, cds_path = None, protein_path = None, cdna_path= None ):
    if cds_path:  
      self._cds_dict = SeqIO.index(cds_path, "fasta")
    if dna_path:
      self._gene_dict = SeqIO.index(dna_path, "fasta")
    if cdna_path:
      self._cdna_dict = SeqIO.index(cdna_path, "fasta")
    if protein_path:
      self._protein_dict = SeqIO.parse(protein_path, "fasta")
    exon = pd.read_csv('exon.txt', sep='\t',  names = 
                                           [ "exon_id",
                                            "seq_region_id",
                                            "seq_region_start",
                                            "seq_region_end",
                                            "seq_region_strand",
                                            "phase",
                                            "end_phase",
                                            "is_current",
                                            "is_constitutive",
                                            "stable_id",
                                            "version",
                                            "created_date",
                                            "modified_date" ]

    )
    transcript = pd.read_csv('transcript.txt',na_values=['\N'], sep='\t',names = 
                                           [  "transcript_id",
                                              "gene_id",
                                              "analysis_id",
                                              "seq_region_id",
                                              "seq_region_start",
                                              "seq_region_end",
                                              "seq_region_strand",
                                              "display_xref_id",
                                              "biotype",
                                              "status",
                                              "description",
                                              "is_current",
                                              "canonical_translation_id",
                                              "stable_id",
                                              "version",
                                              "created_date",
                                              "modified_date" ]

    )
    exon_transcript = pd.read_csv('exon_transcript.txt', sep='\t', names=['exon_id', 'transcript_id', 'rank'])

    exon_w_exon_transcript = pd.merge(exon, exon_transcript, on='exon_id', suffixes=('_exon', '_et'))
    self._exon_w_transcript = pd.merge(exon_w_exon_transcript, transcript, on='transcript_id', suffixes=('_exon', '_transcript'))



  def get_cds(self, transcriptId):
    transcript = self._cds_dict.get(transcriptId, None)
    if transcript:
      (chr, ref, length, start, end, something) = self._parse_description(transcript.description)
      return (start, end, transcript)
    return (-1, -1, transcript)

  def get_cdna(self, transcriptId):
    transcript = self._cdna_dict.get(transcriptId, None)
    if transcript:
      (chr, ref, length, start, end, something) = self._parse_description(transcript.description)
      return (start, end, transcript)
    return (-1, -1, transcript)

  def get_protein(self, transcript_id):
    transcript = self._protein_dict.get(transcript_id, None)
    if transcript:
      (chr, ref, length, start, end, something) = self._parse_description(transcript.description)
      return (start, end, transcript)
    return (-1, -1, transcript)

  def get_exons_from_transcript(self, transcript_id):
    exons = self._exon_w_transcript[self._exon_w_transcript['stable_id_transcript'] == transcript_id]
    return exons[['stable_id_exon', 'seq_region_start_exon', 'seq_region_end_exon']]

  def get_transcript_index_from_pos(self, pos, transcript_id):
    exons = self.get_exons_from_transcript(transcript_id).sort(columns=['seq_region_end_exon'])
    transcript_idx = 0
    for exon in exons:
      if pos > exon['seq_region_end_exon']:
        transcript_idx += exon['seq_region_end_exon'] - exon['seq_region_start_exon']
      elif pos < exon['seq_region_end_exon'] and pos > exon['seq_region_start_exon']:
        return transcript_idx + (pos - exon['seq_region_start_exon'])
      else:
        ## error some
        return -1

  def _parse_description(self, description):
    description_parts = description.split()
    position_information = description_parts[2]

    #chromosome:GRCh37:9:82319718:82333878:1

    (chr, ref, length, start, end, something) = position_information.split(":")
    return (chr, ref, length, int(start), int(end), something)

  def reference_contains_peptide(self, peptide):
    if not self._indexed:
      self._index_peptides()
    k = len(peptide)
    if peptide in self._reference_peptides[k]:
      return True
    return False

  def _generate_polymers(self, sequence, n=9):
    for i in xrange(len(sequence) - n):
      yield sequence[i:i+n]

  def _index_peptides(self):
    self._reference_peptides = {}
    for l in self._lengths:
      self._reference_peptides[l] = set()
      for peptide in self._protein_dict:
        self._reference_peptides[l] = self._reference_peptides[l].union(self._generate_polymers(peptide, l))



