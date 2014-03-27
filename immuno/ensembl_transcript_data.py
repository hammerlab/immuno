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

import ensembl_download

class EnsemblReferenceData(object):
    """
    Singleton class which allows for lazy loading of reference
    cDNA and amino acid sequences of transcripts
    """

    def __init__(self):
        self._cdna_dict = None
        self._protein_dict = None

    def _load_cdna(self):
        self._cdna_dict = SeqIO.index(
            ensembl_download.download_cdna_transcripts(), 'fasta')

    def _load_peptide(self):
        self._protein_dict = SeqIO.index(
            ensembl_download.download_protein_transcripts(), 'fasta')

    def get_cdna(self, transcript_id):
        if self._cdna_dict is None:
            self._load_cdna()
        transcript = self._cdna_dict.get(transcript_id, None)
        return transcript

    def get_protein(self, transcript_id):
        if self._protein_dict is None:
            self._load_peptide()
        transcript = self._protein_dict.get(transcript_id, None)
        return transcript