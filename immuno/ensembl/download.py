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

from os.path import join, exists
from os import environ

import appdirs
import pandas as pd
from epitopes.download import fetch_data, ensure_dir, DATA_DIR

DATA_DIR = environ.get("IMMUNO_DATA_DIR", appdirs.user_cache_dir("immuno"))

STANDARD_CONTIGS = set([
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14',
    '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M'
])

# TODO: Generously describe what all these files are

GENE_HEADER = [
    'gene_id', 'biotype', 'analysis_id',
    'seq_region_id', 'seq_region_start', 'seq_region_end',
    'seq_region_strand', 'display_xref_id',
    'source', 'status', 'description', 'is_current',
    'canonical_transcript_id', 'stable_id',
    'version', 'created_date', 'modified_date'
]

GENE_DATA_URL = \
"ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/gene.txt.gz"

SEQ_REGION_HEADER = ['seq_region_id', 'name', 'coord_system_id']

SEQ_REGION_DATA_URL = \
"ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/seq_region.txt.gz"


EXON_HEADER = [
    "exon_id", "seq_region_id", "seq_region_start",
    "seq_region_end", "seq_region_strand", "phase", "end_phase",
    "is_current", "is_constitutive", "stable_id", "version",
    "created_date", "modified_date"
]


EXON_DATA_URL = \
"ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/exon.txt.gz"


TRANSCRIPT_HEADER = [
    "transcript_id", "gene_id", "analysis_id",
    "seq_region_id", "seq_region_start", "seq_region_end",
    "seq_region_strand", "display_xref_id", "biotype", "status",
     "description", "is_current", "canonical_translation_id", "stable_id",
     "version", "created_date", "modified_date"
]

TRANSCRIPT_DATA_URL = \
"ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/transcript.txt.gz"

EXON_TRANSCRIPT_DATA_URL = \
"ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/exon_transcript.txt.gz"

def download_transcript_metadata(
        output_file = 'gene_exon_transcript.tsv',
        filter_contigs = STANDARD_CONTIGS):
    ensure_dir(DATA_DIR)
    full_path = join(DATA_DIR, output_file)
    if not exists(full_path):
        GENE_DATA_PATH = fetch_data('gene.txt', GENE_DATA_URL)
        SEQ_REGION_DATA_PATH = fetch_data('seq_region.txt', SEQ_REGION_DATA_URL)
        EXON_DATA_PATH = fetch_data('exon.txt', EXON_DATA_URL)
        TRANSCRIPT_DATA_PATH = fetch_data('transcript.txt', TRANSCRIPT_DATA_URL)
        EXON_TRANSCRIPT_DATA_PATH = \
            fetch_data('exon_transcript.txt', EXON_TRANSCRIPT_DATA_URL)

        seqregion = pd.read_csv(
            SEQ_REGION_DATA_PATH,
            sep='\t',
            names = SEQ_REGION_HEADER,
            index_col=False)

        def in_filter_contigs(x):
            return x in filter_configs

        if filter_contigs:
            seqregion = seqregion[seqregion['name'].map(in_filter_contigs)]

        # TODO: Ask Arun what's going on here, leave a good explanation
        gene = pd.read_csv(
            GENE_DATA_PATH,
            sep='\t',
            names = GENE_HEADER,
            index_col=False)
        seqregion_gene = gene.merge(seqregion, on='seq_region_id')
        transcript = pd.read_csv(
            TRANSCRIPT_DATA_PATH, sep='\t',
            names = TRANSCRIPT_HEADER,
            index_col=False)
        gene_transcript = transcript.merge(
            seqregion_gene,
            on='gene_id',
            suffixes = ['', '_gene'])


        exon = pd.read_csv(
            EXON_DATA_PATH,
            sep='\t',
            names=EXON_HEADER,
            index_col=False)
        exon_transcript = pd.read_csv(
            EXON_TRANSCRIPT_DATA_PATH,
            sep='\t',
            names=['exon_id', 'transcript_id', 'rank'])
        exon_w_exon_transcript = pd.merge(
            exon,
            exon_transcript,
            on='exon_id',
            suffixes=('_exon', '_et'))

        exon_w_transcript = pd.merge(
            exon_w_exon_transcript,
            gene_transcript,
            on='transcript_id',
            suffixes=('_exon', '_transcript'))

        exon_cols = [
            'name',
            'stable_id_gene',
            'description_gene',
            'seq_region_start_gene',
            'seq_region_end_gene',
            'stable_id_transcript',
            'seq_region_start_transcript',
            'seq_region_end_transcript',
            'stable_id_exon',
            'seq_region_start_exon',
            'seq_region_end_exon']
        exon_data = exon_w_transcript[exon_cols]
        exon_data.to_csv(full_path, index=False, sep='\t')
    return full_path

PROTEIN_TRANSCIPT_URL = \
'ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.74.pep.all.fa.gz'

PROTEIN_TRANSCRIPT_FILE = 'Homo_sapiens.GRCh37.74.pep.all.fa.gz'

def download_protein_transcripts():
    """
    Downloads a FASTA file containing amino acid sequences of human
    transcripts, return local path.
    """
    return fetch_data(PROTEIN_TRANSCRIPT_FILE, PROTEIN_TRANSCIPT_URL)

CDNA_TRANSCRIPT_URL = \
'ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.74.cdna.all.fa.gz'
CDNA_TRANSCRIPT_FILE = 'Homo_sapiens.GRCh37.74.cdna.all.fa.gz'

def download_cdna_transcripts():
    """
    Download a FASTA file containing cDNA sequences for each known transcript,
    return its local path.
    """
    return fetch_data(CDNA_TRANSCRIPT_FILE, CDNA_TRANSCRIPT_URL)

CDS_TRANSCRIPTS_URL = \
'ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.74.cds.all.fa.gz'

CDS_TRANSCRIPTS_FILE = 'Homo_sapiens.GRCh37.74.cds.all.fa.gz'

def download_cds_transcripts():
    """
    TODO: What are CDS transcripts?
    """
    return fetch_data(CDS_TRANSCRIPTS_FILE, CDS_TRANSCRIPTS_URL)