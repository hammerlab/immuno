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
from epitopes.download import fetch_data, ensure_dir

DATA_DIR = environ.get("IMMUNO_DATA_DIR", appdirs.user_cache_dir("immuno"))

STANDARD_CONTIGS = set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', \
                       '19', '20', '21', '22', 'X', 'Y', 'M'])


GENE_HEADER = [ 'gene_id', 'biotype', 'analysis_id',
        'seq_region_id', 'seq_region_start', 'seq_region_end',
        'seq_region_strand', 'display_xref_id',
        'source', 'status', 'description', 'is_current',
        'canonical_transcript_id', 'stable_id',
        'version', 'created_date', 'modified_date']

SEQ_REGION_HEADER = ['seq_region_id', 'name', 'coord_system_id']

EXON_HEADER = ["exon_id", "seq_region_id", "seq_region_start",
                "seq_region_end", "seq_region_strand", "phase", "end_phase",
                "is_current", "is_constitutive", "stable_id", "version",
                "created_date", "modified_date"]

TRANSCRIPT_HEADER = ["transcript_id", "gene_id", "analysis_id",
                                "seq_region_id", "seq_region_start", "seq_region_end",
                                "seq_region_strand", "display_xref_id", "biotype", "status",
                                "description", "is_current", "canonical_translation_id", "stable_id",
                                "version", "created_date", "modified_date"]

def download_transcript_metadata(output_file = 'gene_exon_transcript.tsv', filter_contigs = STANDARD_CONTIGS):
    ensure_dir(DATA_DIR)
    full_path = join(DATA_DIR, output_file)
    if not exists(full_path):
        GENE_DATA_PATH = fetch_data('gene.txt', "ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/gene.txt.gz")
        SEQ_REGION_DATA_PATH = fetch_data('seq_region.txt', "ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/seq_region.txt.gz")
        EXON_DATA_PATH = fetch_data('exon.txt', "ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/exon.txt.gz")
        TRANSCRIPT_DATA_PATH = fetch_data('transcript.txt', "ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/transcript.txt.gz")
        EXON_TRANSCRIPT_DATA_PATH = fetch_data('exon_transcript.txt', "ftp://ftp.ensembl.org/pub/release-74/mysql/homo_sapiens_core_74_37/exon_transcript.txt.gz")

        seqregion = pd.read_csv(SEQ_REGION_DATA_PATH, sep='\t', names = SEQ_REGION_HEADER, index_col=False)
        if filter_contigs:
            seqregion = seqregion[ seqregion['name'].map(lambda x: x in filter_contigs)]

        gene = pd.read_csv(GENE_DATA_PATH, sep='\t', names = GENE_HEADER, index_col=False)
        seqregion_gene = gene.merge(seqregion, on='seq_region_id')
        transcript = pd.read_csv(TRANSCRIPT_DATA_PATH, sep='\t',names = TRANSCRIPT_HEADER, index_col=False)
        gene_transcript = transcript.merge(seqregion_gene, on='gene_id', suffixes = ['', '_gene'])


        exon = pd.read_csv(EXON_DATA_PATH, sep='\t', names=EXON_HEADER, index_col=False)
        exon_transcript = pd.read_csv(EXON_TRANSCRIPT_DATA_PATH, sep='\t', names=['exon_id', 'transcript_id', 'rank'])
        exon_w_exon_transcript = pd.merge(exon, exon_transcript, on='exon_id', suffixes=('_exon', '_et'))

        exon_w_transcript = pd.merge(exon_w_exon_transcript, gene_transcript, on='transcript_id', suffixes=('_exon', '_transcript'))
        exon_data = exon_w_transcript[['name', 'stable_id_gene', 'description_gene', 'seq_region_start_gene', 'seq_region_end_gene', 'stable_id_transcript',
                        'seq_region_start_transcript', 'seq_region_end_transcript', 'stable_id_exon', 'seq_region_start_exon', 'seq_region_end_exon']]

        exon_data.to_csv(full_path, index=False, sep='\t')

    return full_path

def download_protein_transcripts():
    path = fetch_data('Homo_sapiens.GRCh37.74.pep.all.fa.gz', 'ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.74.pep.all.fa.gz')
    return path

def download_cdna_transcripts():
    path = fetch_data('Homo_sapiens.GRCh37.74.cdna.all.fa.gz', 'ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.74.cdna.all.fa.gz')
    return path

def download_cds_transcripts():
    path = fetch_data('Homo_sapiens.GRCh37.74.cds.all.fa.gz', 'ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.74.cds.all.fa.gz')
    return path

