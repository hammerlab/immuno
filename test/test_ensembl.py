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

import immuno.ensembl_annotation as ensembl
from immuno.ensembl_transcript_data import EnsemblReferenceData
import pandas as pd

# Test case
# Transcript Id: ENST00000453024
# Gene Id : ENSG00000168036
# 17 Exons
# Length 2,841
# Chrom : 3
# Start : 41,240,936
# End :  41,281,227

ref_data = EnsemblReferenceData()

def test_load_transcripts():
    genomic_transcript = "ENST00000453024"
    transcript = ref_data.get_cdna(genomic_transcript)
    print(str(transcript.seq))
    assert(transcript is not None)
    assert(len(transcript.seq) == 2841)


def test_load_protein_transcript():
    protein_transcript = "ENSP00000427553"
    transcript = ref_data.get_protein(protein_transcript)
    print(str(transcript.seq))
    assert(transcript is not None)
    assert(transcript.seq[0] == 'M')
    assert(transcript.seq[-1] == 'E')

def test_load_exon_from_transcript():
    transcript_id = "ENST00000453024"
    exons = ensembl.get_exons_from_transcript(transcript_id)
    assert(exons.shape[0] == 17)

    transcript_id = "ENST00000405570"
    exons = ensembl.get_exons_from_transcript(transcript_id)
    assert(exons.shape[0] == 16)

def test_load_exon_from_transcript_lengt():

    transcript_id = 'ENST00000405570'

    transcript = ref_data.get_cdna(transcript_id)

    exons = ensembl.get_exons_from_transcript(transcript_id)
    exons['length'] = exons['seq_region_end_exon'] - exons['seq_region_start_exon'] + 1

    assert(exons['length'].sum() == len(transcript.seq))

def test_get_gene_from_pos():
    variant = {
        'chr' : '3',
        'pos' : 41250936,
        'ref' : 'A',
        'alt' : 'C'
    }
    vcf = pd.DataFrame.from_records([variant])
    genes = ensembl.annotate_genes(vcf)
    assert( "ENSG00000168036" in set(genes['stable_id_gene']))

def test_get_transcript_from_pos():
    variant = {
        'chr' : '3',
        'pos' : 41250936,
        'ref' : 'A',
        'alt' : 'C'
    }
    vcf = pd.DataFrame.from_records([variant])
    transcripts_ids = ensembl.annotate_transcripts(vcf)
    assert( "ENST00000453024" in set(transcripts_ids['stable_id_transcript']))

def test_get_all_transcript_from_pos():
    variant = {
        'chr' : '3',
        'pos' : 41275636,
        'ref' : 'G',
        'alt' : 'A'
    }
    vcf = pd.DataFrame.from_records([variant])
    transcripts_ids = ensembl.annotate_transcripts(vcf)
    transcript_ids = set(transcripts_ids['stable_id_transcript'])
    assert( "ENST00000405570" in transcript_ids)
    assert( "ENST00000396183" in transcript_ids)
    assert( "ENST00000349496" in transcript_ids)
    assert( "ENST00000453024" in transcript_ids)
    assert( "ENST00000396185" in transcript_ids)

def test_get_transcript_index_from_pos():
    variant = {
        'chr' : '3',
        'pos' : 41275636,
        'ref' : 'G',
        'alt' : 'A'
    }
    transcript_id = 'ENST00000405570'

    idx = ensembl.get_transcript_index_from_pos(41275636, transcript_id)
    assert(idx == 1686)

    transcript = ref_data.get_cdna(transcript_id)
    assert(transcript.seq[idx] == variant['ref'])

def test_interval_search():
    intervals = [ (7,13), (17,19), (21, 24), (35, 45), (47, 50), (60, 70)]
    idx = ensembl.get_idx_from_interval(7, intervals)
    print idx
    assert(idx == 0)

    idx = ensembl.get_idx_from_interval(13, intervals)
    print idx
    assert(idx == 6)

    idx = ensembl.get_idx_from_interval(14, intervals)
    print idx
    assert(idx is None)

    idx = ensembl.get_idx_from_interval(12, intervals)
    print idx
    assert(idx == 5)

    idx = ensembl.get_idx_from_interval(17, intervals)
    print idx
    assert(idx == 7)

    idx = ensembl.get_idx_from_interval(18, intervals)
    print idx
    assert(idx == 8)

    idx = ensembl.get_idx_from_interval(23, intervals)
    print idx
    assert(idx == 12)

    idx = ensembl.get_idx_from_interval(51, intervals)
    print idx
    assert(idx is None)


if __name__ == '__main__':
  from dsltools import testing_helpers
  testing_helpers.run_local_tests()
