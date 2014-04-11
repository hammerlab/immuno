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

from Bio.Seq import Seq
import pandas as pd
from epitopes.mutate import mutate_protein_from_transcript

import immuno.ensembl.annotation as ensembl
from immuno.ensembl.transcript_data import EnsemblReferenceData

# Test case CTNNB1
# Transcript Id: ENST00000453024
# Gene Id : ENSG00000168036
# 17 Exons
# Length 2,841
# Chrom : 3
# Start : 41,240,936
# End :  41,281,227

ref_data = EnsemblReferenceData()

def test_complement_base():
    assert ensembl.complement("G") == "C"

def test_complement_seq():
    assert ensembl.complement("TCTCATCCAGGTACCAGCCAATG") == "AGAGTAGGTCCATGGTCGGTTAC"

def test_get_strand_CASP9():
    genomic_transcript = "ENST00000333868"
    forward = ensembl.is_forward_strand(genomic_transcript)
    assert(forward == False)

def test_get_strand_CTNNB1():
    genomic_transcript = "ENST00000453024"
    forward = ensembl.is_forward_strand(genomic_transcript)
    assert(forward == True)

def test_load_CTNNB1_cdna_transcript():
    genomic_transcript = "ENST00000453024"
    transcript = ref_data.get_cdna(genomic_transcript)
    assert(transcript is not None)
    assert(len(transcript) == 2841), (transcript, len(transcript))


def test_load_CTNNB1_protein_transcript():
    protein_transcript = "ENSP00000427553"
    transcript = ref_data.get_protein(protein_transcript)
    assert(transcript is not None)
    assert(transcript[0] == 'M'), (transcript, len(transcript))
    assert(transcript[-1] == 'E'), (transcript, len(transcript))

def test_load_CTNNB1_exon_from_transcript():
    transcript_id = "ENST00000453024"
    exons = ensembl.get_exons_from_transcript(transcript_id)
    assert(exons.shape[0] == 17)

    transcript_id = "ENST00000405570"
    exons = ensembl.get_exons_from_transcript(transcript_id)
    assert(exons.shape[0] == 16)

def test_load_CTNNB1_exon_from_transcript_length():

    transcript_id = 'ENST00000405570'
    transcript = ref_data.get_cdna(transcript_id)
    exons = ensembl.get_exons_from_transcript(transcript_id)
    exons['length'] = exons['seq_region_end_exon'] - exons['seq_region_start_exon'] + 1
    assert(exons['length'].sum() == len(transcript)), exons

def test_load_SMAD4_cdna_transcript():
    transcript_id = "ENST00000342988"
    transcript = ref_data.get_cdna(transcript_id)
    assert transcript is not None
    assert len(transcript) == 8769, len(transcript)
    assert transcript[0] == 'A', transcript[0]
    assert transcript[-1] == 'T', transcript[-1]

def test_get_gene_from_pos():
    variant = {
        'chr' : '3',
        'pos' : 41250936,
        'ref' : 'A',
        'alt' : 'C'
    }
    vcf = pd.DataFrame.from_records([variant])
    genes = ensembl.annotate_vcf_genes(vcf)
    assert( "ENSG00000168036" in set(genes['stable_id_gene']))

def test_get_transcript_from_pos():
    variant = {
        'chr' : '3',
        'pos' : 41250936,
        'ref' : 'A',
        'alt' : 'C'
    }
    vcf = pd.DataFrame.from_records([variant])
    transcripts_ids = ensembl.annotate_vcf_transcripts(vcf)
    assert( "ENST00000453024" in set(transcripts_ids['stable_id_transcript']))

def test_get_all_transcript_from_pos():
    variant = {
        'chr' : '3',
        'pos' : 41275636,
        'ref' : 'G',
        'alt' : 'A'
    }
    vcf = pd.DataFrame.from_records([variant])
    transcripts_ids = ensembl.annotate_vcf_transcripts(vcf)
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
    idx = ensembl.get_transcript_index_from_pos(
        41275636, transcript_id, skip_untranslated_region = False)
    assert(idx == 1686), idx

    transcript = ref_data.get_cdna(transcript_id)
    assert(transcript[idx] == variant['ref'])

def test_get_5prime_utr_length_RET():

    transcript_id = "ENST00000355710"
    exons = ensembl.get_exons_from_transcript(transcript_id)

    utr_length = ensembl.get_five_prime_utr_length(exons)
    print utr_length
    assert(utr_length == 232)

def test_get_3prime_utr_length_RET():

    transcript_id = "ENST00000355710"
    exons = ensembl.get_exons_from_transcript(transcript_id)
 
    utr_length = ensembl.get_three_prime_utr_length(exons)
    print utr_length
    assert(utr_length == 2082)

def test_get_5prime_utr_length_CTNNB1():

    transcript_id = "ENST00000405570"
    exons = ensembl.get_exons_from_transcript(transcript_id)
 
    utr_length = ensembl.get_five_prime_utr_length(exons)
    print utr_length
    assert(utr_length == 156)

def test_get_3prime_utr_length_CTNNB1():

    transcript_id = "ENST00000405570"
    exons = ensembl.get_exons_from_transcript(transcript_id)

    utr_length = ensembl.get_three_prime_utr_length(exons)
    print utr_length
    assert(utr_length == 11)

def test_get_5prime_utr_length_reverse_strand_CASP9():

    transcript_id = "ENST00000333868"
    exons = ensembl.get_exons_from_transcript(transcript_id)
 
    utr_length = ensembl.get_five_prime_utr_length(exons, forward = False)
    print utr_length
    assert(utr_length == 95)

def test_get_3prime_utr_length_reverse_strand_CASP9():

    transcript_id = "ENST00000333868"
    exons = ensembl.get_exons_from_transcript(transcript_id)

    utr_length = ensembl.get_three_prime_utr_length(exons, forward = False)
    print utr_length
    assert(utr_length == 673)

def test_get_transcript_and_mutate_vcf():
    variant = {
        'chr' : '10',
        'pos' : 43617416,
        'ref' : 'T',
        'alt' : 'C'
    }

    vcf = pd.DataFrame.from_records([variant])
    transcripts_ids = ensembl.annotate_vcf_transcripts(vcf)

    transcript_ids = set(transcripts_ids['stable_id_transcript'])
    assert( "ENST00000355710" in transcript_ids)
    assert( "ENST00000340058" in transcript_ids)

    transcript_id = "ENST00000355710"


    cdna_idx = ensembl.get_transcript_index_from_pos(
        variant['pos'], transcript_id, skip_untranslated_region = False)
    assert cdna_idx is not None
    assert cdna_idx < 5569
    cdna_transcript = ref_data.get_cdna(transcript_id)
    assert(cdna_transcript[cdna_idx] == variant['ref'])

    cds_idx = ensembl.get_transcript_index_from_pos(
        variant['pos'], transcript_id, skip_untranslated_region = True)
    assert cds_idx is not None
    cds_transcript = ref_data.get_cds(transcript_id)
    assert(cds_transcript[cds_idx] == variant['ref'])

    region = mutate_protein_from_transcript(
            cds_transcript,
            cds_idx,
            variant['ref'],
            variant['alt'],
            padding = 10)
    assert region is not None
    assert len(region.seq) == 21, (region.seq, len(region.seq))
    assert region.seq == 'RSQGRIPVKWTAIESLFDHIY'

def test_interval_search():
    intervals = [ (7,13), (17,19), (21, 24), (35, 45), (47, 50), (60, 70)]
    idx = ensembl.get_idx_from_interval(7, intervals)

    assert(idx == 0), idx

    idx = ensembl.get_idx_from_interval(13, intervals)
    assert(idx == 6), idx

    idx = ensembl.get_idx_from_interval(14, intervals)
    assert(idx is None), idx

    idx = ensembl.get_idx_from_interval(12, intervals)
    assert(idx == 5), idx

    idx = ensembl.get_idx_from_interval(17, intervals)
    assert(idx == 7), idx

    idx = ensembl.get_idx_from_interval(18, intervals)
    assert(idx == 8), idx

    idx = ensembl.get_idx_from_interval(23, intervals)
    assert(idx == 12), idx

    idx = ensembl.get_idx_from_interval(51, intervals)
    assert(idx is None), idx

def test_peptide_from_transcript():
    """
    test_peptide_from_transcript:

    """
    transcript_id = 'ENST00000333868'
    cds_transcript = ref_data.get_cds(transcript_id)

if __name__ == '__main__':
  from dsltools import testing_helpers
  testing_helpers.run_local_tests()
