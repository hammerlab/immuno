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

import logging 
from copy import deepcopy

from epitopes.mutate import gene_mutation_description
import pandas as pd

from common import normalize_chromosome_name, is_valid_peptide
from ensembl import annotation, gene_names
from ensembl.transcript_variant import peptide_from_transcript_variant
from vcf import load_vcf
from snpeff import load_snpeff
from maf import load_maf
from fasta import load_fasta

def load_file(
        input_filename, 
        min_peptide_length = 9, 
        max_peptide_length = 31, 
        transcript_log_filename = None):
    """
    Load mutatated peptides from FASTA, VCF, or MAF file. 

    Parameters
    --------

    input_filename : str 

    min_peptide_length : int, optional 

    max_peptide_length : int, optional 

    transcript_log_filename : str, optional
        Write out transcript dataframe to CSV with the given filename 

    Returns a dataframe with columns:
        - chr : chomosome
        - pos : position in the chromosome
        - ref : reference DNA
        - alt : alternate DNA
        - info : gene name and entrez gene ID
        - stable_id_transcript : Ensembl transcript ID
        - SourceSequence : region of protein around mutation
        - MutationStart : first amino acid modified
        - MutationEnd : last mutated amino acid
        - MutationInfo : annotation i.e. V600E
    """

    if input_filename.endswith(".fasta") \
            or input_filename.endswith(".fa"):
        return load_fasta(input_filename, peptide_length = max_peptide_length)

    # VCF and MAF files give us the raw mutations in genomic coordinates
    if input_filename.endswith(".vcf"):
        vcf_df = load_vcf(
            input_filename, 
            min_peptide_length = min_peptide_length, 
            max_peptide_length = max_peptide_length)
    elif input_filename.endswith(".maf"):
        maf_df = load_maf(input_filename, max_peptide_length = max_peptide_length)
        # rename columns from MAF file to make them compatible with VCF names 
        vcf_df = pd.DataFrame({
            'pos' : maf_df['Start_Position'],
            'chr' : maf_df['Chromosome'],
            'id' : maf_df['dbSNP_RS'],
            'ref' : maf_df['Reference_Allele'].str.replace("-", ""),
            'alt' : maf_df['Tumor_Seq_Allele1'].str.replace("-", ""),
            'info' : maf_df['Hugo_Symbol'],
        })
    else:
        assert False, "Unrecognized file type %s" % input_filename

    assert len(vcf_df)  > 0, "No mutation entries for %s" % input_filename 
    
    vcf_df['chr'] = vcf_df.chr.map(normalize_chromosome_name)
    

    # annotate genomic mutations into all the possible known transcripts they might be on
    transcripts_df = annotation.annotate_vcf_transcripts(vcf_df)

    assert len(transcripts_df) > 0, "No annotated mutation entries for %s" % input_filename
    logging.info("Annotated input file %s has %d possible transcripts", input_filename, len(transcripts_df))
    
    new_rows = []

    group_cols = ['chr','pos', 'ref', 'alt', 'stable_id_transcript']

    seen_source_sequences = set([])
    for (chromosome, pos, ref, alt, transcript_id), group in \
            transcripts_df.groupby(group_cols):
        logging.info("---")
        logging.info("GENETIC MUTATION chr%s %s on transcript %s",
            chromosome, gene_mutation_description(pos, ref, alt), transcript_id)
        row = group.irow(0)
        padding = max_peptide_length - 1 
        if transcript_id:
            # logging.info("Getting peptide from transcript ID %s", transcript_id)
            seq, start, stop, annot = \
                peptide_from_transcript_variant(
                    transcript_id, pos, ref, alt,
                    padding = padding)
            assert isinstance(start, int), (start, type(start))
            assert isinstance(stop, int), (stop, type(stop))
        else:
            logging.info("Skipping transcript_id = %s, ref = %s, alt = %s, pos = %s" % (transcript_id, ref, alt, pos))
        if seq:
            if seq in seen_source_sequences:
                logging.info("Skipping %s because already seen sequence %s" % (gene_mutation_description(pos,ref,alt), seq))
                continue
            else:
                seen_source_sequences.add(seq)

            if '*' in seq:
                logging.warning(
                    "Found stop codon in peptide %s from transcript_id %s",
                    seq,
                    transcript_id)
            if not is_valid_peptide(seq):
                logging.warning(
                    "Invalid peptide sequence for transcript_id %s: %s",
                    transcript_id, 
                    seq)
            elif len(seq) < min_peptide_length:
                logging.info(
                    "Truncated peptide too short for transcript %s gene position %s '%s' > '%s' ", 
                    transcript_id, pos, ref, alt)
            else:
                row = deepcopy(row)
                row['SourceSequence'] = seq
                row['MutationStart'] = start
                row['MutationEnd'] = stop
                row['MutationInfo'] = annot
                logging.info("PEPTIDE MUTATION: %s", annot)
                try:
                    gene = gene_names.transcript_id_to_gene_name(transcript_id)
                except:
                    gene = gene_names.transcript_id_to_gene_id(transcript_id)

                row['Gene'] = gene
                new_rows.append(row)
    assert len(new_rows) > 0, "No mutations!"
    peptides = pd.DataFrame.from_records(new_rows)
    peptides['GeneInfo'] = peptides['info']
    peptides['TranscriptId'] = peptides['stable_id_transcript'] 

    transcripts_df = transcripts_df.merge(peptides)
    logging.info("Generated %d peptides from %s",
        len(transcripts_df), input_filename)

    # drop verbose or uninteresting columns from VCF
    for dumb_field in ('description_gene', 'filter', 'qual', 'id', 'name', 'info', 'stable_id_transcript'):
        if dumb_field in transcripts_df.columns:
            transcripts_df = transcripts_df.drop(dumb_field, axis = 1)
    if transcript_log_filename:
        transcripts_df.to_csv(transcript_log_filename, index=False)
    logging.info("---")
    logging.info("FILE LOADING SUMMARY")
    logging.info("---")
    logging.info("# original mutations: %d", len(vcf_df))
    logging.info("# mutations with annotations: %d", len(transcripts_df.groupby(['chr', 'pos', 'ref', 'alt'])))
    logging.info("# transcripts: %d", len(transcripts_df))
    return transcripts_df
