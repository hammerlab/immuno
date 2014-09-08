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
from collections import OrderedDict 

from epitopes.mutate import gene_mutation_description
import pandas as pd

from common import normalize_chromosome_name, is_valid_peptide
from ensembl import annotation, gene_names
from ensembl.transcript_variant import peptide_from_transcript_variant
from vcf import load_vcf
from snpeff import load_snpeff
from maf import load_maf
from fasta import load_fasta

def maf_to_vcf(maf_df):
    """
    Convert DataFrame with columns from MAF file to DataFrame with columns 
    from VCF 
    """
    # rename columns from MAF file to make them compatible with VCF names 
    ref = maf_df['Reference_Allele'].str.replace("-", "")
    alt1 = maf_df['Tumor_Seq_Allele1'].str.replace("-", "")
    vcf_df = pd.DataFrame({
        'pos' : maf_df['Start_Position'],
        'chr' : maf_df['Chromosome'],
        'id' : maf_df['dbSNP_RS'],
        'ref' : maf_df['Reference_Allele'].str.replace("-", ""),
        'alt' : alt1, 
        'info' : maf_df['Hugo_Symbol'],
    })
    # use the second tumor allele if the first matches the reference
    alt2 = maf_df['Tumor_Seq_Allele2'].str.replace("-", "")
    use_alt2 = alt1 == ref
    vcf_df['alt'][use_alt2] = alt2[use_alt2]
    return vcf_df

def tab_to_vcf(tab_df):
    """
    Convert variant file with the following columns:
    - hgncSymbol
    - chrom
    - pos 
    - ref 
    - alt
    - dbsnpId
    """
    return pd.DataFrame({
        'pos' : tab_df['pos'],
        'ref' : tab_df['ref'],
        'alt' : tab_df['alt'],
        'chr' : tab_df['chrom'],
        'info' : tab_df['hgncSymbol'],
        'id' : tab_df['dbsnpId']
    })

def expand_transcripts(
        vcf_df, patient_id, min_peptide_length=9, max_peptide_length=31):
    """
    Applies genomic variants to all possible transcripts. 

    Parameters
    --------

    vcf_df : DataFrame 
        Required to have basic variant columns (chr, pos, ref, alt)

    patient_id : str 

    min_peptide_length : int 

    max_peptide_length : int 
    """

    assert len(vcf_df)  > 0, "No mutation entries for %s" % patient_id 
    
    vcf_df['chr'] = vcf_df.chr.map(normalize_chromosome_name)
    

    # annotate genomic mutations into all the possible 
    # known transcripts they might be on
    transcripts_df = annotation.annotate_vcf_transcripts(vcf_df)

    assert len(transcripts_df) > 0, \
        "No annotated mutation entries for %s" % patient_id
    logging.info(
        "Annotated input %s has %d possible transcripts",
         patient_id,
         len(transcripts_df))
    
    new_rows = []

    group_cols = ['chr','pos', 'ref', 'alt', 'stable_id_transcript']

    seen_source_sequences = set([])

    # for each genetic variant in the source file, 
    # we're going to print a string describing either the resulting 
    # protein variant or whatever error prevented us from getting a result 
    variant_report = OrderedDict()

    for (chromosome, pos, ref, alt, transcript_id), group in \
            transcripts_df.groupby(group_cols):
        
        mutation_description = "chr%s %s" % (
            chromosome,
            gene_mutation_description(pos, ref, alt),  
        )
        key = (mutation_description, transcript_id) 

        def skip(msg, *args):
            msg = msg % args
            logging.info(
                "Skipping %s on %s: %s" , 
                    mutation_description, 
                    transcript_id, 
                    msg)
            variant_report[key] = msg

        def error(msg, *args):
            msg = msg % args
            logging.warning(
                "Error in %s on %s: %s" , 
                    mutation_description, 
                    transcript_id, 
                    msg)
            variant_report[key] = msg

        def success(row):
            new_rows.append(row)
            msg = "SUCCESS: Gene = %s, Mutation = %s" % \
                (row['Gene'], row['PeptideMutationInfo'])
            variant_report[key] = msg

        
        if chromosome.upper().startswith("M"):
            skip("Mitochondrial DNA is insane, don't even bother")
            continue 
        elif ref == alt:
            skip("Not a variant, since ref %s matches alt %s", ref, alt)
            continue 

        padding = max_peptide_length - 1 
        if transcript_id:
            seq, start, stop, annot = \
                peptide_from_transcript_variant(
                    transcript_id, pos, ref, alt,
                    padding = padding)
        else:
            error("Skipping due to invalid transcript ID")
            continue

        if not seq:
            error(annot)
        else:
            if any(s.startswith(seq) for s in seen_source_sequences):
                skip("Already seen sequence starting with %s", seq)
                continue
            else:
                seen_source_sequences.add(seq)

            if '*' in seq:
                error(
                    "Found stop codon in peptide %s from transcript_id %s",
                    seq,
                    transcript_id)
            elif not is_valid_peptide(seq):
                error(
                    "Invalid peptide sequence for transcript_id %s: %s",
                    transcript_id, 
                    seq)
            elif len(seq) < min_peptide_length:
                skip(
                    "Truncated peptide (len %d) too short for transcript %s", 
                    len(seq), 
                    transcript_id)
            else:
                row = deepcopy(group.irow(0))
                row['SourceSequence'] = seq
                row['MutationStart'] = start
                row['MutationEnd'] = stop
                gene_mutation_info = "chr%s %s" % (
                    chromosome, 
                    gene_mutation_description(pos, ref, alt)) 
                row['GeneMutationInfo'] = gene_mutation_info
                row['PeptideMutationInfo'] = annot
                try:
                    gene = gene_names.transcript_id_to_gene_name(transcript_id)
                except:
                    gene = gene_names.transcript_id_to_gene_id(transcript_id)

                row['Gene'] = gene
                success(row)

    assert len(new_rows) > 0, "No mutations!"
    peptides = pd.DataFrame.from_records(new_rows)
    peptides['GeneInfo'] = peptides['info']
    peptides['TranscriptId'] = peptides['stable_id_transcript'] 

    transcripts_df = transcripts_df.merge(peptides)
    logging.info(
        "Generated %d peptides from %s",
        len(transcripts_df),
        patient_id
    )

    # drop verbose or uninteresting columns from VCF
    dumb_fields = (
        'description_gene', 
        'filter', 
        'qual', 
        'id', 
        'name', 
        'info', 
        'stable_id_transcript'
    )
    for dumb_field in dumb_fields:
        if dumb_field in transcripts_df.columns:
            transcripts_df = transcripts_df.drop(dumb_field, axis = 1)
    
    return transcripts_df, vcf_df, variant_report 

def load_variants(input_filename):
    """
    Read the input file into a DataFrame containing (at least) 
    the basic columns of a VCF:
        - chr 
        - pos 
        - ref 
        - alt 
    """
    # VCF and MAF files give us the raw mutations in genomic coordinates
    if input_filename.endswith(".vcf"):
        vcf_df = load_vcf(input_filename)
    elif input_filename.endswith(".maf"):
        maf_df = load_maf(input_filename)
        vcf_df = maf_to_vcf(maf_df)
    elif input_filename.endswith("tab"):
        tab_df = pd.read_csv(input_filename, sep='\t', header=0)
        vcf_df = tab_to_vcf(tab_df)
    else:
        assert False, "Unrecognized file type %s" % input_filename
    return vcf_df 

def load_file(input_filename, min_peptide_length=9, max_peptide_length=31):
    """
    Load mutatated peptides from FASTA, VCF, or MAF file. 
    For the latter two formats, expand their variants across all 
    annotated transcripts.  

    Parameters
    --------

    input_filename : str 

    min_peptide_length : int 

    max_peptide_length : int 

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
        - GeneMutationInfo : original genetic variant e.g. chr3 g.484899 C>T
        - PeptideMutationInfo : annotation e.g. V600E
    """

    if input_filename.endswith(".fasta") \
            or input_filename.endswith(".fa"):
        return load_fasta(input_filename, peptide_length = max_peptide_length)

    vcf_df = load_variants(input_filename)
    return expand_transcripts(
        vcf_df,
        input_filename, 
        min_peptide_length = min_peptide_length,
        max_peptide_length = max_peptide_length)
