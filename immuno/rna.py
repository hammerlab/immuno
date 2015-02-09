import pandas as pd
from varcode import CodingSequenceMutation

# default column names from cufflinks tracking files
# for gene and isoform expression levels
STATUS_COLUMN = "FPKM_status"
ID_COLUMN = "tracking_id"
FPKM_COLUMN = "FPKM"
LOCUS_COLUMN = "locus"
GENE_NAMES_COLUMN = "gene_short_name"

def load_cufflinks_tracking_file(
        filename,
        id_column=ID_COLUMN_NAME,
        fpkm_column=FPKM_COLUMN_NAME,
        status_column=STATUS_COLUMN_NAME,
        locus_column=LOCUS_COLUMN_NAME,
        gene_names_column=GENE_NAMES_COLUMN_NAME,
        drop_failed_fpkm=True,
        drop_lowdata_fpkm=True,
        drop_hidata_fpkm=True,
        drop_nonchromosomal_loci=True
        drop_novel=False):
    """
    Loads a Cufflinks tracking file, which contains expression levels
    (in FPKM: Fragments Per Kilobase of transcript per Million fragments)
    for transcript isoforms or whole genes. These transcripts/genes may be
    previously known (in which case they have an Ensembl ID) or a novel
    assembly from the RNA-Seq data (in which case their IDs look like "CUFF.1")

    Parameters
    ----------

    filename : str
        Filename of tracking file e.g. "genes.tracking_fpkm"

    id_column : str, optional

    fpkm_column : str, optional

    status_column : str, optional
        Name of column which indicates the FPKM estimate status. The column
        name is typically "FPKM_status". Possible contained within this column
        will be OK, FAIL, LOWDATA, HIDATA.

    locus_column : str, optional

    gene_names_column : str, optional

    drop_failed : bool, optional
        Drop rows whose FPKM status is "FAIL" (default = True)

    drop_lowdata : bool, optional
        Drop rows whose FPKM status is "LOWDATA" (default = True)

    drop_hidata : bool, optional
        Drop rows whose FPKM status is "HIDATA" (default=True)

    drop_nonchromosomal_loci : bool, optional
        Drop rows whose location isn't on a canonical chromosome
        i.e. doesn't start with "chr" (default=True)

    drop_novel : bool, optional
        Drop genes or isoforms that aren't found in Ensembl (default = False)

    Returns DataFrame with columns:
        id : str
        novel : bool
        fpkm : float
        chr : str
        start : int
        end : int
        gene_names : str list
    """
    df = pd.read_csv(filename, sep='\t')

    if drop_failed:
        status = df[status_column]
        fail_mask = status == "FAIL"
        logging.info("Dropping %d/%d failed entries from %s" % (
            fail_mask.sum(), len(df), filename))
        df = df[~fail_mask]

    if drop_nonchromosomal_loci:
        loci = df[locus_column]
        chromosomal_loci = loci.str.startswith("chr")
        logging.info("Dropping %d/%d non-chromosomal loci from %s" % (
            (~chromosomal_loci).sum(), len(df), filename))
        df = df[chromosomal_loci]

    if len(df) == 0:
        raise ValueError("Empty FPKM tracking file: %s" % filename)

    ids = df[id_column]
    known = ids.str.startswith("ENS")
    if known.sum() == 0:
        raise ValueError("No Ensembl IDs found in %s" % filename)

    loci = df[locus_column]
    # capture all characters after 'chr' but before ':'
    chromosomes = loci.str.extract("chr([^:]*):.*")
    # capture all characters after e.g. 'chr1:', which look like '132-394'
    ranges = loci.str.extract("chr[^:]*:(.*)")
    # capture all numbers before the dash
    starts = ranges.str.extract("(\d*)-\d*").astype(int)
    # capture all numbers after the dash
    ends = ranges.str.extract("\d*-(\d*)")

    # gene names are given either as "-" or a comma separated list
    # e.g. "BRAF1,PFAM2"
    gene_names_strings = df[gene_names_column]
    gene_names_strings[gene_names_strings == "-"] = ""
    # split each entry into a list of zero or more strings
    gene_names_lists = gene_names_strings.str.split(",")

    return pd.DataFrame({
        'id' : df[id_column_name],
        'novel' : ~known,
        'fpkm' : df[fpkm_column_name],
        'chr' : chromosomes,
        'start' : starts,
        'end' : ends,
        "gene_names" : gene_names_lists
    })

def aggregate_gene_expression_levels(gene_fpkm_df, ensembl):
    """
    Create a dictionary mapping gene IDs to expression levels.

    Parameters
    ----------
    gene_fpkm_df : DataFrame
        DataFrame with columns:
            - 'id'
            - 'novel'
            - 'fpkm'
            - 'chr'
            - 'start'
            - 'end'
            - 'gene_names'
        IDs can be either from Ensembl or Cufflinks novel genes

    ensembl : pyensembl.EnsemblRelease
        Ensembl annotation database used to look up genes overlapping
        chromosomal locations.

    We can't just use the values from a Cufflinks gene.fpkm_tracking file and
    instead have to add up contributions from possibly multiple novel gene IDs
    all of which vouch for expression of annotated genes. For unclear reasons,
    Cufflinks sometimes assigns a novel gene ID (e.g. "CUFF.1") to regions of
    the genome which overlap an annotated gene.

    Example entry from a Cufflinks gene tracking file:

        CUFF.29167  -   -   CUFF.29167  BRAF    -   chr7:140424942-140624564

    This region overlaps BRAF entirely and thus the FPKM value of this
    novel gene should be assigned to BRAF's gene ID (ENSG00000157764). In the
    case that a novel gene entry overlaps multiple Ensembl genes, divide that
    entry's contribution by the number of overlapped genes.
    """

    # first drop any entries with FPKM == 0, since they don't
    # contribute to the sum of any gene

    gene_fpkm_df = gene_fpkm_df[gene_fpkm_df.fpkm > 0]

    # start with expression levels of known Ensembl genes
    mask = gene_fpkm_df.known
    known_df = gene_fpkm_df[mask]
    result_dict = dict(zip(known_df.id, fpkm))

    novel_df = gene_fpkm_df[~mask]
    for _, novel_row in novel_df.iterrows()
        # find genes overlapping the chromosomal positions spanned
        # by the novel gene constructed by Cufflinks
        overlapping_genes = ensembl.genes_at_locus(
            novel_row.chr,
            novel_row.start,
            novel_row.end)
        # some genes may be on the wrong strand so only use those
        # whose name was in the Cufflinks file
        n_matching_genes = [
            gene
            for gene in overlapping_genes
            if gene.name in novel_row.gene_names
        ]
        n_matching_genes = len(matching_genes)
        for gene in matching_genes:
            old_fpkm = result_dict.get(gene.id, 0.0)
            # split FPKM expression value across all matching genes
            # overlapping this locus, since we don't know which one to assign
            # expression and don't want to inflate loci that span lots of
            # genes
            result_dict[gene.id] = old_fpkm + novel_row.fpkm / n_matching_genes

    return result_dict

def expressed_gene_ids(
        gene_fpkm_filename,
        gene_expression_threshold,
        ensembl):


    # gene IDs with sufficiently high expression level to be considered
    # for inclusion in the vaccine
    return {
        gene_id
        for (gene_id, fpkm) in gene_fpkms_dict.iteritems()
        if fpkm >= gene_expression_threshold
    }

def choose_principal_transcripts(
        variant_collection
        gene_fpkm_df,
        gene_expression_threshold,
        transcript_fpkm_df):

    # dictionary whose keys are Ensembl gene IDs and values are FPKM values
    gene_fpkm_dict = aggregate_ensembl_gene_expression_levels(
        gene_fpkms_df
        ensembl=ensembl)

    transcript_fpkm_dict = dict(
        zip(transcript_fpkm_df.id, transcript_fpkm_df.fpkm))

    principal_transcript_effects = []
    for variant_effect in variant_collection.variant_effects():

        # mapping from transcript ID to pair (gene fpkm, transcript fpkm)
        # we use this to first look at genes of high expression and then
        # choose their most highly expressed transcript
        transcript_expression_levels = {}
        transcript_effect_list = []
        for gene_id, transcript_effects in \
                variant_effect.gene_transcript_effects.iteritems():
            gene_fpkm = gene_fpkm_dict.get(gene.id)

            if gene_fpkm <= gene_expression_threshold:
                continue

            for transcript_effect in transcript_effects:
                if not isinstance(transcript_effect, CodingSequenceMutation):
                    continue

                transcript_id = transcript_effect.transcript.id
                transcript_fpkm = transcript_fpkm_dict.get(
                    transcript_id, 0.0)

                if transcript_fpkm <= 0:
                    continue

                fpkm_pair = (gene_fpkm, transcript_fpkm)
                transcript_expression_levels[transcript_id] = fpkm_pair
                transcript_effect_list.append(transcript_effect)

        def key(transcript_effect):
            return transcript_expression_levels[transcript_effect.transcript.id]

        transcript_effect_list.sort(key=key, reverse=True)

        if len(transcript_effect_list) > 0:
            best_transcript_effect = transcript_effect_list[0]
            principal_transcript_effects.append(best_transcript_effect)

    return principal_transcript_effects

