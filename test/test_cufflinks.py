"""
test_cufflinks : Test that we can correctly load Cufflinks tracking files which
contain the estimated expression levels of genes and isoforms (computed from
RNA-Seq reads).
"""
from immuno.rna import load_cufflinks_tracking_file
from nose.tools import eq_


def test_load_cufflinks_genes():
    genes_df = load_cufflinks_tracking_file("data/genes.fpkm_tracking")
    gene_ids = set(genes_df.id)
    expected_gene_ids = {
        "ENSG00000240361",
        "ENSG00000268020",
        "ENSG00000186092",
        "ENSG00000269308",
        "CUFF.1",
        "CUFF.2",
        "CUFF.3",
        "CUFF.4",
        "CUFF.5"
    }
    eq_(gene_ids, expected_gene_ids)

def test_load_cufflinks_genes_drop_novel():
    genes_df = load_cufflinks_tracking_file(
        "data/genes.fpkm_tracking",
        drop_novel=True)
    gene_ids = set(genes_df.id)
    expected_gene_ids = {
        "ENSG00000240361",
        "ENSG00000268020",
        "ENSG00000186092",
        "ENSG00000269308",
    }
    eq_(gene_ids, expected_gene_ids)


def test_load_cufflinks_isoforms():
    transcripts_df = load_cufflinks_tracking_file("data/isoforms.fpkm_tracking")
    transcript_ids = set(transcripts_df.id)
    expected_transcript_ids = {
        "ENST00000492842",
        "ENST00000594647",
        "ENST00000335137",
        "ENST00000417324",
        "ENST00000461467",
        "ENST00000518655",
        "CUFF.7604.1",
    }
    eq_(transcript_ids, expected_transcript_ids)

def test_load_cufflinks_isoforms_drop_novel():
    transcripts_df = load_cufflinks_tracking_file(
        "data/isoforms.fpkm_tracking", drop_novel=True)
    transcript_ids = set(transcripts_df.id)
    expected_transcript_ids = {
        "ENST00000492842",
        "ENST00000594647",
        "ENST00000335137",
        "ENST00000417324",
        "ENST00000461467",
        "ENST00000518655",
    }
    eq_(transcript_ids, expected_transcript_ids)