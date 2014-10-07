from immuno.maf import load_maf
from immuno.load_file import load_variants

from nose.tools import eq_

def test_load_maf():
	filename = 'data/SKCM.maf'
	maf_df = load_maf(filename)
	assert len(maf_df) > 0
	cols = set(maf_df.columns)
	assert 'Tumor_Sample_Barcode' in cols
	assert 'Reference_Allele' in cols
	assert 'Tumor_Seq_Allele1' in cols
	assert 'Tumor_Seq_Allele2' in cols
	assert 'Chromosome' in cols
	assert 'Start_Position' in cols
	assert 'End_Position' in cols
	df = load_variants(filename)
	assert len(maf_df) == len(df)