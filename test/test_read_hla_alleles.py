from cStringIO import StringIO
from immuno.hla_file import read_hla_file

def test_read_hla():
    # contains two C identical C alleles
	hla_filename = 'data/SKCM.hla'
	alleles = read_hla_file(hla_filename)
	print alleles
	assert alleles == \
		[
		  'HLA-A*01:01',
		  'HLA-A*03:01',
		  'HLA-B*08:01',
		  'HLA-B*07:02',
		  'HLA-C*07:02',
		  'HLA-C*07:02'
		]
	assert len(set(alleles)) == 5