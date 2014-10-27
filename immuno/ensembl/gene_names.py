
import logging
import pandas as pd
import datacache


# use genenames.org to make a table mapping between HUGO gene names and ensembl
# genen ID with columns
#  'Approved Symbol'
#  'Approved Name'
#  'Ensembl ID(supplied by Ensembl)'
_HUGO_URL = \
"http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=gd_app_name&col=md_ensembl_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit"


def load_hugo_table(_table_cache = [None]):
	if _table_cache[0] is None:
		print "Downloading %s" % _HUGO_URL
		df = datacache.fetch_csv_dataframe(_HUGO_URL, sep="\t")
		_table_cache[0] = df
	return _table_cache[0]

def gene_id_to_name(gene_id, _lookup_cache = [None]):
	if _lookup_cache[0] is None:
		# dataframe mapping between Hugo Names and Ensembl Gene IDs
		hugo = load_hugo_table()
		# Ensembl IDs are sometimes missing, drop those rows
		bad_rows = hugo['Ensembl ID(supplied by Ensembl)'].isnull()
		hugo = hugo[~bad_rows]
		hugo_names = hugo['Approved Symbol']
		gene_ids = hugo['Ensembl ID(supplied by Ensembl)']
		d = dict(zip(gene_ids, hugo_names))
		_lookup_cache[0] = d
	else:
		d = _lookup_cache[0]

	assert isinstance(gene_id, str), \
		"Ensembl gene ID %s must be a string" % gene_id
	assert gene_id.startswith("ENSG"), \
		"Ensembl gene ID must start with characters 'ENSG'" % gene_id
	assert gene_id in d, "Ensembl gene ID %s not found" % gene_id
	return d[gene_id]


_BIOMART_QUERY_TRANSCRIPT_ID_TO_GENE_ID = \
"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1"
	uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
<Attribute name = "ensembl_gene_id" />
<Attribute name = "ensembl_transcript_id" />
</Dataset>
</Query>
""".replace("\n", "")
_BIOMART_URL_TRANSCRIPT_ID_TO_GENE_ID = \
	"http://feb2014.archive.ensembl.org/biomart/martservice/result?query=%s" % \
	_BIOMART_QUERY_TRANSCRIPT_ID_TO_GENE_ID

def transcript_id_to_gene_id(transcript_id, _table_cache = [None]):
	if _table_cache[0] is None:
		print ("Fetching Ensembl ID mappings from BioMart %s"
			) % _BIOMART_URL_TRANSCRIPT_ID_TO_GENE_ID
		biomart_filename = \
			datacache.fetch_file(_BIOMART_URL_TRANSCRIPT_ID_TO_GENE_ID,
				"biomart_transcript_gene.tsv")
		df = pd.read_csv(biomart_filename, sep='\t')
		gene_ids = df['Ensembl Gene ID']
		transcript_ids = df['Ensembl Transcript ID']
		mapping = dict(zip(transcript_ids, gene_ids))
		_table_cache[0] = mapping
	mapping = _table_cache[0]
	return mapping[transcript_id]

def transcript_id_to_gene_name(transcript_id):
	gene_id = transcript_id_to_gene_id(transcript_id)
	return gene_id_to_name(gene_id)

_BIOMART_QUERY_TRANSCRIPT_ID_TO_TRANSCRIPT_NAME = \
"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1"
	uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
<Attribute name = "ensembl_transcript_id" />
<Attribute name = "external_transcript_id" />
</Dataset>
</Query>
""".replace("\n", "")
_BIOMART_URL_TRANSCRIPT_ID_TO_TRANSCRIPT_NAME = \
	"http://feb2014.archive.ensembl.org/biomart/martservice/result?query=%s" % \
	_BIOMART_QUERY_TRANSCRIPT_ID_TO_TRANSCRIPT_NAME

def transcript_id_to_transcript_name(transcript_id, _table_cache = [None]):
	if _table_cache[0] is None:
		print ("Fetching Ensembl ID mappings from BioMart %s"
			) % _BIOMART_URL_TRANSCRIPT_ID_TO_TRANSCRIPT_NAME
		biomart_filename = \
			datacache.fetch_file(_BIOMART_URL_TRANSCRIPT_ID_TO_TRANSCRIPT_NAME,
				"biomart_transcript_name.tsv")
		df = pd.read_csv(biomart_filename, sep='\t')
		transcript_ids = df['Ensembl Transcript ID']
		transcript_names = df['Associated Transcript Name']
		mapping = dict(zip(transcript_ids, transcript_names))
		_table_cache[0] = mapping
	mapping = _table_cache[0]
	return mapping[transcript_id]
