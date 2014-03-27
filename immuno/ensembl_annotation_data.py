
from ensembl_download import download_transcript_metadata

class EnsemblAnnotationData(object):

    """
    Singleton class which allows for lazy loading of
    exon/transcript annotations
    """

    def __init__(self):
        pass

    @property
    def transcript_metadata_path(self):
        if not hasattr(self, '_transcript_metadata_path'):
            self._transcript_metadata_path = download_transcript_metadata()
        return self._transcript_metadata_path

    @property
    def exon_data(self):
        """
        Dataframe containing exon data
        """
        if not hasattr(self, '_exon_data'):
            path = self.transcript_metadata_path
            self._exon_data = pd.read_csv(path, sep='\t')
        return self._exon_data

    @property
    def transcript_data(self):
        """
        Subset columns for transcript data only
        """
        if not hasattr(self, '_transcript_data'):
            transcript_cols = [
                'name', 'stable_id_gene', 'description_gene',
                'seq_region_start_gene', 'seq_region_end_gene',
                'stable_id_transcript', 'seq_region_start_transcript',
                'seq_region_end_transcript'
            ]
            self._transcript_data = self.exon_data[transcript_cols]
        return self._transcript_data

    @property
    def gene_data(self):
        """
        Subset columns for gene data only
        """
        if not hasattr(self, '_gene_data'):
            gene_cols = [
                'name', 'stable_id_gene', 'description_gene',
                'seq_region_start_gene', 'seq_region_end_gene'
            ]
            self._gene_data = self.transcript_data[gene_cols].drop_duplicates()
