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


import pandas as pd
from transcript_metadata import download_transcript_metadata

def cached_property(fn):
    field_name = '_cached_%s' % fn.__name__
    @property
    def _cached_property(self):
        if not hasattr(self, field_name):
            setattr(self, field_name, fn(self))
        return getattr(self, field_name)
    return _cached_property


class EnsemblAnnotationData(object):

    """
    Singleton class which allows for lazy loading of
    exon/transcript annotations
    """

    def __init__(self):
        pass

    @cached_property
    def transcript_metadata_path(self):
        return download_transcript_metadata()
        

    @cached_property
    def exons_dataframe(self):
        """
        Dataframe containing exon data with columns:
            - stable_id_transcript  
            - stable_id_gene
            - seq_region_start_gene
            - seq_region_end_gene
            - seq_region_strand_gene
            - seq_region_start_transcript
            - seq_region_end_transcript
            - seq_start
            - start_exon_id
            - seq_end
            - end_exon_id 
            - stable_id_translation  
            - stable_id_exon 
            - exon_id 
            - rank
        """
        path = self.transcript_metadata_path
        return pd.read_csv(path, sep='\t', low_memory = False)
        
    @cached_property
    def transcript_exons_groups(self):
        """
        Mapping from stable_id_transcript to groupby objects of the exons dataframe
        """
        exons_df = self.exons_dataframe
        return exons_df.groupby('stable_id_transcript', sort = False)
    
    @cached_property
    def transcript_exons_dict(self):
        """
        Mapping from stable_id_transcript to groupby objects of the exons dataframe
        """
        transcript_groups = self.transcript_exons_groups
        return dict((k,v) for k,v in transcript_groups)
        

    @cached_property
    def start_exons_dataframe(self):
        """
        Subset of the exon table but only for first exon in each gene
        """
        all_exons = self.exons_dataframe
        mask = all_exons['exon_id'] == all_exons['start_exon_id']
        return all_exons[mask]
        
    @cached_property
    def transcripts_dataframe(self):
        """
        Subset columns for transcript data only
        """
        transcript_cols = [
            'name', 
            'stable_id_gene', 
            'description_gene',
            'seq_region_start_gene', 
            'seq_region_end_gene',
            'seq_region_strand_gene', 
            'stable_id_transcript', 
            'seq_region_start_transcript', 
            'seq_region_end_transcript'
        ]
        return self.exons_dataframe[transcript_cols].drop_duplicates()

    @cached_property
    def gene_dataframe(self):
        """
        Subset columns for gene data only
        """
        gene_cols = [
            'name', 
            'stable_id_gene', 
            'description_gene',
            'seq_region_start_gene', 
            'seq_region_end_gene'
        ]
        return self.transcripts_dataframe[gene_cols].drop_duplicates()
        