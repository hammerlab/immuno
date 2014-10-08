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
    """
    Run the given function `fn` the first time this property is accessed,
    save the result.
    """
    field_name = '_cached_%s' % fn.__name__
    def _cached_property(self):
        if not hasattr(self, field_name):
            setattr(self, field_name, fn(self))
        return getattr(self, field_name)
    _cached_property.__name__ = fn.__name__
    return property(_cached_property)


class EnsemblAnnotationData(object):

    """
    Singleton class which allows for lazy loading of
    exon/transcript annotations
    """

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
    def transcript_exons_dict(self):
        """
        Mapping from stable_id_transcript to DataFrame of that
        transcript's exons
        """
        class DelayedDict(object):
            def __init__(self, exons_df):
                self.exons_df = exons_df
                self.groupby_object = exons_df.groupby(['stable_id_transcript'])
                self.group_indices = self.groupby_object.groups
                self._d = {}

            def __getitem__(self, transcript_id):
                # have we retrieved this transcript_id before?
                if transcript_id in self._d:
                    return self._d[transcript_id]
                # is this transcript in the original exons DataFrame?
                if transcript_id not in self.group_indices:
                    raise KeyError(transcript_id)
                # which rows had the given transcript id?
                indices = self.group_indices[transcript_id]
                subset = self.exons_df.ix[indices]
                self._d[transcript_id] = subset
                return subset

            def __contains__(self, transcript_id):
                return transcript_id in self.group_indices

        return DelayedDict(self.exons_dataframe)

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
