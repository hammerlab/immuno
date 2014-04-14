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

import datetime
import numpy as np

page_template = \
"""
<html>
<style>
    body { padding: 2em; font-family: sans-serif; }
    table { padding: 0em; border: 0px solid black; }
    table, td, th
    {

        text-align:center;

    }
    td, th {
        padding: 0.2em;
        border:1px solid gray;
    }

    .seq td {
        height: 2em;
        width:1.5em;
        background-color: rgb(220,220,220);
        padding: 0em;
        color: rgb(0,0,0);
    }

    .seq td.wildtype { color: rgb(30, 30, 30); }
    .seq td.mutant  {
        color: rgb(190, 30, 30);
        background-color: rgb(220, 240, 220);
        font-weight: bold;
    }

    .seq td.near_mutant {
        color: rgb(100, 30, 30);
        background-color: rgb(235, 235, 220);
    }

    .seq td.row_title {
        background-color: rgb(190,190,190);
    }

    th { background-color: rgb(90, 190, 240); }
</style>
<head><title>Immune Pipeline Results (%s)</title></head>
<body>
<h2>Vaccine Peptides</h2>
%s
<hr>

<h2>All Peptide Scores</h2>
<center>
%s
</center> 
<hr>

<h2>All Epitope Scores</h2>
<center>
%s
</center>
</body>
</html>
"""

table_template = \
"""
<table>
<center>
<tr>
<td  class='row_title'>Sequence</td>
%s
</tr>
<tr>

<td class='row_title'>MHC Binding</td>
%s
</tr>
<tr>

<td class='row_title'>Immunogenicity</td>
%s
</tr>
</center>
</table>
"""

peptide_div_template = \
"""
<div
    style="border-bottom: 1px solid gray; margin-bottom: 1em;"
    class="seq">
<h3>Peptide Score = %0.4f (%s, transcript=%s)</h3>
%s
<br>
</div>
"""

color_cell = \
"""
<td style="background-color: %s;">&nbsp;</td>
"""

def build_html_report(scored_epitopes, scored_peptides, unique_transcripts = True):
    scored_epitopes = scored_epitopes.sort(
        columns=('combined_score',), ascending=False)

    # take each source sequence and shade its amino acid letters
    # based on the average score of each epitope containing that letter
    peptide_divs = []
    peptide_scores = []

    # use size of longest possible epitope to figure out if a letter
    # is 'near' a mutation
    max_epitope_len = scored_epitopes.Epitope.str.len().max()

    seen_transcript_ids = set([])

    group_cols = [
        "Peptide",
        "PeptideStart",
        "PeptideEnd",
        "PeptideMutationStart",
        "PeptideMutationEnd",
        "GeneInfo",
        "TranscriptId",
        "Score",
        "SourceSequence"]

    for (
            peptide,
            peptide_start, peptide_end,
            mut_start, mut_end,
            gene_info,
            transcript_id,
            peptide_score, src_seq
        ), _ in scored_peptides.groupby(group_cols):

        if unique_transcripts:
            if transcript_id in seen_transcript_ids:
                continue
            else:
                seen_transcript_ids.add(transcript_id)  
                      
        n = len(peptide)

        scores = np.zeros(n, dtype=float)
        imm_scores = np.zeros(n, dtype=float)
        mhc_scores = np.zeros(n, dtype=float)
        score_counts = np.ones(n, dtype=int)

        mask = (scored_epitopes.SourceSequence == src_seq)
        mask &= scored_epitopes.EpitopeStart >= peptide_start
        mask &= scored_epitopes.EpitopeEnd <= peptide_end

        rowslice = scored_epitopes[mask]

        for _, row in rowslice.iterrows():
            epitope_start = int(row['EpitopeStart'])
            assert epitope_start >= 0, \
                "Expected epitope start %d >= 0" %  epitope_start
            assert epitope_start >= peptide_start, \
                "Expected epitope start %d >= peptide start %d" % \
                (epitope_start, peptide_start)
            epitope_end = int(row['EpitopeEnd'])
            assert epitope_end > epitope_start, epitope_end
            assert epitope_end <= peptide_end, epitope_end

            start = epitope_start - peptide_start
            stop = epitope_end - peptide_end

            scores[start:stop] += row['combined_score']
            if 'immunogenicity' in row:
                imm_scores[start:stop] += row['immunogenicity']
            if 'percentile_rank' in row:
                mhc_scores[start:stop] += (100 - row['percentile_rank']) / 100.0
            score_counts[start:stop] += 1

        # default background for all letters of the sequence is gray
        # but make it more red as the score gets higher
        letters = []
        colors = []
        imm_colors = []
        mhc_colors = []
        scores /= score_counts
        imm_scores /= score_counts
        mhc_scores /= score_counts
        for i in xrange(n):
            letter = peptide[i]
            score = scores[i]
            # is the amino acid at position i inside a mutated region?
            if i >= mut_start and i < mut_end:
                letter_class = 'mutant'
            # if not, is it at least within an epitope's length near the
            # mutated region?
            elif (i < mut_start and mut_start - i < max_epitope_len) or \
                    (i >= mut_end and i - mut_end < max_epitope_len - 1):
                letter_class = 'near_mutant'
            else:
                letter_class = 'wildtype'

            letter_td = "<td class = '%s'>%s</td>" % (letter_class, letter)
            letters.append(letter_td)

            imm = imm_scores[i]
            mhc = mhc_scores[i]
            maxval = 256
            mhc_intensity = int(mhc**1.5*maxval)
            mhc_rgb = "rgb(%d, %d, %d)" % \
                (mhc_intensity/3, mhc_intensity, mhc_intensity/2)
            imm_intensity =  int(imm**2*maxval)
            imm_rgb = "rgb(%d, %d, %d)" % \
                (imm_intensity, imm_intensity/2, imm_intensity/3)


            mhc_color_cell = color_cell %  mhc_rgb
            imm_color_cell = color_cell % imm_rgb
            imm_colors.append(imm_color_cell)
            mhc_colors.append(mhc_color_cell)

        letters_cols = '\n\t'.join(letters)
        mhc_color_cols = '\n\t'.join(mhc_colors)
        imm_color_cols = '\n\t'.join(imm_colors)
        colored_letters_table = \
            table_template % (letters_cols, mhc_color_cols, imm_color_cols)

        div = peptide_div_template % \
            (peptide_score, gene_info, transcript_id, colored_letters_table)
        peptide_divs.append(div)
        peptide_scores.append(peptide_score)

    seq_order = reversed(np.argsort(peptide_scores))
    peptide_divs_html = "\n".join(peptide_divs[i] for i in seq_order)

    epitope_columns = [
        'Epitope',
        'info',
        'combined_score'
    ]

    optional_columns = [
        'percentile_rank',
        'ann_rank',
        'ann_ic50',
        'immunogenicity',
        'mhc_score',
        'imm_score',
        'stable_id_transcript',
        'ref', 'alt',  'chr', 'pos', "MutationInfo"
    ]
    for col_name in optional_columns:
        if col_name in scored_epitopes:
            epitope_columns.append(col_name)
    peptide_table = scored_peptides.to_html(
        index = False, 
        na_rep = "-",
    )
    
    epitope_table = scored_epitopes.to_html(
        index=False,
        na_rep="-",
        columns = epitope_columns)

    page = page_template % \
        (datetime.date.today(), peptide_divs_html, peptide_table, epitope_table)
    return page
