import datetime
import numpy as np

def build_html_report(scored_data):
    scored_data = scored_data.sort(columns=('combined_score',), ascending=False)
    table = scored_data.to_html(
        index=False,
        na_rep="-",
        columns = [
            'SourceSequence', 'info', 'ref', 'alt', 'pos',
            'Epitope', 'EpitopeStart', 'EpitopeEnd',
            'percentile_rank',
            'ann_rank',
            'ann_ic50',
            'immunogenicity',
            'mhc_score', 'imm_score',
            'combined_score'
        ])

    # take each source sequence and shade its amino acid letters
    # based on the average score of each epitope containing that letter
    seq_divs = []
    seq_scores = []

    for seq in scored_data.SourceSequence.unique():
        scores = np.zeros(len(seq), dtype=float)
        imm_scores = np.zeros(len(seq), dtype=float)
        mhc_scores = np.zeros(len(seq), dtype=float)
        score_counts = np.ones(len(seq), dtype=int)
        rowslice = scored_data[scored_data.SourceSequence == seq]
        gene_info = None
        for _, row in rowslice.iterrows():
            gene_info = row['info']
            start = int(row['EpitopeStart'] - 1)
            assert start >= 0, start
            stop = int(row['EpitopeEnd'])
            assert stop > start, stop
            scores[start:stop] += row['combined_score']
            imm_scores[start:stop] += row['immunogenicity']
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
        for i in xrange(len(seq)):
            letter = seq[i]
            score = scores[i]
            letter_td = "<td>%s</td>" % letter
            letters.append(letter_td)

            imm = imm_scores[i]
            mhc = mhc_scores[i]
            maxval = 256
            mhc_intensity = int(mhc**2*maxval)
            mhc_rgb = "rgb(%d, %d, %d)" % \
                (mhc_intensity/3, mhc_intensity, mhc_intensity/2)
            imm_intensity =  int(imm**2*maxval)
            imm_rgb = "rgb(%d, %d, %d)" % \
                (imm_intensity, imm_intensity/2, imm_intensity/3)


            color_cell = \
            """
            <td style="background-color: %s;">&nbsp;</td>
            """
            mhc_color_cell = color_cell %  mhc_rgb
            imm_color_cell = color_cell % imm_rgb
            imm_colors.append(imm_color_cell)
            mhc_colors.append(mhc_color_cell)

        median_score = np.median(scores)
        letters_cols = '\n\t'.join(letters)
        mhc_color_cols = '\n\t'.join(mhc_colors)
        imm_color_cols = '\n\t'.join(imm_colors)
        colored_letters_table = \
            """

            <table>
            <center>
            <tr>
            <td style='background-color: rgb(190,190,190);'>Sequence</td>
            %s
            </tr>
            <tr>

            <td style='background-color: rgb(190,190,190);'>MHC Binding</td>
            %s
            </tr>
            <tr>

            <td style='background-color: rgb(190,190,190);'>Immunogenicity</td>
            %s
            </tr>
            </center>
            </table>
            """ % (letters_cols, mhc_color_cols, imm_color_cols)

        div = """
            <div
                style="border-bottom: 1px solid gray; margin-bottom: 1em;"
                class="seq">
            <h3>Median Epitope Score = %0.4f (%s)</h3>
            %s
            <br>
            </div>
            """ % (median_score, gene_info, colored_letters_table)
        seq_divs.append(div)
        seq_scores.append(median_score)

    seq_order = reversed(np.argsort(seq_scores))
    seq_divs_html = "\n".join(seq_divs[i] for i in seq_order)

    page = """
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
            }
            th { background-color: rgb(90, 190, 240); }
        </style>
        <head><title>Immune Pipeline Results (%s)</title></head>
        <body>
        <h2>Mutation Regions</h2>
        %s
        <hr>
        <h2>Sorted Scores Results</h2>
        <center>
        %s
        </center>
        </body>
        </html>
    """ % (datetime.date.today(), seq_divs_html, table)
    return page
