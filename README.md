## Immunotherapy Pipeline


### Usage

#### From peptide string
```sh
python cancer_pipeline.py --string EDLTVKIGDFGLATEKSRWSGSHQFEQLS
```

#### From input file string
```sh
python cancer_pipeline.py --input <.maf, .vcf, .eff.vcf, .fasta> --output <output_file> 

```

If  VCF file is input, it will be annotated w/ gene information and Ensembl transcripts, otherwise input an annotated SnpEFF vcf file.

If no output_file is specific the results will print to stdout


### Adding a new scorer

```python

from pipeline import PipelineElement

class NewScorer(PipelineElement):
  def __init__(name):

  def apply(self, data):
    return transform(data)

```
Then in `cancer_pipeline.py`

```python
pipeline.add_scorer(NewScorer(name="new scoring mech"))

```


