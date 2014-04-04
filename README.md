## Immunotherapy Pipeline


### Usage

#### From peptide string
```sh
python cancer_pipeline.py --string EDLTVKIGDFGLATEKSRWSGSHQFEQLS
```

#### From input file
```sh
python cancer_pipeline.py --input <.maf, .vcf, .eff.vcf, .fasta> --output <output_file> 

```

If  VCF file is input, it will be annotated w/ gene information and Ensembl transcripts, otherwise input an annotated SnpEFF vcf file.

If no output_file is specific the results will print to stdout

