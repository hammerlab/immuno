## Immunotherapy Pipeline

### Usage
First you'll have to setup the epitopes package.  You can download that here: https://github.com/hammerlab/epitopes

```sh
git clone https://github.com/hammerlab/epitopes.git
cd epitopes
python setup.py install
```

#### From peptide string
```sh
python cancer_pipeline.py --string EDLTVKIGDFGLATEKSRWSGSHQFEQLS
```

#### From input file
```sh
python cancer_pipeline.py --input <.maf, .vcf, .eff.vcf, .fasta> --allele_file <allele-file> --output <output_file> 

```

If  VCF file is input, it will be annotated w/ gene information and Ensembl transcripts, otherwise input an annotated SnpEFF vcf file.

If no output_file is specific the results will print to stdout

If you don't input an allele_file it will default to HLA:A*02:01 and you can also score a single peptide sequence with


