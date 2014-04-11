## Immunotherapy Pipeline


### Usage 
#### From peptide string
```sh
python cancer_pipeline.py --string EDLTVKIGDFGLATEKSRWSGSHQFEQLS --hla HLA-B*35:01
```

#### From input file
```sh
python cancer_pipeline.py --input <.maf, .vcf, .eff.vcf, .fasta> --hla-file <allele-file> 

```

If  VCF file is input, it will be annotated w/ gene information and Ensembl transcripts, otherwise input an annotated SnpEFF VCF file.

If you don't either an `hla` or `hla-file` then the default "HLA:A*02:01" will be used. 

### Options:
* `--input`: Input file containing somatic variants (one of MAF, VCF, or SnpEff annotated VCF)
* `--string`: String of amino acids, with optional annotation of a mutated region between underscores (i.e. QYS\_LL\_Q)
* `--hla`: Comma separated list of HLA alleles. 
* `--hla-file`: Text file containing one HLA allele per line. 
* `--peptide-length`: Length of vaccine peptide (window around mutation, default 31)
* `--epitopes-output`: CSV output file for dataframe containing scored epitopes
* `--peptides-output`: CSV output file for dataframe containing scored vaccine
* `--print-epitopes`: Print dataframe with epitope scores
* `--print-peptides`: Print dataframe with vaccine peptide scores
* `--html-report`: Path to HTML report containing scored peptides and epitopes
* `--skip-mhc`: Don't predict MHC binding


### Requirements
First you'll have to install the epitopes package.  You can download that from: [https://github.com/hammerlab/epitopes](https://github.com/hammerlab/epitopes/)

```sh
git clone https://github.com/hammerlab/epitopes.git
cd epitopes
python setup.py install
```

Additionally, you'll need to install:

* [nose](https://nose.readthedocs.org/en/latest/)
* [NumPy](http://www.numpy.org/)
* [pandas](http://pandas.pydata.org/)
* [BioPython](http://biopython.org/wiki/Main_Page)
