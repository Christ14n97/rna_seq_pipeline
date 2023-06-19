# RNA-seq pipeline

This repository contains essential files to cover RNA-seq pipeline. In particular from FASTQ files mapping with STAR, to count matrix with R. This pipeline depends on `STAR`,
`samtools`, `R-4.3`, `BiocManager=3.17`, `GenomicAlignments` & `GenomicFeatures`.

The Makefile under the directory contains all the instructions of the pipeline. However, files in `bin/` directory are needed to compute the count matrix.

Pipeline:
- First it performs mapping of the FASTQ using `STAR` to a reference genome.
- Secondly, computes statistics with `samtools`.
- Ultimately, generates the count matrix using R scripts in `bin/` folder.

## How to use it
Assuming that you have your reference genome ready to be used for mapping purposes, these are the steps to follow:

`STEP 1`: place Make file in the folder containing FASTQ.gz files and navigate to it

```bash
cp Makefile path/to/fastq.gz_files/.
cd path/to/fastq.gz_files
```

`STEP 2`: create an array

```bash
files=($(basename --suffix=.fastq.gz -- *.fastq.gz))
```

`STEP 3`: loop over the array and do `make`

```bash
for n in ${files[@]}; do make BASE_NAME=$n \
  GENOME_DIR=../reference_genomes/dicty_myco_merged_star_index \
  QUANT_MODE='gene' \
  BAM_FILE=../mapped/${n}.DdMm.star/Aligned.sortedByCoord.out.bam ;
  done
```
