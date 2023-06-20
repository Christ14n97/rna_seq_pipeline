# RNA-seq pipeline

This repository contains essential files to cover RNA-seq pipeline. In particular from FASTQ files mapping with STAR, to count matrix with R. This pipeline depends on `STAR`,
`samtools`, `R-4.3`, `BiocManager=3.17`, `GenomicAlignments` & `GenomicFeatures`.

The Makefile under the directory contains all the instructions of the pipeline. However, files in `bin/` directory are needed to compute the count matrix.

Pipeline:
- First it performs mapping of the FASTQ using `STAR` to a reference genome.
- Secondly, computes statistics with `samtools`.
- Ultimately, generates the count matrix using R scripts in `bin/` folder.

## Build docker image

```bash
docker build --progress=plain --cpuset-cpus 8 --rm -f "container/Dockerfile" -t "unigebsp/rnaseq_pipe" "container"
```

## Docker container directory structure

Inside the container the structure is as follows:

```
|-----home/user/pipeline
                |
                |---fastq/ - folder containing .fastq.gz files.
                |---ref/ - folder containing reference genome fasta files.
                |---Makefile - file containing pipeline.
```

## How to use the container

Run the docker container:

```bash
docker run -it --rm -v $PWD/test/reference_genomes/:/home/user/pipeline/ref/ \
       -v $PWD/test/fastq/:/home/user/pipeline/fastq/ unigebsp/rnaseq_pipe
```

In this case we mount:

1. `$PWD/test/reference_genomes/`, our local folder containing a reference genome of our interest to docker container folder `/home/user/pipeline/ref/` to work with it.
2. `$PWD/test/fastq/`, our local folder containing our fastq.gz files to docker container folder `home/user/pipeline/fastq/` to work with it.