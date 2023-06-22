#!/usr/bin/bash

files=($(ls fastq/ | grep -E "*.fastq.gz$" | sed 's/\.fastq.gz$//'))

for n in ${files[@]}; do make BASE_NAME=$n \
       GENOME_DIR=ref/dicty_myco_merged_star_index \
       QUANT_MODE='gene' \
       BAM_FILE=fastq/${n}.DdMm.star/Aligned.sortedByCoord.out.bam; done
