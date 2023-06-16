
STAR_FLAGS += --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAA,CTGTCTCTTATACACATCT


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Define paths to the softwares
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
STAR := ml GCC/8.3.0 STAR/2.7.4a && STAR
SAMTOOLS := ml GCC/8.3.0 SAMtools/1.10 && samtools



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STAR mapping
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Paired-end version
.PRECIOUS:%.GRCm39.star/Aligned.sortedByCoord.out.bam
%.GRCm39.star/Aligned.sortedByCoord.out.bam:%.fastq.gz
	mkdir -p '$(@D)' && \
        $(STAR) $(STAR_FLAGS) \
          --runThreadN 8 \
          --outSAMtype BAM SortedByCoordinate \
          --limitBAMsortRAM 4000000000 \
          --outMultimapperOrder Random \
          --outSAMmultNmax 1 \
          --outReadsUnmapped Fastx \
          --readFilesCommand gunzip -c \
          --quantMode GeneCounts \
          --genomeDir './data/ref/GRCm39_vM27/star_index' \
          --outFileNamePrefix '$(@D)/' \
          --readFilesIn $^




#-#-#-#-#-#-#-#-#
# SAMTOOLS 
#-#-#-#-#-#-#-#-#
%.bam.bai:%.bam
	$(SAMTOOLS) index $<

%.bam.flagstat:%.bam
	$(SAMTOOLS) flagstat $< > $@

%.bam.idxstat:%.bam
	$(SAMTOOLS) idxstat $< > $@




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Compute genome coverage
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
%.fwd.wig:%.bam
	./src/bin/bam_coverage --strand=forward --out '$@' $<
%.rev.wig:%.bam
	./src/bin/bam_coverage --strand=reverse --out '$@' $<
%.wig:%.bam
	./src/bin/bam_coverage --strand=both --out '$@' $<



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Quantify Reads into each gene
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
data/fastq/quantif_GRCm39.%.rds:
	./local/bin/bam_quantify \
	  --out '$@' \
	  --mode '$*' \
	  --gtf='data/ref/GRCm39_vM27/gencode.vM27.annotation.gtf' \
	  data/fastq/*.GRCm39.star/Aligned.sortedByCoord.out.bam

data/fastq/quantif_spneum.%.rds:
	./local/bin/bam_quantify \
	  --out '$@' \
	  --mode '$*' \
	  --gtf='data/ref/atcc6303_s_pneumoniae/Streptococcus_pneumoniae_ATCC_6303.gff' \
	  data/fastq/*.spneum.star/Aligned.sortedByCoord.out.bam

data/fastq/quantif_lm.%.rds:
	./local/bin/bam_quantify \
	  --out '$@' \
	  --mode '$*' \
	  --gtf='data/ref/l_murinis/GCF_003288115.1_ASM328811v1_genomic.gtf' \
	  data/fastq/*.lm.star/Aligned.sortedByCoord.out.bam

data/fastq/quantif_influenza.%.rds:
	./local/bin/bam_quantify \
	  --out '$@' \
	  --mode '$*' \
	  --gtf='data/ref/influenza_A_virus/influenza_A_virus_segments.gff3' \
	  data/fastq/*.influenza.star/Aligned.sortedByCoord.out.bam










