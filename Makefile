
STAR_FLAGS += --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAA,CTGTCTCTTATACACATCT


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Define paths to the softwares
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
STAR := STAR
SAMTOOLS := samtools



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STAR mapping
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Paired-end version
.PRECIOUS:%.DdMm.star/Aligned.sortedByCoord.out.bam
%.DdMm.star/Aligned.sortedByCoord.out.bam:%.fastq.gz
	mkdir -p 'test/$(@D)' && \
        $(STAR) $(STAR_FLAGS) \
          --runThreadN 8 \
          --outSAMtype BAM SortedByCoordinate \
          --limitBAMsortRAM 4000000000 \
          --outMultimapperOrder Random \
          --outSAMmultNmax 1 \
          --outReadsUnmapped Fastx \
	  --readFilesIn 'test/fastq/' \
          --readFilesCommand gunzip -c \
          --quantMode GeneCounts \
          --genomeDir 'test/reference_genomes/dicty_myco_merged_star_index' \
          --outFileNamePrefix 'test/$(@D)/' \
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
	./bin/bam_coverage --strand=forward --out '$@' $<
%.rev.wig:%.bam
	./bin/bam_coverage --strand=reverse --out '$@' $<
%.wig:%.bam
	./bin/bam_coverage --strand=both --out '$@' $<



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Quantify Reads into each gene
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
data/fastq/quantif_GRCm39.%.rds:
	./bin/bam_quantify \
	  --out '$@' \
	  --mode '$*' \
	  --gtf='data/ref/GRCm39_vM27/gencode.vM27.annotation.gtf' \
	  data/fastq/*.GRCm39.star/Aligned.sortedByCoord.out.bam




