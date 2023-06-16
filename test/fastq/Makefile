
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
          --readFilesCommand gunzip -c \
          --quantMode GeneCounts \
          --genomeDir '../reference_genomes/dicty_myco_merged_star_index' \
          --outFileNamePrefix '../mapped/$(@D)/' \
          --readFilesIn $^




#-#-#-#-#-#-#-#-#
# SAMTOOLS 
#-#-#-#-#-#-#-#-#

# storage
DIR := $(dir $(BAM_FILE)))
# expansion, otherwise make is not getting the value
OUTPATH := $(DIR)

%.bam.bai:$(BAM_FILE)
	$(SAMTOOLS) index $<

%.bam.flagstat:$(BAM_FILE)
	$(SAMTOOLS) flagstat $(abspath $<) > $(abspath $(OUTPATH)/$@

%.bam.idxstat:$(BAM_FILE)
	$(SAMTOOLS) idxstat "$(abspath $<)" > "$(abspath $(OUTPATH))/$@"




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Compute genome coverage
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
%.fwd.wig:%.bam
	.././bin/bam_coverage --strand=forward --out '$@' $<
%.rev.wig:%.bam
	.././bin/bam_coverage --strand=reverse --out '$@' $<
%.wig:%.bam
	.././bin/bam_coverage --strand=both --out '$@' $<



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Quantify Reads into each gene
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
quantif_GRCm39.%.rds:
	.././bin/bam_quantify \
	  --out '$@' \
	  --mode '$*' \
	  --gtf='test/reference_genomes/dicty_myco_merged_annotation.gtf' \
	  *.GRCm39.star/Aligned.sortedByCoord.out.bam




