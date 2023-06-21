

all: $(BASE_NAME).DdMm.star/Aligned.sortedByCoord.out.bam $(BASE_NAME).bam.bai \
	$(BASE_NAME).bam.flagstat $(BASE_NAME).bam.idxstat quantif_$(BASE_NAME).$(QUANT_MODE).rds

.PHONY: all

STAR_FLAGS += --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAA,CTGTCTCTTATACACATCT

BASE_NAME := $(BASE_NAME)
vpath %.fastq.gz fastq/
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Define paths to the softwares
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
STAR := STAR
SAMTOOLS := samtools



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STAR mapping
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

GENOME_DIR := $(GENOME_DIR)

# Paired-end version
.PRECIOUS:%.DdMm.star/Aligned.sortedByCoord.out.bam
$(BASE_NAME).DdMm.star/Aligned.sortedByCoord.out.bam:$(BASE_NAME).fastq.gz
	mkdir -p 'fastq/$(@D)' && \
        $(STAR) $(STAR_FLAGS) \
          --runThreadN 8 \
          --outSAMtype BAM SortedByCoordinate \
          --limitBAMsortRAM 4000000000 \
          --outMultimapperOrder Random \
          --outSAMmultNmax 1 \
          --outReadsUnmapped Fastx \
          --readFilesCommand gunzip -c \
          --quantMode GeneCounts \
          --genomeDir $(GENOME_DIR) \
          --outFileNamePrefix 'fastq/$(@D)/' \
          --readFilesIn $^




#-#-#-#-#-#-#-#-#
# SAMTOOLS 
#-#-#-#-#-#-#-#-#

# storage
DIR := $(dir $(BAM_FILE))
# expansion, otherwise make is not getting the value
OUTPATH := $(DIR)

$(BASE_NAME).bam.bai:$(BAM_FILE)
	$(SAMTOOLS) index $< 

$(BASE_NAME).bam.flagstat:$(BAM_FILE)
	$(SAMTOOLS) flagstat "$(abspath $<)" > "$(abspath $(OUTPATH))/$@"

$(BASE_NAME).bam.idxstat:$(BAM_FILE)
	$(SAMTOOLS) idxstat "$(abspath $<)" > "$(abspath $(OUTPATH))/$@"




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Compute genome coverage
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#%.fwd.wig:%.bam
#	.././bin/bam_coverage --strand=forward --out '$@' $<
#%.rev.wig:%.bam
#	.././bin/bam_coverage --strand=reverse --out '$@' $<
#%.wig:%.bam
#	.././bin/bam_coverage --strand=both --out '$@' $<



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Quantify Reads into each gene
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
QUANT_MODE := $(QUANT_MODE)
quantif_$(BASE_NAME).$(QUANT_MODE).rds:$(BAM_FILE)
	  ../pipeline/bin/./bam_quantify \
	  --out '$(abspath $(OUTPATH))/$@' \
	  --mode $(QUANT_MODE) \
	  --gtf='reference_genomes/dicty_myco_merged_annotation.gff' \
	  $(BAM_FILE)




