


NCPU := 6
GENOME := DdMm
STAR_FLAGS := --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAA,CTGTCTCTTATACACATCT
STAR_INDEX_FLAGS :=


# Paths to the softwares
STAR := STAR
SAMTOOLS := samtools


usage:
	@echo "usage: rnaseq_pipeline <option> target"
	@echo "options:" 
	@echo "    -k               Continue on error" 
	@echo "    -j <int>         Parallel execution on <int> jobs"
	@echo "    -n               Dry run: don't run anything, just show commands"
	@echo "targets:" 
	@echo "    file.bam         Find RepA gene perform circular rotation of a FASTA"
	@echo "Example local usage"
	@echo "    docker run unigebsp/rnaseq_pipeline"
	@echo "Example usage on Baobab:"
	@echo "    # ask for a 16 cpu node for 8 hours"
	@echo "    salloc -p shared-cpu --time=8:00:00 --mem=64G -n1 --cpus-per-task=16"
	@echo "    # run the pipeline"
	@echo "    singularity run docker://unigebsp/rnaseq_pipeline"


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# STAR mapping
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

.PRECIOUS:%.$(GENOME).star/Aligned.sortedByCoord.out.bam
%.$(GENOME).star/Aligned.sortedByCoord.out.bam:%.fastq.gz
	mkdir -p '$(@D)' \
	&& $(STAR) $(STAR_FLAGS) \
	        --runThreadN $(NCPU) \
          --outSAMtype BAM SortedByCoordinate \
          --limitBAMsortRAM 4000000000 \
          --outMultimapperOrder Random \
          --outSAMmultNmax 1 \
          --outReadsUnmapped Fastx \
          --readFilesCommand gunzip -c \
          --quantMode GeneCounts \
          --genomeDir '/ref/$(GENOME).star_index/' \
          --outFileNamePrefix '$(@D)/' \
          --readFilesIn $^ \
   && pigz $(@D)/Unmapped.out.*


%.star_index/:%.fasta %.gtf
	$(STAR) --runMode genomeGenerate \
			$(STAR_INDEX_FLAGS) \
			--runThreadN $(NCPU) \
	    --genomeDir '$@/' \
	    --genomeFastaFiles '$(word 1,$^)' \
	    --sjdbGTFfile '$(word 2,$^)'

%.star_index/:%.fasta %.gff
	$(STAR) --runMode genomeGenerate \
	    $(STAR_INDEX_FLAGS) \
	    --runThreadN $(NCPU) \
	    --genomeDir '$@/' \
	    --genomeFastaFiles '$(word 1,$^)' \
	    --sjdbGTFfile '$(word 2,$^)' \
	    --sjdbGTFtagExonParentTranscript Parent


#-#-#-#-#-#-#-#-#
# SAMTOOLS rules
#-#-#-#-#-#-#-#-#
%.fasta.fai:%.fasta
	$(SAMTOOLS) faidx "$<"

%.bam.bai:%.bam
	$(SAMTOOLS) index "$<" 

%.bam.flagstat:%.bam
	$(SAMTOOLS) flagstat "$<" > $@

%.bam.idxstat:%.bam
	$(SAMTOOLS) idxstat "$<" > $@



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
	  bin/./bam_quantify \
	  --out '$(abspath $(OUTPATH))/$@' \
	  --mode $(QUANT_MODE) \
	  --gtf='ref/dicty_myco_merged_annotation.gff' \
	  $(BAM_FILE)



