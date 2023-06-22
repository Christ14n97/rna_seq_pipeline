
# Build the container
docker build --platform linux/amd64 -t rnaseq_pipeline container/


# Index a reference genome
docker run --rm \
  -v "$(pwd):/data" \
  rnaseq_pipeline \
    STAR_INDEX_FLAGS="--genomeSAindexNbases 11" \
    test/ref/DdMm/DdMm.star_index/

# Generate BAM files
docker run --rm \
  -v "$(pwd):/data" \
  -v "$(pwd)/test/ref/DdMm/:/ref" \
  rnaseq_pipeline \
    GENOME=DdMm \
    test/fastq/sample1_100k.DdMm.star/Aligned.sortedByCoord.out.bam



