


library(BiocParallel)
register(MulticoreParam(workers=8))
library(GenomicAlignments)
library(Rsamtools)
library(SummarizedExperiment)

# Read given STAR output file
read_final_out <- function(f) {
  x <- read.table(f,header = FALSE, fill = TRUE, sep="\t", row.names=1,stringsAsFactors=FALSE)
  rownames(x) <- sub(" *[|] *","",sub("^ *","",rownames(x)))
  x <- setNames(x$V2,rownames(x))
  c(Tot_reads = as.integer(x["Number of input reads"]),
    Mapped_reads = as.integer(x["Uniquely mapped reads number"]),
    Mapped_reads_including_multimappers = sum(as.integer(x[c("Uniquely mapped reads number","Number of reads mapped to multiple loci","Number of reads mapped to too many loci")]))
  )
}

