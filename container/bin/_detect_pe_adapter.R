


library(ShortRead) # BiocManager::install("ShortRead")


R1 <- FastqStreamer("data/fastq/Lung16_S11_L001_R1.fastq.gz",n = 1000000)
R2 <- FastqStreamer("data/fastq/Lung16_S11_L001_R2.fastq.gz",n = 1000000)

R1 <- yield(R1)
R2 <- yield(R2)
R1 <- sread(R1)
R2 <- sread(R2)

n <- 50
msk <- subseq(R1,1,n) == reverseComplement(subseq(R2,1,n))

table(subseq(R1[msk]))
table(subseq(R2[msk],n+1))

R1[msk]
R2[msk]




