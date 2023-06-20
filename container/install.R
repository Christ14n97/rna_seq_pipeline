install.packages(c('optparse', 'BiocManager'))

BiocManager::install(version='3.17')
BiocManager::install(c('BiocParallel', 'GenomicAlignments', 'GenomicFeatures', 'rtracklayer'))
