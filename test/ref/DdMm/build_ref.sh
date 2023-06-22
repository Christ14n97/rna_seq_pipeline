


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Get Mycobacterium marinum genome from NCBI
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000018345.1
wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/345/GCF_000018345.1_ASM1834v1/GCF_000018345.1_ASM1834v1_genomic.fna.gz'
wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/018/345/GCF_000018345.1_ASM1834v1/GCF_000018345.1_ASM1834v1_genomic.gff.gz'

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Get Dicty genome from dictybase.org
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#FASTA taken from here: http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl
curl -kL -o dicty_gff3.zip 'http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=gff3&ID=dicty_gff3.zip'
unzip dicty_gff3.zip

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Merge the two genomes into DdMm.fasta DdMm.gff
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
gzip -dc GCF_000018345.1_ASM1834v1_genomic.fna.gz dicty_chromosomal.gz > DdMm.fasta
gzip -dc GCF_000018345.1_ASM1834v1_genomic.gff.gz | sed -E 's/(\tCDS\t)/\texon\t/' > DdMm.gff
cat usr/local/dicty/data/gff3/*.gff >> DdMm.gff

