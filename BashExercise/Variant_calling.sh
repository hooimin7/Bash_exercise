ssh-copy-id -p 2222 -i ~/.ssh/id_ed25519.pub inf-51-2023@bioinf-serv2.cob.lu.se
ssh inf-51-2023@bioinf-serv2.cob.lu.se
scp -r /Users/med-snt/Downloads/bowtie2-2.5.3-sra-linux-x86_64/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin
mkdir Variant_calling 
mv bowtie2-2.5.3-sra-linux-x86_64/bowtie2* ./ # in bin folder (should install in bashrc, check in the bash lecture)
bowtie2-build -h
cp /resources/binp28/Data/yeastGenome.fa.gz .
cp /resources/binp28/Data/yeastReads.fastq.gz .
gunzip yeastGenome.fa.gz
gunzip yeastReads.fastq.gz
genome size yeastGenome.fa
# 12157105
program bow1e2-build for indexing
bowtie2-build yeastGenome.fa yeastGenome
# genome file should have an extension of .fa and that the index name (yeastGenome) 
# for convenience should have the same name as the genome file (without the .fa) extension
# Mapping
# he read file and the genome index files
bowtie2 -x yeastGenome -U yeastReads.fastq -S readsAligned.sam
# 100000 reads; of these:
#   100000 (100.00%) were unpaired; of these:
#     8420 (8.42%) aligned 0 times
#     84908 (84.91%) aligned exactly 1 time
#     6672 (6.67%) aligned >1 times
# 91.58% overall alignment rate
# save what is output to the screen 
bowtie2 -x yeastGenome \
-U yeastReads.fastq \
-S readsAligned.sam \
2> bow1e2.log
less bow1e2.log
# The SAM format
# Take a look at the output file:
less readsAligned.sam
# The first line tells us that the output is unsorted. There are two possible ways to sort
# the data (more on this later):
# • On genomic posi1ons.
# • On read names
# Then comes the name of the sequences in the genome file. All lines star1ng with an
# at sign (@) contain meta data
cat readsAligned.sam \
| grep -v "^@" \
| head -1 \
| tr "\t" "\n" \
| cat -n
# explain the output: https://broadinstitute.github.io/picard/explain-flags.html
# 1  read1
#      2  16
#      3  chrI
#      4  134097
#      5  42
#      6  50M
#      7  *
#      8  0
#      9  0
#     10  TGCTTTAACGCCATTTTGATTATACACATTGTATTACTTATTTTTTAACC
#     11  22222222222222222222222222222222222222222222222222
#     12  AS:i:-3
#     13  XN:i:0
#     14  XM:i:1
#     15  XO:i:0
#     16  XG:i:0
#     17  NM:i:1
#     18  MD:Z:15A34
#     19  YT:Z:UU
# Check that point 17 is correct if you know awk
cat yeastGenome.fa \ |awk ' /^>chrII/ {exit} NR==1 {next} {seq=seq $0} END {print seq} ' \
| cut -b 134097-`echo 134097+49 | bc`
cat readsAligned.sam | grep -v "^@" | head -1 | cut -f 10
# Download htslib and samtools to bin 
scp -r /Users/med-snt/Downloads/samtools-1.19.1/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin
# install htslib first
tar -xvjf htslib-1.19.1.tar.bz2
cd htslib-1.19.1
# install C or C++
./configure # Not always required
make
make install
# install on our own directory
./configure
make
make prefix=. install
# same procedure for samtools and include it in your PATH
tar -xvjf samtools-1.19.1.tar.bz2
cd samtools-1.19.1
./configure
make
make prefix=. install
# in our VariantCalling directory
samtools view -b -S -o readsAligned.bam readsAligned.sam # -b output in BAM format, -S input in SAM format
# Compare the file sizes of the original read file, the SAM file and the BAM file
ls -lh yeastReads.fastq readsAligned.sam readsAligned.bam
# -rw-r--r--. 1 inf-51-2023 students 3.2M Jan 24 13:12 readsAligned.bam
# -rw-r--r--. 1 inf-51-2023 students  19M Jan 24 11:25 readsAligned.sam
# -rw-r--r--. 1 inf-51-2023 students  11M Jan 24 10:58 yeastReads.fastq
# convert back to SAM format using a BAM file
samtools view readsAligned.bam | less
# view the data is to use the -f option. output data for sequences that mapped on the reverse strand
samtools view -f 16 readsAligned.bam > reverse.sam
# reverse.sam
# count the number, add a -c
samtools view -c -f 16 readsAligned.bam
# 45770
wc -l reverse.sam
# 45770 reverse.sam
# You can also filter out sequences that don’t meet the criteria by using -F. Let’s say that
# you want a SAM file with only those reads that mapped
samtools view -F 4 readsAligned.bam > mappedReads.sam
samtools view -q 40 -c readsAligned.bam # Here we count these reads
# 84081 (number include forward and reverse, but with good quality)
# Default sort on positions:
samtools sort -o readsAligned.sortedPosi1ons.bam readsAligned.bam # Sort on names
samtools sort -n -o readsAligned.sortedNames.bam readsAligned.bam # -n sort on names
# Compare the different sorts:
samtools view readsAligned.sortedNames.bam | head -5
samtools view readsAligned.sortedPosi1ons.bam | head -5
# samtools can take input from STDIN instead from a specified file
cat readsAligned.bam | samtools sort | samtools view
cp /Users/med-snt/Downloads/bcftools-1.19 . 
./configure
make
make prefix=. install
# needs to copy all files in 
bcftools mpileup -Ob -f yeastGenome.fa \ # -Ob output in BCF format, -f reference genome
readsAligned.sortedPosi1ons.bam \
> variants.bcf
bcftools call -m -v variants.bcf > variants.vcf # -m multiallelic calling, -v output variant sites only # -m multiallelic calling, -v output variant sites only
# understand the VCF format:
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
less variants.vcf
# For how many sites does the sequenced strain differ from the original strain reference genome
grep -v "^#" variants.vcf | wc -l
# 466
# Examining VCFs
# The number of variable sites
grep -v "^#" variants.vcf | wc -l
# The number of chromosomes these variable sites are distributed over
grep -v "^#" variants.vcf | cut -f 1 | uniq | wc -l
# Whether data is phased or not
grep -v "^#" variants.vcf | cut -f 10 | head -1
# How many reads there are per position (coverage)
grep -v "^#" variants.vcf | cut -f 8 | head -1
cp /resources/binp28/Data/birds.Z_and_5.vcf.gz .
gunzip birds.Z_and_5.vcf.gz
# look at the vcf file
grep -v ^\#  birds.Z_and_5.vcf | cut -f 1 | uniq | head
grep -v ^\# birds.Z_and_5.vcf | grep "chrZ" | wc -l
# 2012086
# A) Find out how many SNPs there are in the file through combining grep and wc to
count the relevant lines.
grep -v ^\# birds.Z_and_5.vcf | wc -l
# 4960898
# Across how many chromosomes are these SNPs distributed? Use awk t
grep -v ^\# birds.Z_and_5.vcf | cut -f 1 | uniq | wc -l
# 2
# Download vcftools
scp -r /Users/med-snt/Downloads/vcftools_0.1.13/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin
# --site-mean-depth, Calculate the mean depth per site, producing a file with two columns: the first is the site and the second is the mean depth.
vcftools --vcf birds.Z_and_5.vcf --site-mean-depth
#  This output file has the suffix ".ldepth".
# more info: https://vcftools.sourceforge.net/man_latest.html
# use awk to estimate the mean coverage across all sites and individuals
awk '{sum+=$3} END {print sum/NR}' out.ldepth.mean
# 8.1535
# how many column in the file
awk '{print NF}' out.ldepth.mean
# how many individuals the file contains
grep -v ^\# birds.Z_and_5.vcf | cut -f 10 | uniq | wc -l
# 2935612
# how many SNPs and chromosomes
grep -v ^\# birds.Z_and_5.vcf | wc -l
# 4960898
# remove the header
grep -v ^\CHROM out.ldepth.mean | head
grep -v ^\CHROM out.ldepth.mean | cut -f 1 | uniq | wc -l
# d) what the mean coverage per site is for chromosome Z
grep -v ^\CHROM out.ldepth.mean | grep "chrZ" | awk '{sum+=$3} END {print sum/NR}'
# 8.82825
# d) what the mean coverage per site is for chromosome 5
grep -v ^\CHROM out.ldepth.mean | grep "chr5" | awk '{sum+=$3} END {print sum/NR}'
# 7.6931