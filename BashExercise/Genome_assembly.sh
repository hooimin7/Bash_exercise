# velvet (~/miniconda3/envs/assembly)
sra-toolkit (~/miniconda3/envs/assembly)
spades (~/miniconda3/envs/assembly)
samtools (~/miniconda3/envs/assembly)
perl (~/miniconda3/envs/assembly)
ncbi-blast (~/miniconda3/envs/assembly)
bbmap (~/miniconda3/envs/assembly)
htslib (~/miniconda3/envs/assembly)
ragtag (~/miniconda3/envs/assembly)
# Install miniconda for bash
# https://docs.conda.io/projects/miniconda/en/latest/
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
# Genome assembly
conda -V
conda info
conda update -- help
conda update conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install mamba
conda install conda-forge::mamba
conda env list
# Create a new Conda environment
conda create -n assembly
# environment location: /home/inf-51-2023/anaconda3/envs/assembly
# activate 
conda activate assembly
conda env list
# base                     /home/inf-51-2023/anaconda3
# assembly              *  /home/inf-51-2023/anaconda3/envs/assembly
conda install velvet
https://www.overleaf.com/project/5ffed362d7828d64b0f31630
cd
velvetg
velveth
########### spades
conda install spades
spades.py --test
tree
# Install the SRA Toolkit 
conda install sra-tools
# check if sra tools installed
fastq-dump
# data download
cd
mkdir BacterialGenomeAssembly
cd BacterialGenomeAssembly
# NCBI’s Short Read Archive (SRA)
# http://www.ncbi.nlm.nih.gov/sra/SRX377522
# 1. What is the bacterium called and from where was it isolated?
# Geobacillus thermopakistaniensis MAS1, isolated from Pakistan, isolation source	hot springs 
# https://www.ncbi.nlm.nih.gov/biosample/SAMN02371536
# 2. The strategy used is WGS. What does shotgun mean (S in WGS)?
# whole genome sequencing sequences the genome in a highly ordered manner, while shotgun sequencing sequences the genome in a more random and fragmented way before assembly. 
# Shotgun sequencing is a method of sequencing used for sequence large genome targets like Whole Genome
# Shotgun sequencing is a laboratory technique for determining the DNA sequence of an organism's genome. The method involves randomly breaking up the genome into small DNA fragments that are sequenced individually.
# 3. What sequencing technology has been used?
# https://www.ncbi.nlm.nih.gov/genbank/wgs/
# 4. What kind of libraries (layout) have been used?
# Name: MAS1
# Instrument: Illumina HiSeq 2500
# Strategy: WGS
# Source: GENOMIC
# Selection: unspecified
# Layout: PAIRED
# 5. How long are the reads (you can double check this later)?
# length = 2*150 basepairs, https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR1029680&display=metadata
# 6. How many reads are there? One spot represents both the forward and
# the corresponding reverse read.
echo $((325997 * 2))
# 651994
# 7. How many bases have been sequenced (take the results from the two
# questions above and multiply them)? Compare with what is stated
# on the page. 
echo $((300*325997)) 
# 97799100 bases
# 325,997 (spots),  97.8 M(bases), 38.7Mb (file size)
# 8. This bacterium has a genome size of about 3.5 Mb. What is the sequencing
# depth (in average)
echo $((97800/3500))
# 27
# 9. how many sequences there are
# FASTQ file consists of four lines (header, sequence, separator, and quality scores), 
# dividing the total number of lines by 4 gives the number of sequences.
echo $(cat SRR1029680.fastq | wc -l)/4 | bc -l 
# 651994.00000000000000000000
# 10. check how many base pairs there are (calculate one by one base pair, total length)
awk 'NR%4==2 {sum += length($0)} END {print sum}' SRR1029680.fastq
# 97799100
# 11. check how long a read is, only difference between the two lines is that a ’1’ is substituted
# to a ’2’. ’1’ denotes forward read and ’2’ reverse read
# FASTQ format
# https://en.wikipedia.org/wiki/FASTQ_format
head -8 SRR1029680.fastq | grep @SRR
# 150
fastq-dump --origfmt --split-files SRR1029680
# Compare the new files with the old one
# check 1. Why is the size of the original file larger than the sum of the two
# new ones? Remember that normally one character corresponds to one byte.
# 1. The headers are different with more information (e.g. read length)included in
# the manually downloaded file.
# check 2. What is the reason that there is one line more in the original file
# as compared to the two new ones? (Don’t spend time on this one).
# 2. The last read of first file needs an extra newline before the first read of the
# second file is added.
ls -lh # Note that M means megabytes and G means gigabytes .
wc -l *.fastq
# 1303988 SRR1029680_1.fastq
# 1303988 SRR1029680_2.fastq
# 2607976 total
head -2 *.fastq
# run velveth (velveth require a lot of memory (RAM).)
# velveth which prepares a data set for assembly
# velvetg which makes the actual assembly
# Required RAM = (-109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads
# - 51092*K)/1024
# Read size is in bases.
# • Genome size is in millions of bases (Mb)
# • Number of reads is in millions
# • K is the kmer hash value used in velveth (default 21, but use 31)
# • The answer is in MB.
# 3. What is the memory requirement?
# 4. Optional: Imagine that you want to assemble a bird genome. The
# genome size is 1 Gb. You will sequence the genome to 50 times coverage
# using 2x150 paired end reads. How many reads are required and what
# will the memory requirement be? 8.5 GB
# 5. How much memory do you have on your computer?
# 6 Why does velvet have two steps, velvetg and velveth, instead of one? velveth prepares a data set for assembly, velvetg makes the actual assembly
# Setting up of directories and indexing is only needed once. The structure and
# index can be reused many times (for the same) if needed.
# 7.Compare the results with the assembly from other assemblers once
# you have generated them. The velvetg assembly output file is contigs.fa.
# For example, to find the longest contig.
# 7. Short contigs are generated (largest is 3708 basepairs only with
# kmer 31). A larger kmer size such as 51 gives larger contigs (largest= 15733).

velveth Velvet 31 -fastq -shortPaired SRR1029680_1.fastq SRR1029680_2.fastq
# To understand the options :
velveth -h | less
# assemble the reads
velvetg Velvet
# find the longest contig (output file is contigs.fa.)
cat contigs.fa | grep \> | cut -d _ -f4 | sort -n
# 8. Improve the velvetg assembly by optimizing some of the
# options. You may take a look at for some tips on which settings to change.
# 8. Running velveth with a range instead of one randomly picked kmer will
# improve the results, but take more time if the number of kmers tested is large
# check *** http://evomics.org/learning/assembly-and-alignment/velvet/
# Assembly statistics
conda install bbmap
stats.sh Velvet/contigs.fa 
A       C       G       T       N       IUPAC   Other   GC      GC_stdev                                                
0.2414  0.2599  0.2594  0.2393  0.0000  0.0000  0.0000  0.5193  0.0853                                                  
                                                                                                                        
Main genome scaffold total:             18135                                                                           
Main genome contig total:               18135                                                                           
Main genome scaffold sequence total:    4.384 MB                                                                        
Main genome contig sequence total:      4.384 MB        0.000% gap                                                      
Main genome scaffold N/L50:             2355/499                                                                        
Main genome contig N/L50:               2355/499                                                                        
Main genome scaffold N/L90:             12305/93                                                                        
Main genome contig N/L90:               12305/93                                                                        
Max scaffold length:                    4.306 KB                                                                        
Max contig length:                      4.306 KB                                                                        
Number of scaffolds > 50 KB:            0                                                                               
% main genome in scaffolds > 50 KB:     0.00%


Minimum         Number          Number          Total           Total           Scaffold
Scaffold        of              of              Scaffold        Contig          Contig  
Length          Scaffolds       Contigs         Length          Length          Coverage
--------        --------------  --------------  --------------  --------------  --------
    All                 18,135          18,135       4,384,088       4,384,088   100.00%
     50                 18,135          18,135       4,384,088       4,384,088   100.00%
    100                 10,968          10,968       3,837,565       3,837,565   100.00%
    250                  4,454           4,454       2,940,851       2,940,851   100.00%
    500                  2,345           2,345       2,190,318       2,190,318   100.00%
   1 KB                    726             726       1,052,825       1,052,825   100.00%
 2.5 KB                     28              28          81,989          81,989   100.00%
# Improve the velvetg assembly
# https://evomics.org/learning/assembly-and-alignment/velvet/

# how to check how many cores bash
# https://www.cyberciti.biz/faq/check-how-many-cpus-are-there-in-linux-system/
# Spades alters the k-mer size depending on the coverage of the genome.
# For low coverage regions k-mer size should be low and for high coverage
# regions high Assembly with spades.
# The values were chosen by using the help on:
# check https://cab.spbu.ru/files/release3.15.3/manual.html#sec3.4
spades.py --pe1 -1 SRR1029680_1.fastq --pe1 -2 SRR1029680_2.fastq -k 21,33,55,77 --careful -t 16 -o Spades
# Statistics Guide
# https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/statistics-guide/
# to check the statistic 
# 9. How big is the genome? Use the scaffolds.fasta file.
# 9. The genome size is 3520016 ~3.5 Mbps
stats.sh Spades/scaffolds.fasta
A       C       G       T       N       IUPAC   Other   GC      GC_stdev
0.2387  0.2565  0.2652  0.2396  0.0001  0.0000  0.0000  0.5217  0.0934

Main genome scaffold total:             336
Main genome contig total:               339
Main genome scaffold sequence total:    3.520 MB
Main genome contig sequence total:      3.520 MB        0.008% gap
Main genome scaffold N/L50:             6/206.115 KB
Main genome contig N/L50:               7/171.476 KB
Main genome scaffold N/L90:             19/57.572 KB
Main genome contig N/L90:               21/52.899 KB
Max scaffold length:                    447.446 KB
Max contig length:                      390.948 KB
Number of scaffolds > 50 KB:            21
% main genome in scaffolds > 50 KB:     93.40%


Minimum         Number          Number          Total           Total           Scaffold
Scaffold        of              of              Scaffold        Contig          Contig  
Length          Scaffolds       Contigs         Length          Length          Coverage
--------        --------------  --------------  --------------  --------------  --------
    All                    336             339       3,520,252       3,519,965    99.99%
     50                    336             339       3,520,252       3,519,965    99.99%
    100                    311             314       3,518,078       3,517,791    99.99%
    250                    153             156       3,487,870       3,487,583    99.99%
    500                     58              61       3,460,672       3,460,385    99.99%
   1 KB                     44              47       3,450,932       3,450,645    99.99%
 2.5 KB                     36              39       3,438,051       3,437,764    99.99%
   5 KB                     31              34       3,422,386       3,422,099    99.99%
  10 KB                     28              31       3,401,823       3,401,536    99.99%
  25 KB                     22              25       3,320,520       3,320,233    99.99%
  50 KB                     21              24       3,287,889       3,287,602    99.99%
 100 KB                     11              14       2,648,875       2,648,588    99.99%
 250 KB                      4               6       1,465,100       1,464,913    99.99%
 how to copy file from local to server
scp -r /Users/med-snt/Downloads/ncbi_dataset/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/
# 10. How big is the genome?
# 10. Genome is 3,592,666 basepairs long. Note that this is for Geobacillus kaustophilus HTA426
# Installation Minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.26_x64-linux/minimap2
# Installation Ragtag
# install with conda
conda install -c bioconda ragtag
# find my file
ls ncbi_dataset/ncbi_dataset/data/GCF_000009785.1/
# correct a query assembly
ragtag.py correct ref.fasta query.fasta
ragtag.py correct ncbi_dataset/ncbi_dataset/data/GCF_000009785.1/GCF_000009785.1_ASM978v1_genomic.fna Spades/scaffolds.fasta
# scaffold a query assembly
ragtag.py scaffold ref.fasta query.fasta
ragtag.py scaffold ncbi_dataset/ncbi_dataset/data/GCF_000009785.1/GCF_000009785.1_ASM978v1_genomic.fna Spades/scaffolds.fasta 
# scaffold with multiple references/maps
ragtag.py scaffold -o out_1 ref1.fasta query.fasta
ragtag.py scaffold -o out_2 ref2.fasta query.fasta
ragtag.py merge query.fasta out_*/*.agp other.map.agp

# use Hi-C to resolve conflicts
ragtag.py merge -b hic.bam query.fasta out_*/*.agp other.map.agp

# make joins and fill gaps in target.fa using sequences from query.fa
ragtag.py patch target.fa query.fa
# check ragtag output
cat ragtag_output/ragtag.scaffold.stats 
# check the statistic ragtag at right environment (assembly)
stats.sh ragtag_output/ragtag.scaffold.fasta

