ssh inf-51-2023@bioinf-serv2.cob.lu.se
# fastqc, multiqc, trimmomatic, hisat2 in QcEnv, HTSeq in htseq_env, samtools in assembly
conda create -n QcEnv
# install multiqc
conda install multiqc  # Install 
multiqc . 
conda install -c bioconda fastqc  
# check what is installed in the environment
conda list
# install trimmomatic
scp -r /home/inf-51-2023/Downloads/Trimmomatic-0.39/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin

environments=$(conda env list | awk '{print $1}' | tail -n +3)

# samtools in assembly env
# install Hisat2
export RNA_HOME=$(pwd)
# make a new directory for the RNA_seq_analysis
mkdir RNA_seq_analysis
# copy file to local server
cp /resources/binp28/Data/fungalTranscriptome.tgz . # (local server bin)
# unzip file
tar -xvzf fungalTranscriptome.tgz
# 2 Quality control
# 2.1 Genome and annotation
export reference=$RNA_HOME/01_Refs/genome.fna
export gff_file=$RNA_HOME/01_Refs/genes.gff
# Create the directory 01_Refs
mkdir 01_Refs
# Move genome.fna and genes.gff to 01_Refs
mv genome.fna 01_Refs
mv genes.gff 01_Refs
# 1. Is the gff version 2 or 3?
# version 3  (ID=gene1;Name=example_gene)
# 2. How big is the genome?
cat genome.fna | grep -v ">" | wc -c
# 45946985
3. How many scaffolds are there?
cat genes.gff | grep -v "#" | wc -l
# 131027
# 4. What proportion of the scaffolds do not have any annotated genes?
cat genes.gff | cut -f 1 | sort | grep "0" | wc -l
# 14642
# 14642/131027 how to calculate the proportion?
echo "scale=4; 14642/131027" | bc
# 0.1117
# 2.2 Reads
# Create the directory 02_Quality
mkdir 02_Quality
# Move all of the raw sample reads ending with .fastq to 02_Quality
mv *.fastq 02_Quality
# 5. How many reads are there for each replicate
cd 02_Quality
# to run all lines at once
for file in *.fastq
do
  echo $file
  cat $file | wc -l
done
# to pick only one line in each read , as each read has 4 lines
for file in *.fastq; do
    total_lines=$(wc -l < "$file")
    let "total_reads = $total_lines / 4"
    echo "$file has $total_reads reads"
done
# fh1.fastq has 300000 reads
# fh2.fastq has 100000 reads
# fh3.fastq has 100000 reads
# ref1.fastq has 300000 reads
# ref2.fastq has 100000 reads
# ref3.fastq has 100000 reads
# 6. What is the read length?
for file in *.fastq; do
    read_length=$(awk 'NR%4==2 {print length; exit}' "$file")
    echo "$file has a read length of $read_length"
done
# fh1.fastq has a read length of 50
# fh2.fastq has a read length of 50
# fh3.fastq has a read length of 50
# ref1.fastq has a read length of 50
# ref2.fastq has a read length of 50
# ref3.fastq has a read length of 50
# 2.3 Quality assessment with FastQC
# Per base sequence quality
# • Per tile sequence quality
# • Per sequence quality scores
# • Per base sequence content
# • Per sequence GC content
# • Per base N content
# • Sequence Length Distribution
# • Sequence Duplication Levels
# Overrepresented sequences
# • Adapter Content
mkdir fastQC_raw
fastqc --threads 6 fh1.fastq --outdir fastQC_raw # --threads 6 to use 6 cores, --outdir to specify the output directory
# To run the command on the rest of the files
ls *fastq | grep -v fh1.fastq | \
while read file; do fastqc --threads 6 $file --outdir fastQC_raw; \
echo $file processed;
done
# use a wildcard (*) to run fastqc
fastqc * fastq --threads 6 --outdir fastQC_raw
# run MultiQC on all of the FastQC output files 
multiqc fastQC_raw/ .
# 7. Do the reads have adapter sequences?
# yes
# 3 Trimming
# When should you trim?
# variant calling, genome annotation or transcriptome assembly, you
# can trim adapters with tools like BBduk, Trimmomatic, flexbar, TrimGalore (wrapper for CutAdapt, you
# need to know your adapter sequences). In these cases, you should also quality trim your reads. The UC Davis
# genome center recommends, “combine gentle quality trimming with a threshold quality score of Q15 with a
# read length filter retaining only reads longer than 35 bp in length” using the tools BBduk or Trimmomatic
# trimming your reads if you are getting a low mapping rate
# Trimming tools:
# • Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
# • TrimGalore, a wrapper of CutAdapt (https://github.com/FelixKrueger/TrimGalore)
# • BBduk (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
# • flexbar (https://github.com/seqan/flexbar)
# Trimmomatic: used to lightly trim and
# discard low quality fastq data as well as to remove adapters. It supports multi-threading and the use of
# compressed input files (gzip or bzip2).
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
# Trimmomatic is run using a JAR file (Java ARchive) containing java code. JAR files are executed using the
# command java -jar
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE fh1.fastq fh1.trim.fastq \
AVGQUAL:28 MINLEN:46 CROP:46 SLIDINGWINDOW:8:20
# Explanation of options:
# • SE stands for single end reads, as opposed to paired end reads
# • fh1.fastq is the input file
# • fh1.trim.fastq is the output file.
# • AVGQUAL is the average quality threshold. Sequences with a quality below this will be discarded.
# • MINLEN is the minimum length threshold. Shorter sequences will be dis- carded (that is the length after
# clipping).
# • CROP means that the output sequences are clipped at the 3’ end to a fixed length (46 in this case,
# excluding the last 4 bases). Not recommended without clear evidence that the excluded bases are
# clearly of lower quality.
# • SLIDINGWINDOW here means that a window size of 8 is used and if the average quality of this window
# is below 20, the sequences will be clipped from where the windows started
# return output
# Automatically using 1 threads
# Quality encoding detected as phred33
# Input Reads: 300000 Surviving: 291686 (97.23%) Dropped: 8314 (2.77%)
# TrimmomaticSE: Completed successfully
# run all files using a wildcard (*)
ls *fastq | \
while read file; do java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE $file ${file}.trim.fastq \
AVGQUAL:28 MINLEN:46 CROP:46 SLIDINGWINDOW:8:20;
echo $file processed;
done
# ls -lh to check the size of the files
# trimmed should be smaller than the untrimmed files
ls -l * .trim.fastq
# -rw-r--r--. 1 inf-51-2023 students 47829707 Apr  6  2014 fh1.fastq
# -rw-r--r--. 1 inf-51-2023 students 44122153 Jan 29 14:04 fh1.trim.fastq
# -rw-r--r--. 1 inf-51-2023 students 15924772 Apr  7  2014 fh2.fastq
# -rw-r--r--. 1 inf-51-2023 students 14547810 Jan 29 14:20 fh2.fastq.trim.fastq
# -rw-r--r--. 1 inf-51-2023 students 15922479 Apr  7  2014 fh3.fastq
# -rw-r--r--. 1 inf-51-2023 students 14575140 Jan 29 14:20 fh3.fastq.trim.fastq
# -rw-r--r--. 1 inf-51-2023 students  1348060 Jan 29 11:53 multiqc_report.html
# -rw-r--r--. 1 inf-51-2023 students 47832002 Apr  6  2014 ref1.fastq
# -rw-r--r--. 1 inf-51-2023 students 44245038 Jan 29 14:20 ref1.fastq.trim.fastq
# -rw-r--r--. 1 inf-51-2023 students 15917226 Apr  6  2014 ref2.fastq
# -rw-r--r--. 1 inf-51-2023 students 14760682 Jan 29 14:20 ref2.fastq.trim.fastq
# -rw-r--r--. 1 inf-51-2023 students 15911839 Apr  6  2014 ref3.fastq
# -rw-r--r--. 1 inf-51-2023 students 14740252 Jan 29 14:20 ref3.fastq.trim.fastq
# run fastqc on the trimmed files
# To run the command on the rest of the files
ls *trim.fastq | \
while read file; do fastqc --threads 6 $file --outdir fastQC_raw; \
echo $file processed;
done
# run multiqc on the trimmed files (in fastQC_raw)
multiqc ./*trim_fastqc*
# copy multiqc_report.html to local computer
scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/RNA_seq_analysis/02_Quality/multiqc_report.html .
# 8. Within each trimmed.fastq file, do the reads have the same or different lengths?
# different
for file in *fastq.trim.fastq; do
    total_lines=$(wc -l < "$file")
    let "total_reads = $total_lines / 4"
    echo "$file has $total_reads reads"
done
# fh1.fastq.trim.fastq has 291686 reads
# fh2.fastq.trim.fastq has 96357 reads
# fh3.fastq.trim.fastq has 96557 reads
# ref1.fastq.trim.fastq has 292459 reads
# ref2.fastq.trim.fastq has 97740 reads
# ref3.fastq.trim.fastq has 97642 reads
for file in *fastq.trim.fastq; do
    read_length=$(awk 'NR%4==2 {print length; exit}' "$file")
    echo "$file has a read length of $read_length"
done
# fh1.fastq.trim.fastq has a read length of 46
# fh2.fastq.trim.fastq has a read length of 46
# fh3.fastq.trim.fastq has a read length of 46
# ref1.fastq.trim.fastq has a read length of 45
# ref2.fastq.trim.fastq has a read length of 46
# ref3.fastq.trim.fastq has a read length of 46
# html file: If there are additional adapters, you can use the ILLUMINACLIP option to specify a fasta file containing
ILLUMINACLIP:~/bin/Trimmomatic -0.38/adapters/TruSeq3-SE.fa:2:10:3:1
# 4 Alignment (aka mapping)
# Combinations of tools and methods matter (Corchete et al. (2020), Liu et al. (2022))
# • Alignment and adapter trimming (optional) are less important to the results than counting and differ-
# ential expression methods (Corchete et al. (2020)) (but see the previous point)
# Reference-guided aligners: - STAR (https://github.com/alexdobin/STAR) - hisat2 (http://daehwankimlab.
# github.io/hisat2/) - StringTie (https://ccb.jhu.edu/software/stringtie
# 4.1 Genome indexing
# move to 01_Refs
cd 01_Refs
# build index
echo $reference
export reference=genome.fna
echo $reference
echo $gff_file 
export gff_file=genes.gff
echo $gff_file 
/home/inf-51-2023/bin/hisat2_2.2.1/hisat2-build $reference ${reference%.fna} 2> ${reference%.fna}_${gff_file%.gff}_hisat-build.log
ls ${reference%.fna}.* # To see the index files
# 4.2 Aligning reads to the index
# create the directory 03_Mapping
mkdir 03_Mapping
ln -s ../02_Quality/fh1.trim.fastq # ln -s to create a symbolic link
# or for all files
ls ../02_Quality/*.fastq.trim.fastq | while read line;
do ln -s $line;
done
# map the trimmed reads for the first replicate (fh1) to the genome index:
/home/inf-51-2023/bin/hisat2_2.2.1/hisat2 -p 4 --max-intronlen 5000 -U fh1.fastq.trim.fastq -x /home/inf-51-2023/RNA_seq_analysis/01_Refs/${reference%.fna} -S fh1.sam --summary-file fh1.sam.summary
# alternative ways to run the command
/home/inf-51-2023/bin/hisat2_2.2.1/hisat2 -p 4 --max-intronlen 5000 -U /home/inf-51-2023/RNA_seq_analysis/02_Quality/fh1.fastq.trim.fastq -x /home/inf-51-2023/RNA_seq_analysis/01_Refs/genome -S fh1_1.sam --summary-file fh1_1.sam.summary
# -p Number of threads
# • -max-intronlen Maximum intron size
# • -U Input file (unpaired reads)
# • -x Base name of index
# • -S Output file with aligned reads in SAM format
# • –summary-file Summary information about how mapping went
# run the rest of the files
ls *fastq.trim.fastq |  grep -v fh1.fastq | while read line;
do /home/inf-51-2023/bin/hisat2_2.2.1/hisat2 -p 4 --max-intronlen 5000 -U $line -x /home/inf-51-2023/RNA_seq_analysis/01_Refs/${reference%.fna} 
-S ${line%.fastq.trim.fastq}.sam --summary-file ${line%.fastq.trim.fastq}.sam.summary;
done
# 4.2.1 Use SAM flags to assess the alignment
ls fh1.sam*
less fh1.sam.summary
# 9. Was the mapping successful? What percentage of reads mapped successfully (overall alignment rate)
# in the fh1.sam?
# 97.60%
# SAM file starts with a series of lines beginning with @. These are header lines that contain information
# focus on the sequence lines, and count the number of sequences with different SAM flags
grep -v "^@" fh1.sam | cut -f 2 | sort -n | uniq -c 
# 143137 0 # 0: forward strand
#    7002 4 # 4: read unmapped
#  141547 16 # 16: reverse strand
#    1952 256 # 256: not primary alignment
#    1923 272 # 272: reverse strand, not primary alignment
# 10. Use the flags to count how many reads didn’t align to the reference according to the SAM file?
# 7002 + 1952 + 1923 = 10877
awk ' $2 == 4 || $2 == 256 || $2 == 272 ' fh1.sam | wc -l
# 10877
# 11. What does the 256 code for the SAM file mean?
# not primary alignment
# 12. How many reads do we have in forward and reverse direction? Don’t count reads with codes 256 and 272.
# 143137 + 141547 = 284684
awk ' $2 == 0 || $2 == 16 ' fh1.sam | wc -l
# 284684
# check remaining samples, if they are similar file size to fh1.sam
ls -l
# the answer is no
# 4.2.2 Map the remaining samples (I have done earlier)
for reads in *trim.fastq;
  do
    file=$(basename $reads .trim.fastq);
~/path/to/hisat/hisat2 -p 4 --max-intronlen 5000 -U ${file}.trim.fastq -x genome -S ${file}.sam --summary-file ${file}.sam.summary;
echo ${file} mapped;
done
# 13. Why have we used an intron value (--max-intronlen 5000) which is a small fraction of the default?
# The default is 500,000 bp, which is too large for our genome
# genome size
cat genome.fna | grep -v ">" | wc -c
# 45946985
# 4.3 Sorting and indexing mapped reads
# se samtools to convert the mapped reads to BAM format (view), sort them by read position (sort) and create an index (index) for each file
# change environment to assembly
conda activate assembly
# convert and pipe directly to sort
# Add import statement for samtools 
THREADS=6
samtools view -@ $THREADS -bS fh1.sam | samtools sort -@ $THREADS - -o fh1.sorted.bam 
# run the rest of the files
THREADS=6
ls *sam | grep -v fh1.sam | while read line;
do samtools view -@ $THREADS -bS $line | samtools sort -@ $THREADS - -o ${line%.sam}.sorted.bam;
done
# index 
samtools index -@ $THREADS fh1.sorted.bam
# run the rest of the files
ls *sorted.bam | grep -v fh1.sorted.bam | while read line;
do samtools index -@ $THREADS $line;
done
# 5 Quantification
# featureCounts is very flexible, however make mistakes with large effects on the counts
# StringTie can identify, assign reads to, and quantify transcripts, but appears to have a high error rate
# HTSeq counts one alignment per read, If there are two or more alignments for a read HTSeq takes the one
# with the highest score, primary alignment
# count with samtools
samtools view -c -F 256 fh1.sam # -c count, -F 256 not primary alignment
# 291686
# Compare you result with the HISAT2 summary file output (cat fh1.sam.summary)
cat fh1.sam.summary
# # 291686, same as above
# 5.1 Install HTSeq
# install HTSeq
conda create -n htseq_env
conda activate htseq_env
conda install -c bioconda htseq
# Create an output directory for your read count files (04_Counting) in your RNAseq analysis directory
mkdir 04_Counting
# run HTSeq on fh1.sam within the 03_Mapping directory. 
cd 03_Mapping
# Look at the options
htseq-count -h
# specify which genomic feature htseq-count should use to evaluate aligned reads.
# Features are identified in the 3rd column of the genes.gff (use head -10 genes.gff to check it out)
# genes.gff (use head -10 genes.gff)
cat genes.gff | head -10
# third column genes.gff
cat genes.gff | cut -f 3 | sort | uniq -c
# perform a differential whole gene expression analysis, reads that map to CDS features
cat genes.gff | grep "CDS" | wc -l
# option -t to specify the feature type in the htseq-count command (exon is the default)
# we ultimately want to count all of the reads mapping to each gene.
# genes are identified by name in the last column of genes.gff
# option -i to specify this identifier in the command below
# View the output before saving (run in 03_Mapping): 
htseq-count -f sam -r pos -s no -t CDS -i name -m intersection-nonempty fh1.sorted.bam /home/inf-51-2023/RNA_seq_analysis/01_Refs/genes.gff
# 100000 GFF lines processed.
# 200000 GFF lines processed.
# 209798 GFF lines processed.
# 100000 alignment records processed.
# 200000 alignment records processed.
# 295561 alignment records processed.
# View the output before saving (run in 03_Mapping)
htseq-count -f sam -r pos -s no -t CDS -i name -m intersection-nonempty fh1.sorted.bam /home/inf-51-2023/RNA_seq_analysis/01_Refs/genes.gff 
# Save results to file (run in 03_Mapping), output saved in 04_Counting
htseq-count -f sam -r pos -s no -t CDS -i name -m intersection-nonempty fh1.sam /home/inf-51-2023/RNA_seq_analysis/01_Refs/genes.gff | grep -v "^__" > ../04_Counting/fh1.count # Two underscores
# Count the reads for the remaining samples
ls *.sam | grep -v fh1.sam | while read line;
do htseq-count -f sam -r pos -s no -t CDS -i name -m intersection-nonempty $line /home/inf-51-2023/RNA_seq_analysis/01_Refs/genes.gff | grep -v "^__" > ../04_Counting/${line%.sam}.count;
done
# 100000 GFF lines processed.
# 200000 GFF lines processed.
# 209798 GFF lines processed.
# 100000 alignment records processed.
# 200000 alignment records processed.
# 295561 alignment records processed.
# 14. What did grep -v do to our final count file?
# remove lines starting with __
# • -f sam or bam format
# • -r whether the alignment is sorted by position (pos) or read name (name)
# • -s specifies whether the reads are forward stranded (yes, default), reverse stranded (reverse), or not
# stranded (no)
# • -t sub feature type (3rd column in GFF file) to be used
# • -i GFF attribute to be used as identifier of the meta-feature of interest (what we want to count.)
# • -m mode to handle reads overlapping more than one feature (here, CDS).
# Counting reads in features with htseq-count
# https://htseq.readthedocs.io/en/release_0.11.1/count.html

# 5.3 Generate a count matrix
# produce a read count matrix from these individual count files using paste
#make sure the gene IDs are all in the same order
cd ../04_Counting
ls *.count | while read file; do cat $file | cut -f 1 > $file.id; done
# Confirm the IDs are in the same order
md5sum *.id
# create the matrix file
paste *.count | cut -f 1,2,4,6,8,10,12 > count.matrix.txt
less count.matrix.txt
# Final matrix with column headers, use echo and cat
echo -e "#gene_id\tfh1\tfh2\tfh3\tref1\tref2\tref3" > header.txt #-e to enable interpretation of backslash escapes
cat header.txt count.matrix.txt | tr -d "#" > final.matrix.txt
# removed # because they are # special characters in R and are used for commenting the code
# move to RStudio
scp -r inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/RNA_seq_analysis/04_Counting .
# create a directory for your differential expression analysis results, 05_DiffExpr.
mkdir 05_DiffExpr





