
# study of Benjamino et al. (2016) : xplores the bacterial
# community in the hindgut of the termite Reticulitermes flavipes. The microbes in the
# hindgut help the termite digest the lignocellulose from the “wood food.
conda create -n amplicon_env
# install the tools in the environment
# install pandaseq using bioconda in the amplicon_env (merged paired-end reads)
# test pandaseq
pandaseq -h
conda install bioconda::pandaseq
# test pandaseq
pandaseq -h
# install Vsearch
conda install bioconda::vsearch
# https://github.com/torognes/vsearch
# make directory for the project
mkdir -p Amplicon
cp /resources/binp28/Data/AmpliconData.tar.gz .
tar -xvzf AmpliconData.tar.gz
# unzip the files fastq_files.zip
unzip fastq_files.zip
# unzip the files XML.tgz
tar -xvzf XML.tgz
# 4. How much space do the files take after decompression (use disk usage du)
du
# Answer: 3920316 KB
# This number represents the disk space usage in kilobytes. To convert it to megabytes, you can divide by 1024 (since 1 MB = 1024 KB). 
# So, in megabytes, it would be approximately 3824 MB (3920316 KB / 1024)
scp -r /Users/med-snt/Downloads/Colonies.tar.gz/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon
tar -xvzf Colonies.tar.gz
# 7. How many sequences do you get in total? 
# Answer: Total reads = 7850134
# script to count the number of reads

files=("CT_A_1.fastq" "CT_A_2.fastq" "CT_B_1.fastq" "CT_B_2.fastq" "CT_C_1.fastq" "CT_C_2.fastq" "CT_D_1.fastq" "CT_D_2.fastq" "MA_A_1.fastq" "MA_A_2.fastq" "MA_B_1.fastq" "MA_B_2.fastq" "MA_C_1.fastq" "MA_C_2.fastq")
total_reads=0

for file in "${files[@]}"
do
  lines=$(wc -l < "$file")
  reads=$((lines / 4))
  total_reads=$((total_reads + reads))
done

echo "Total reads: $total_reads"

# 8. What is the average read length?
awk 'NR==2 {print length; exit}' MA_C_1.fastq # single read counted
awk '{if(NR%4==2) {count++; sum+=length($0)}} END {print sum/count}' CT_A_2.fastq # average read length for all reads in a file
# script to count the average read length for all reads in all files
files=("CT_A_1.fastq" "CT_A_2.fastq" "CT_B_1.fastq" "CT_B_2.fastq" "CT_C_1.fastq" "CT_C_2.fastq" "CT_D_1.fastq" "CT_D_2.fastq" "MA_A_1.fastq" "MA_A_2.fastq" "MA_B_1.fastq" "MA_B_2.fastq" "MA_C_1.fastq" "MA_C_2.fastq")

for file in "${files[@]}"
do
  echo "Processing $file"
  awk '{if(NR%4==2) {count++; sum+=length($0)}} END {print sum/count}' "$file"
done
ls -l, wc -l and md5sum # to check the files
# md5sum to check the files
files=("CT_A_1.fastq" "CT_A_2.fastq" "CT_B_1.fastq" "CT_B_2.fastq" "CT_C_1.fastq" "CT_C_2.fastq" "CT_D_1.fastq" "CT_D_2.fastq" "MA_A_1.fastq" "MA_A_2.fastq" "MA_B_1.fastq" "MA_B_2.fastq" "MA_C_1.fastq" "MA_C_2.fastq")

for file in "${files[@]}"
do
  md5sum "$file"
done
# Create a folder called 1_fastqc and use a bash while loop to generate FastQC reports for each sample (in Colonies folder)
mkdir 1_fastqc
# script to generate FastQC reports for each sample
ls *fastq | \
while read file; do fastqc --threads 6 $file --outdir 1_fastqc; \
echo $file processed;
done
scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon/Colonies/1_fastqc/CT_A_1_fastqc.html .
scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon/Colonies/1_fastqc/CT_A_2_fastqc.html .
# make directory 2_trimming
mkdir 2_trimming
# trimmomatic 
java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -trimlog CT_A.log -baseout CT_A.fastq /home/inf-51-2023/Amplicon/Colonies/CT_A_1.fastq \
/home/inf-51-2023/Amplicon/Colonies/CT_A_2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:8:15 MINLEN:140 AVGQUAL:20
# take a look at the file
ls -l
wc -l *fastq
# looping through all the files and trimming them
ls ../Colonies/*_1.fastq | sed s/_1.fastq// | while read line; do id=$(echo $line | sed "s/.*\///"); java -jar ~/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -trimlog $id.log -baseout $id.fastq ${line}_1.fastq ${line}_2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:8:15 MINLEN:140 AVGQUAL:20; done

# 14. What percentage of the reads was discarded?
total_reads=$(awk 'END {print NR/4}' *P.fastq *U.fastq) | echo $total_reads # count the number of reads in the files and assign the result to total_reads variable
surviving_reads=$(awk 'END {print NR/4}' *P.fastq) echo $surviving_reads  # count the number of reads in the files and assign the result to surviving_reads variable
discarded_percentage=$(echo "scale=5; (1 - $surviving_reads / $total_reads) * 100" | bc) # calculate the percentage of reads discarded
echo "Discarded percentage: $discarded_percentage" # print the result
# 1.05% of the reads was discarded

# Make directory 3_mergereads
mkdir 3_mergereads
# script to generate FastQC reports for each sample
ls *fastq | \
while read file; do fastqc --threads 6 $file --outdir 3_fastqc_trimmed; \
echo $file processed;
done    
#
scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon/2_trimming/3_fastqc_trimmed/CT_A_1P_fastqc.html .
scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon/2_trimming/3_fastqc_trimmed/CT_A_2P_fastqc.html .
conda activate amplicon_env
# pandaseq
pandaseq -f 2_trimming/CT_A_1P.fastq -r 2_trimming/CT_A_2P.fastq -g 3_mergereads/CT_A.log -u 3_mergereads/CT_A.unaligned.fasta \
-w 3_mergereads/CT_A.fastq -F
# make a pandasq loop to merge all the files
ls 2_trimming/*_1P.fastq | sed s/_1P.fastq// | while read line; do id=$(echo $line | sed "s/.*\///"); pandaseq -f ${line}_1P.fastq -r ${line}_2P.fastq -g 3_mergereads/${id}.log -u 3_mergereads/${id}.unaligned.fasta -w 3_mergereads/${id}.fastq -F; done
# 17 What is the average length of the merged sequences?
awk 'NR==2 {print length; exit}' CT_A.fastq # single read counted
# for all reads in a file
awk '{if(NR%4==2) {count++; sum+=length($0)}} END {print sum/count}' *.fastq
# 253.015
# Create a new folder called 4_precluster
mkdir 4_precluster
# Vsearch (dereplication)
vsearch --fastx_uniques 3_mergereads/CT_A.fastq --strand plus --sizeout --relabel CT_A --fasta_width 0 \
--fastaout 4_precluster/CT_A.derep.fasta
# make a vsearch loop to dereplicate all the files
ls 3_mergereads/*.fastq | sed s/.fastq// | while read line; do id=$(echo $line | sed "s/.*\///"); vsearch --fastx_uniques ${line}.fastq --strand plus --sizeout --relabel ${id} --fasta_width 0 --fastaout 4_precluster/${id}.derep.fasta; done
# 21. What percentage of the reads did you discard during dereplication (not how many your kept)? 
# Make sure you understand how the reads are represented before and after dereplication.
# number of reads in the files before dereplication
awk 'NR==2 {print length; exit}' 3_mergereads/CT_A.fastq # single read counted
# for all reads in a file
awk '{if(NR%4==2) {count++}} END {print count}' 3_mergereads/CT_A.fastq
# 914006
grep -v ">" 4_precluster/CT_A.derep.fasta | wc -l # count the number of reads in the files and assign the result to surviving_reads variable
# 80214
echo "scale=5; (1 - 80214 / 914006) * 100" | bc
# 91.22400
# 22. calculate for the rest of the files
grep -v ">" 4_precluster/*.derep.fasta | wc -l # count the number of reads in the files and assign the result to surviving_reads variable
# calculate the percentage of reads discarded
# 580341
# Use the cat command to produce a single file containing all the dereplicated FASTA-sequences called all_combined.fasta
cat 4_precluster/*.derep.fasta > all_combined.fasta
# 22 Verify that the resulting file contains the same number of reads and nucleotides as the original files added together. 
# How many reads are in the resulting file?
grep -v ">" all_combined.fasta | wc -l # count the number of reads in the files and assign the result to surviving_reads variable
# 580341
# Dereplicate the reads using Vsearch
vsearch --derep_fulllength 4_precluster/all_combined.fasta --sizein --sizeout --fasta_width 0 --uc 4_precluster/all_derep.uc --output 4_precluster/all_derep.fasta
# 24. How many reads does the most frequent sequence from the Vsearch output represent?
cat 4_precluster/all_derep.fasta | grep 'size' | cut -d "=" -f 2 | sort -nr | head -n 1
# 181549
# 25 Calculate the number of reads from each sample. This can be achieved by looping over the sample names and extracting corresponding size counts for each. Make sure they correspond to the previous counts. Write down the number of reads per sample.
# script to count the number of reads
grep ">" 4_precluster/all_derep.fasta | grep "CT_A" | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum}'
# 3012197
grep ">" 4_precluster/all_derep.fasta | grep "CT_B" | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum}'
# 406046
grep ">" 4_precluster/all_derep.fasta | grep "CT_C" | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum}'
# 29568
grep ">" 4_precluster/all_derep.fasta | grep "CT_D" | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum}'
# 47582
grep ">" 4_precluster/all_derep.fasta | grep "MA_A" | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum}'
# 315366
grep ">" 4_precluster/all_derep.fasta | grep "MA_B" | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum}'
# 70644
grep ">" 4_precluster/all_derep.fasta | grep "MA_C" | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum}'
# 7516
# Precluster the reads using Vsearch
vsearch --cluster_size 4_precluster/all_derep.fasta --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --uc 4_precluster/all_preclust.uc --centroids 4_precluster/all_preclust.fasta
# 27. How many FASTA entries do we have after preclustering?
grep ">" 4_precluster/all_preclust.fasta | wc -l
# 48148
# 28. What is the number of centroids in the file?
grep '^S' 4_precluster/all_preclust.uc | wc -l
# 48148
# creating a directory called 5_chimera
mkdir 5_chimera
# Identify chimeras using Vsearch
vsearch --uchime3_denovo 4_precluster/all_preclust.fasta --threads 4 --sizein --sizeout --fasta_width 0 --nonchimeras 5_chimera/all.denovo.nonchimeras.fasta \
--chimeras 5_chimera/all.denovo.chimeras.fasta --uchimeout 5_chimera/all.denovo.uchime
# 32. What percentage of the reads was classified as chimeric?
# Count the total number of reads (sequences) in the input file
total_reads=$(grep -c '^>' 4_precluster/all_preclust.fasta)
grep -c '^>' 4_precluster/all_preclust.fasta
# 48148
# Count the number of chimeric reads
chimeric_reads=$(grep -c '^>' 5_chimera/all.denovo.chimeras.fasta)
grep -c '^>' 5_chimera/all.denovo.chimeras.fasta
# 15376
# Count the number of nonchimeric reads
nonchimeric_reads=$(grep -c '^>' 5_chimera/all.denovo.nonchimeras.fasta)
grep -c '^>' 5_chimera/all.denovo.nonchimeras.fasta
# 32722
# Calculate the percentage of chimeric reads
echo "15376/48148 * 100" | bc -l
# 31.937
scp -r /Users/med-snt/Downloads/silva.gold.bacteria.zip/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon
# unzip the file, trim away ‘-’ and ‘.’ characters
unzip silva.gold.bacteria.zip
unzip -p silva.gold.bacteria.zip silva.gold.align | sed -e "s/[.-]//g" > gold.fasta
# 33. Inspect the database. How many sequences does it contain? Approximately how long are the sequences?
cat gold.fasta | grep '^>' | wc -l
# 5181
cat gold.fasta | grep '^>' | wc -c
# 66450
# we are ready to perform the reference-based checking
vsearch --uchime_ref 5_chimera/all.denovo.nonchimeras.fasta --db gold.fasta --sizein --sizeout --fasta_width 0 --nonchimeras 5_chimera/all.ref.nonchimeras.fasta --chimeras 5_chimera/all.ref.chimeras.fasta --uchimeout 5_chimera/all.ref.uchime
less all.ref.uchime
# Check if the numbers are correct (for both chimera checking steps):
# Find out the column number for chimeric (Y) or not (N)
cat all.ref.uchime | head -1 | tr "\t" "\n" | cat -n 
cat all.ref.uchime | cut -f 18 | sort | uniq -c 
#   5 ?
# 28702 N
#  4065 Y
cat all.ref.chimeras.fasta | grep -c \>
# 4065
cat all.ref.nonchimeras.fasta | grep -c \>
# 28702
# 35 How many reads were discarded (in percentage) in each of the chimera checking steps?
echo "4065/32722 * 100" | bc -l
# 12.422
# retrieve all reads that are part of preclusters classified as not chimeric and save them in a directory called 6_clustering.E
mkdir 6_clustering
scp -r /Users/med-snt/Downloads/map.pl/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon
# Extract all reads belonging to non-chimeric preclusters
# all_derep.fasta: all reads dereplicated across colonies
# all_preclust.uc: information about the clustering step
# all.ref.nonchimeras.fasta: preclusters classified
#   as non chimeric
# make map.pl executable
chmod +x map.pl
perl map.pl 4_precluster/all_derep.fasta 4_precluster/all_preclust.uc 5_chimera/all.ref.nonchimeras.fasta > 6_clustering/all.nonchimeras.derep.fasta
# Next, we extract all non-chimeric reads
# all_combined.fasta: the starting set of reads
# all_derep.uc: information about the dereplication step
# all.nonchimeras.derep.fasta: dereplicated sequences
#   classified as non-chimeric (from previous step)
perl map.pl 4_precluster/all_combined.fasta 4_precluster/all_derep.uc 6_clustering/all.nonchimeras.derep.fasta > 6_clustering/all.nonchimeras.fasta
# 39. How many non-chimeric preclusters, dereplicated reads and raw reads do we have?
# Counting non-chimeric preclusters 
grep -c "^>" 5_chimera/all.ref.nonchimeras.fasta
# 28702
# Counting dereplicated reads
grep -c "^>" 6_clustering/all.nonchimeras.derep.fasta
# 387970
# Counting raw reads
grep -c "^>" 6_clustering/all.nonchimeras.fasta
# 473532
# 40. For each step, compare the percentage of the total number of reads, dereplicated reads and preclusters that are classified as chimeric.
# preclusters % for chimera (4_precluster/all_preclust.fasta:5_chimera/all.ref.nonchimeras.fasta)
grep -c "^>" 4_precluster/all_preclust.fasta
# 48148
echo "28702/48148 * 100" | bc -l
echo "100 -59.6" | bc -l
# 40.4
# dereplicated reads % for chimera (6_clustering/all.nonchimeras.derep.fasta:4_precluster/all_derep.fasta)
grep -c "^>" 4_precluster/all_derep.fasta
# 480652
echo "387970/480652 * 100" | bc -l
echo "100 -80.72" | bc -l
# 19.28
# raw reads % for chimera (6_clustering/all.nonchimeras.fasta:4_precluster/all_combined.fasta)
grep -c "^>" 4_precluster/all_combined.fasta
# 580341
echo "473532/580341 * 100" | bc -l
echo "100 -81.6" | bc -l
# 18.4
# rename the FASTA headers in the final output
sed "s/[0-9]\+;/;/" 6_clustering/all.nonchimeras.fasta > 6_clustering/all.nonchimeras.renamed.fasta
# Perform clustering using Vsearch
vsearch --cluster_size 6_clustering/all.nonchimeras.renamed.fasta --threads 4 --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --relabel OTU_ --mintsize 3 --uc 6_clustering/final.uchime --centroids 6_clustering/otus.fasta --otutabout 6_clustering/otus.tsv
# 42 Calculate the total number of OTUs and how many sequences they represent together. Does this correspond to the number of input sequences?
grep -c "^>" 6_clustering/otus.fasta
# 54304
awk 'NR>1 {for(i=2; i<=NF; i++) sum+=$i} END {print sum}' 6_clustering/otus.tsv
# 3679150
# 43. How many OTUs were identified in each of the seven colonies? (Hint: investigate the otus.tsv file)
awk 'NR>1 && $2!=0 {count++} END {print count}' 6_clustering/otus.tsv # CT_A 
# 9786
awk 'NR>1 && $3!=0 {count++} END {print count}' 6_clustering/otus.tsv # CT_B
# 5971 
awk 'NR>1 && $4!=0 {count++} END {print count}' 6_clustering/otus.tsv # CT_C
# 2403
awk 'NR>1 && $5!=0 {count++} END {print count}' 6_clustering/otus.tsv # CT_D
# 3247
awk 'NR>1 && $6!=0 {count++} END {print count}' 6_clustering/otus.tsv # MA_A
# 28349
awk 'NR>1 && $7!=0 {count++} END {print count}' 6_clustering/otus.tsv # MA_B
# 10357
awk 'NR>1 && $8!=0 {count++} END {print count}' 6_clustering/otus.tsv # MA_C
# 1347
# 44 What was the average number of sequences per OTU in each of the seven colonies?
awk 'NR>1 && $2!=0 {sum+=$2} END {print sum}' 6_clustering/otus.tsv # CT_A
# 871553
# average number of sequences per OTU in CT_A
echo "871553/5971" | bc -l
# 146
awk 'NR>1 && $3!=0 {sum+=$3} END {print sum}' 6_clustering/otus.tsv # CT_B
# 667450
# average number of sequences per OTU in CT_B
echo "667450/5971" | bc -l
# 112
awk 'NR>1 && $4!=0 {sum+=$4} END {print sum}' 6_clustering/otus.tsv # CT_C
# 360091
# average number of sequences per OTU in CT_C 
echo "360091/2403" | bc -l
# 150
awk 'NR>1 && $5!=0 {sum+=$5} END {print sum}' 6_clustering/otus.tsv # CT_D
# 277965
# average number of sequences per OTU in CT_D
echo "277965/3247" | bc -l
# 86
awk 'NR>1 && $6!=0 {sum+=$6} END {print sum}' 6_clustering/otus.tsv # MA_A
# 1022364
# average number of sequences per OTU in MA_A
echo "1022364/28349" | bc -l
# 36
awk 'NR>1 && $7!=0 {sum+=$7} END {print sum}' 6_clustering/otus.tsv # MA_B
# 438329
# average number of sequences per OTU in MA_B
echo "438329/10357" | bc -l
# 42
awk 'NR>1 && $8!=0 {sum+=$8} END {print sum}' 6_clustering/otus.tsv # MA_C
# 41398
# average number of sequences per OTU in MA_C
echo "41398/1347" | bc -l
# 31
# Create a directory called 7_classify.
mkdir 7_classify
scp -r /Users/med-snt/Downloads/rdp_classifier_2.14.zip/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon
unzip rdp_classifier_2.14.zip 
# 46. How many sequences are in it?
/home/inf-51-2023/Amplicon/rdp_classifier_2.14/data
grep -c "^>" rdp_classifier_2.14/data/ref_seqs.fa
# around three millions
# install the RDP classifier in Amplicon folder (amplicon_env)
conda install -c bioconda rdp_classifier
rdp_classifier -h
which rdp_classifier
# Classify the OTUs using the RDP classifier
rdp_classifier -c 0.8 -f fixrank -o 7_classify/all.fixedRank -h 7_classify/all.significant 6_clustering/otus.fasta
head -1 all.fixedRank | tr "\t" "\n" | cat -n # to see the columns
# 48 How many cluster sequences could not be classified as bacterial?
# 18806
cat all.fixedRank | awk -F "\t" '$3 != "Bacteria" {split($1,a,"="); sum+=a[2]} END {print sum}'

cat all.fixedRank | cut -f 9 | sort -n | uniq -c
head all.fixedRank | tr "\t" @ # to see the columns
# be determined on the phylum level (on the 0.8 confidence level)?
awk -F "\t" ' {if ($8>=0.8) ++i} END {print i} ' all.fixedRank
# 31061

cat all.significant | grep -P "\tphylum\t" | sed "s/.*\t//" | tr "\n" + | sed "s/+$/\n/" | bc
# 31061
# 49 How many cluster sequences were determined on the genus level (using 0.8)?
awk -F "\t" ' {if ($20>=0.8) ++i} END {print i} ' all.fixedRank
# 12027
# 50 How would it change if you had used thresholds of 0.7 and 0.9?
awk -F "\t" ' {if ($20>=0.7) ++i} END {print i} ' all.fixedRank
# 0.7 => 15326
awk -F "\t" ' {if ($20>=0.9) ++i} END {print i} ' all.fixedRank
# 0.9 => 8376
scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon/7_classify/all.fixedRank .
# 51. plot graph 
# print all the sequences that classified as bacterial
cat all.fixedRank | awk -F "\t" '$3 == "Bacteria" {print $0}' > all.bacterial
# be determined on the phylum level (on the 0.8 confidence level) and print out, no counting?
awk -F "\t" ' {if ($8>=0.8) print $0 }' all.bacterial | cut -f 1,6 > all.bacterial.phylum1
# remove the size information from all.bacterial.phylum
sed -r 's/;size=[0-9]+//g' all.bacterial.phylum1 > new_file1

scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon/7_classify/new_file1 .
scp inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/Amplicon/6_clustering/otus.tsv .

