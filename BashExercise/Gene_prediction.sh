ssh inf-51-2023@bioinf-serv2.cob.lu.se
# Installation Prodigal
# https://github.com/hyattpd/Prodigal
# file downloaded at local 51-2023, paste on server
scp -r /home/inf-51-2023/Downloads/Prodigal-GoogleImport/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin
# On server
cd ~/bin
# unzip Prodigal-GoogleImport.zip
cd Prodigal-GoogleImport
make # to compile the program
# Installation Aragon
# http://www.ansikte.se/ARAGORN/
# download file using curl 
curl -O http://www.ansikte.se/ARAGORN/Downloads/aragorn1.2.41.c
# compiling
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.41.c
# Download manual
curl -O http://www.ansikte.se/ARAGORN/Downloads/aragorn.1
# make directory
mkdir -p ~/GenePrediction/Bacteria
# change directory in server
cd /
# copy file from server to local 
scp inf-51-2023@bioinf-serv2.cob.lu.se:/resources/binp28/Data/MyFirstAssembly_LargeContigs_out_AllStrains.unpadded.fasta .
# rename a file
mv MyFirstAssembly_LargeContigs_out_AllStrains.unpadded.fasta geo.fna
# move file to Bacteria directory in server
scp -r /Users/med-snt/GenePrediction/Bacteria/geo.fna/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/BacterialGenomeAssembly/
# location of the prodigal
/home/inf-51-2023/bin/Prodigal-GoogleImport/prodigal 
./prodigal --help # at the same folder

# how to predict genes in a genome using prodigal
prodigal -i geo.fna -a geo.faa -d geo.fna -f gff -o geo.gff -p meta
Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
                 [-p mode] [-q] [-s start_file] [-t training_file] [-v]

# 2. What is the purpose of the -c and -m options?
# -c: Closed ends. Do not allow genes to run off edges. -m: Treat runs of N as masked sequence; don't build genes across them.
# -g: Specify a translation table to use (default 11).
# To run prodigal (two times)
/home/inf-51-2023/bin/Prodigal-GoogleImport/prodigal -i geo.fna -a geoProteins.faa -d geoGenes.fna -f gff -o geo.gff -t geoTraining -g 11 -c -m
# 3. check to find the nonstandard nuecleotide
cat geoGenes.fna | grep -v \> | tr -d "\n" | tr -d "ACTG" | wc -c # 16
cat geoGenes.fna | grep -v \> | tr -d "\n"acgtACGT | wc -c #16
# get rid of the # in the file
grep -v ^# geo.gff
# 4. frequency of the start codon (ATG, TTG, GTG, CTG)
grep -c '>' geoGenes.fna # 3521
g=$(cat geo.gff | grep -v "^#" | cut -f 9 | cut -d ";" -f 3 | sort | uniq)
for f in $g; do
    echo "$f"
    grep -v "^#" geo.gff | grep -c "$f"
done
# You can calculate the frequency by dividing each number by the total number of genes 3521
cat geo.gff | grep -v ^# | grep -c "start_type=ATG" #2779
cat geo.gff | grep -v ^# | grep -c "start_type=TTG" # 284
cat geo.gff | grep -v ^# | grep -c "start_type=GTG" # 458
cat geo.gff | grep -v ^# | grep -c "start_type=CTG" # 0
# 5. Is the number of genes similar in your gene prediction to that of the reference genome
cat geoGenes.fna | grep -c ">My" | tr -d "\n" # 3521, not the same as the reference genome (Mine: 3650) 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009785.1/
# Teacher: Check The reference genome, https://www.ncbi.nlm.nih.gov/genome/1659,
# has 3584 genes, while the gene prediction shows 3521 genes.

# 6. How many sequences are there in the forward strand and reverse strand
cat geo.gff | grep -v "^#" | cut -f7 | sort | uniq -c
awk '$7 == "+" ' geo.gff | wc -l # single ' is used for the command line by awk # 1639 forward strand
awk '$7 == "-"' geo.gff | wc -l # single ' is used for the command line by awk # 1882 reverse strand
# 7. Do the genes have higher or lower GC content as compared to the genome?
# In the same .gff file, it says that the GC content of the genome is
# 52.04%. The average GC content of the genes is 52.37%. This can
# be obtained by counting the sum of all the values in the GC content
# column of the .gff file, then dividing by the number of genes (3521). 
# The genes appear to have a slightly higher GC content than the genome.
cat geo.gff | grep -v "#" | cut -f9 | cut -d ";" -f6 | cut -d "=" -f2 | paste -sd+ | bc # 1843.837
# 1843,837÷3521 = 0,523668560
cat geo.gff | grep -v "^#" | cut -d ";" -f 6 | cut -d "=" -f 2 | awk '{sum+=$1} END {print sum/NR}' # 0.523669
# 8. What is the gene density (how much of the genome consists of genes)?
# Notice that we do not include non-protein coding genes yet.
# In the Genome Assembly exercise the genome was calculated as being 3520016 bp long. 
# The total length of the predicted genes is 2999490 bp
cat geoGenes.fna | grep -v \> | tr -d "\n" | wc -c # 2999490
stats.sh geo.fna # 3.529 MB (total genome base pairs) conda activate assembly
stats.sh geoGenes.fna # 2.999 MB (total gene base pairs) conda activate assembly
echo "$( echo 2999000/3529000 | bc -l )" # 0.84981581184471521677

Installation of Aragon
scp -r /Users/med-snt/Downloads/aragorn1.2.41.c/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin
scp -r /Users/med-snt/Downloads/aragorn.1/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin
# compiling Install the tRNA predictor aragorn from:
# http://www.ansikte.se/ARAGORN/
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.41.c
aragorn -h
# https://manpages.ubuntu.com/manpages/trusty/man1/aragorn.1.html (manual for aragorn)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC373265/

# 9. How many tRNAs are predicted in your assembly? Look at the last lines in the output files.
# tail geo.trna, 170 sequences searched. Total tRNA genes = 90
aragorn -t geo.fna -o geo.tRNA # (tRNA mode) = aragorn -w geo.fna -o geo1.tRNA -t # 90
aragorn -w geo.fna -o geo1.tRNA # (batch mode)
# 10. Do you find at least one tRNA for each amino acid?
cat geo1.tRNA | grep -v \> | grep "tRNA-" | cut -d "-" -f 2 | cut -d " " -f1 | sort | uniq -c
cat geo1.tRNA | grep tRNA- | awk '{print $2 }'| sort | uniq -c | wc -l #20
# 11. How many anticodons do you find? Why not 61?
cat geo1.tRNA | grep -v \> | grep "tRNA-" | cut -d "(" -f 2 | sort | uniq | wc -l # 41
cat geo1.tRNA | grep tRNA- | awk '{print $5 }'| sort | uniq -c | wc -l # 41
# 12. Do you find any tmRNA (use -m)?
aragorn -w geo.fna -o geoM.tRNA -m # 1
# GeneMark
# output file gff
/resources/binp28/Data/genemark.gtf 
# copy file from server to server
cp /resources/binp28/Data/Paxin1_AssemblyScaffolds_Repeatmasked.fasta.gz . # (local server EukaryoticGenePrediction)
# unzip file
gunzip Paxin1_AssemblyScaffolds_Repeatmasked.fasta.gz
# rename file
mv Paxin1_AssemblyScaffolds_Repeatmasked.fasta paxillusGenome.fna
# run stats
stats.sh paxillusGenome.fna # 1.000 MB (total genome base pairs) at assembly environment
# 15. How big is the genome and what GC content does it have?
# GC content 0.5019 from stats.sh
cat paxillusGenome.fna | grep -v \> | tr -cd cgCG | wc -c 
# 24642235 (total GC base pairs)
# 16. 454 paired end sequencing has been used (in Illumina corresponding to
# mate pairs). What percentage of the assembled genome sequence is not known (containing Ns)?
# 15.781% gap
# Working with GFF files
# http://gmod.org/wiki/GFF2
# 17. How many genes were predicted in scaffold 1?
cat /resources/binp28/Data/genemark.gtf | cut -d$'\t' -f 3 | grep "gene" -c # 758
cat /resources/binp28/Data/genemark.gtf | grep "scaffold_1" | cut -f 3 | grep -c "gene"
# 758
cat /resources/binp28/Data/genemark.gtf | grep -v ^# | grep -c "scaffold_1" 
cat /resources/binp28/Data/genemark.gtf | grep -c "scaffold_1" | cut -f 3 

# 18. How many exons are there in average in scaffold 1?
grep -v '^#' /resources/binp28/Data/genemark.gtf | awk '$3=="gene"{gene++} $3=="CDS"{cds++} END{printf "%.2f\n", cds/gene}'
# 6.99
cat /resources/binp28/Data/genemark.gtf | grep "scaffold_1" | cut -f 3 | grep -c "CDS"
# 5301
 echo "$( echo 5301/758 | bc -l )" # 6.9934036939
 19. Optional (difficult): What is the average intron size in scaffold 1?
awk '$3 == "intron" {intron_length+=$5-$4+1; intron_count++} END {print intron_length/intron_count}' /resources/binp28/Data/genemark.gtf
# 97.5036
cat /resources/binp28/Data/genemark.gtf | grep "scaffold_1" | cut -f 3 | grep "intron" | cut -f 4 5 |  |awk '{sum+=$1} END {print sum/NR}'
# calculate the difference between the two columns 4 and 5
cat /resources/binp28/Data/genemark.gtf | grep "scaffold_1" | grep "intron" | cut -f 4-5 | awk '{print $2-$1}' | awk '{sum+=$1} END {print sum/NR}'
# 96.5036
# 20. Run ARAGORN on the complete genome file and compare the number of tRNAs to that of “your” Geobacillus. Discuss gene dosage and expression
aragorn -w paxillusGenome.fna -o paxillus.tRNA -t
# 94 tRNAs
# copy gffParse.pl to bin
cp /resources/binp28/Data/gffParse.pl . # (local server bin)
# add chmod to make it executable
chmod +x gffParse.pl
gffParse.pl -h
#./gffParse.pl -h
#perl gffParse.pl -h
grep -m 2 -n scaffold paxillusGenome.fna # To get line numbers
head -56224 paxillusGenome.fna > scaffold1.fna
# Alternatively with awk:
awk ' /^>/ {++i} {if (i==2) {exit}} {print} ' paxillusGenome.fna > scaffold1.fna
# -m 2 tells grep to only report the two first matching lines while -n makes grep
# to also output the current line number.

gffParse.pl -i scaffold1.fna -g /resources/binp28/Data/genemark.gtf -f CDS \ -b paxillus -p -a gene_id
# How to know Ribosomal binding site sequence of Geobacillus' predicted gene c23_25 is AGGAGG?
cat geoGenes.fna |  grep c23_25
# how to check operating system of the server
uname -a
# Linux bioinf-serv2.cob.lu.se 6.6.6-100.fc38.x86_64 #1 SMP PREEMPT_DYNAMIC Mon Dec 11 17:27:04 UTC 2023 x86_64 GNU/Linux
# Installation of GeneMark
scp -r /Users/med-snt/Downloads/gm_key_64/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/
# check the key 
ls ~/.gm_key
# simply calling ‘gmes_petap.pl
gmes_petap.pl
# In this program -h is not understood
nohup gmes_petap.pl --ES --sequence scaffold1.fna &
# Then hit return and you get the prompt back
top # check that it is running
u user # to list only processes that the user runs
c # display the entire command line
q # quit top
# To get the process in foreground again
fg
# To get it back to background
Ctrl-z # Pauses execution
bg # Resumes execution in background
# /home2/resources/binp28/Programs/Genemark/gmes_linux_64/ProtHint/bin/prothint.py:503: SyntaxWarning: invalid escape sequence '\.'
#  systemCall('sed \"s/\.//\" ' + args.proteins + ' | sed \"s/|/_/g\" > ' +
# error, file not found /usr/local/bin/gmes_petap.pl: scaffold1.fna
