ssh inf-51-2023@bioinf-serv2.cob.lu.se
q# Download BLAST
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
scp -r /home/inf-51-2023/Downloads/Prodigal-GoogleImport/ inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/bin
ls ncbi-blast-2.15.0+/bin/
# In the server bin directory
ls bin/ncbi/ncbi-blast-2.15.0+/bin
# In the home directory
~/.bashrc # permission denied
less .bashrc # check the content of the file
# don't know if we can modify things in the bashrc file
# we create an environment variable for the blast, using the conda
conda create -n blast_env
conda activate blast_env
conda search -c bioconda blast # check the channel
conda install -c bioconda blast # install blast
# Database installation
cd
mkdir BlastDB
cd BlastDB
update_blastdb.pl --showall # what NCBI offers and path to the database
# which: no gsutil in (/home/inf-51-2023/miniconda3/envs/blast_env/bin:/home/inf-51-2023/miniconda3/condabin:/home/inf-51-2023/.local/bin:/home/inf-51-2023/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin)
# which: no gcloud in (/home/inf-51-2023/miniconda3/envs/blast_env/bin:/home/inf-51-2023/miniconda3/condabin:/home/inf-51-2023/.local/bin:/home/inf-51-2023/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin)
# Download pdbaa
update_blastdb.pl pdbaa
ls
# the md5sum of this file should be the same as in pdbaa.tar.gz.md5. Test
md5sum -c pdbaa.tar.gz.md5
# extract the database
tar -xvzf pdbaa.tar.gz
# fasta file retrieved
blastdbcmd -outfmt %f -db pdbaa -entry all -out pdbaa.fas # -outfmt %f: fasta format
less pdbaa.fas # have a look
# Make your own BLAST database
# copy paxillus file from resources to BlastDB
cp /resources/binp28/Data/paxillus.faa.gz .
cp /resources/binp28/Data/paxillus.fna .
# unzip the file
gunzip paxillus.faa.gz
# create your own database from any fasta file
makeblastdb -in paxillus.fna -out paxillusNt -dbtype nucl
ls -l paxillusNt.n* # n stands for nucleotide
# do it for the protein file
makeblastdb -in paxillus.faa -out paxillusAa -dbtype prot
ls -l paxillusAa.p* # p stands for protein
makeblastdb -h # For short help
makeblastdb -help # For extensive help
# 2.5 BLAST example usage
cd .. # Move up to the parent directory
mkdir BlastExercise
cd BlastExercise
# cope the file from resources to BlastExercise
cp /resources/binp28/Data/yeast.faa.gz . 
# unzip the file
gunzip yeast.faa.gz
# -evalue: expect value 
blastp -query yeast.faa -db ../BlastDB/paxillusAa -evalue 1e-10 -out yeast_vs_Paxillus.blastp -num_descriptions 10  -num_alignments 5 -num_threads 2
# find the threads
nproc
# compare the input and output files, query files (>) and output files (“Query=”) has to be matched
# Check how many queries there are.
cat yeast.faa | grep -c '>'
# 6717
# Check processed queries, in the output,
# each query will start with 'Query='
cat yeast_vs_Paxillus.blastp | grep -c 'Query='
# 6717 
# look at the output file
less yeast_vs_Paxillus.blastp
# No hits found
# BLAST Command Line Applications User Guide
# https://www.ncbi.nlm.nih.gov/books/NBK279690/
# 1. How many of the yeast proteins had hits in the Paxillus proteome?
# 0
# Paxillus involutus is a basidiomycete. Another basidiomycete is Phanerochaete chrysosporium.
# copy the file from resources to BlastExercise
cp /resources/binp28/Data/Phchr2.faa.gz . 
# unzip the file
gunzip Phchr2.faa.gz
# copy the file from resources to BlastExercise
cp /resources/binp28/Data/cytB.faa.gz . 
cp /resources/binp28/Data/cytB.fna.gz .
# unzip the file
gunzip cytB.faa.gz
gunzip cytB.fna.gz
# segment of cytochrome B is encoded on the mitochondrial genome ‘mtDNA’ not the nuclear genome
cat Phchr2.faa | grep -c '>'
# 11777
cat Phanero_vs_Paxillus.blastp | grep -c 'Query='
# 11777
less Phanero_vs_Paxillus.blastp
# 2. How many Phanerochaete proteins produce significant hits against the
# Paxillus database? Motivate your answer.
cat Phanero_vs_Paxillus.blastp | grep -c 'No hits found' # 4464
cat yeast_vs_Paxillus.blastp | grep -c 'No hits found' # 3175
# calculate the difference
echo "$( echo 11777-4464 | bc -l )" # 7313 hits
# 3. Which of the two proteomes (yeast and Phanerochaete) had most hits
# against Paxillus. (Use an evalue of 1e-10)
# Phanerochaete has most hits against Paxillus 
cat Phanero_vs_Paxillus.blastp | grep -c 'Expect = 1e-10' # 113
cat yeast_vs_Paxillus.blastp | grep -c 'Expect = 1e-10' # 48
# 4. Is the result expected? 
# yes. 
# two databases consisting of one sequence each, DNA in one and protein in one
# Nucleotide:
makeblastdb -in cytB.fna -out nt -dbtype nucl
# Protein
makeblastdb -in cytB.faa -out aa -dbtype prot
# Take a look
ll nt* aa*
# https://sequenceserver.com/blog/choosing-blast-algorithms
blastn -query cytB.fna -db nt # Prints to stdout which is ok for now (nucleotide to nucleotide)
blastp -query cytB.faa -db aa # Prints to stdout which is ok for now (protein to protein)
tblastn -query cytB.fna -db nt # Prints to stdout which is ok for now (translated nucleotide to nucleotide)
tblastx -query cytB.fna -db nt # Prints to stdout which is ok for now (translated nucleotide to protein with translated query)
blastx -query cytB.fna -db aa # Prints to stdout which is ok for now (nucleotide to protein)
# 5. Explain the result for tblastx above

#InterPro on the web
# https://www.ebi.ac.uk/interpro/
# 1. Which domains did you find?
# RNA polymerase (8)
# 2. What is the function of the domains?
# RNA polymerase Rpb1, domain 1
# DNA binding (GO:0003677)
# DNA-directed 5'-3' RNA polymerase activity (GO:0003899)
# 3. Click on the first domain to go to its details page.
# https://www.ebi.ac.uk/interpro/entry/InterPro/IPR007080/

# 3.3 InterPro from the shell using InterProScan
mkdir Annotation
# copy file to current Annotation directory
cp /home/inf-51-2023/BlastDB/paxillus.faa . 
# look at the first 10 sequences (first 78 lines)
head -78 paxillus.faa > 10.faa
# Interproscan does not like * (stop codons)
# so you have to get rid of them before
# annotate domains for ten proteins from the Paxillus involutus proteome
cat 10.faa | sed 's/*//g' > 10_nostop.faa
# Launch the InterPro annotation analysis
interproscan.sh -i 10_nostop.faa -cpu 2 -d paxillus_interpro -goterms # -cpu 2: number of threads, -d: output directory, -goterms: GO terms
# not able to run (command line above), so look at the output (command line below)
cp -r /resources/binp28/Data/paxillus_interpro .
# Check the log file for details
less paxillus_interpro/10_nostop.faa.log
# copy the file from resources to Annotation
cp -r /resources/binp28/Data/pi.blastp.gz .
cp -r /resources/binp28/Data/pi.pfam.gz .
# unzip the file
gunzip pi.blastp.gz
gunzip pi.pfam.gz
# 4. Looking at the beginning of the blastp file, which database was searched?
less pi.blastp
# Database: UniProt
# Sometimes, gene predictors can over-predict and identify ORFs that are NOT in
# fact protein-coding genes. If a protein/gene does not have a significant hit using
# BLAST or PFAM, it might be a false positive protein/gene prediction.
# 5. How many proteins did not have a significant hit in the BLASTP file?
# (Remember how to see if a query protein did NOT have a hit? Revisit the BLAST exercise.)
cat pi.blastp | grep -c 'No hits found' | cut -d ";" -f 6 
# 4698
cat pi.blastp | grep "No hits" -B 5  | grep "Query=" | cut -d "=" -f 2 > no_hits.txt
# 6. Using bash, determine how many of the proteins without a BLASTP hit had a PFAM hit.
less pi.pfam
cat pi.pfam | grep -v \# | grep -w -f no_hits.txt | awk '{print $1} '| sort | uniq | wc -l # 239 (sort only the first field) -w: match whole word, -f: file, \# : # is a special character, so we need to escape it
cat pi.pfam | grep -v \# | grep -w -f no_hits.txt | sort | uniq | wc -l # 254 (may have duplicates)
# This is my own playing around (not correct)
awk '{print $1, $6, $13, $12}' pi.pfam | wc -l
# 14026 (total number of hits)
awk '$13 <= 1e-10 {print $1, $6, $13, $12}' pi.pfam | wc -l # 8624 (total number of significant hits)
echo "$( echo 14026-8624 | bc -l )" # 5402 (total number of insignificant hits)
awk '$13 > 1e-10 {print $1, $6, $13, $12}' pi.pfam | wc -l # 5402 (total number of insignificant hits)
cat pi.blastp | grep -c 'No hits found' | cut -d ";" -f 6 
cat pi.blastp | grep -E 'Query=|No hits found' | grep -c 'No hits found' | cut -d ";" -f 6
# Extract protein IDs without a BLASTP hit
grep 'No hits found' pi.blastp | cut -d ";" -f 6 > no_blastp_hits.txt
# pfam website
# https://www.google.com/search?q=how+to+look+at+pfam+hit+&sca_esv=600117354&rlz=1C5CHFA_enSE890SE892&biw=1273&bih=738&tbm=vid&sxsrf=ACQVn09z21Jl12gft929ka07Y1TgFlGJQA%3A1705784592366&ei=EDWsZZL9FerOwPAP0L6CsAs&udm=&ved=0ahUKEwiSxcPl7uyDAxVqJxAIHVCfALYQ4dUDCA0&uact=5&oq=how+to+look+at+pfam+hit+&gs_lp=Eg1nd3Mtd2l6LXZpZGVvIhhob3cgdG8gbG9vayBhdCBwZmFtIGhpdCAyBBAjGCdI9QVQyQNYyQNwAHgAkAEAmAFqoAG0AaoBAzEuMbgBA8gBAPgBAYgGAQ&sclient=gws-wiz-video#fpstate=ive&vld=cid:bd6222ed,vid:unExXijjMmI,st:0
# https://pfam-docs.readthedocs.io/en/latest/faq.html#what-is-a-clan

# Extract protein IDs with a PFAM hit
awk '$13 <= 1e-10 {print $1}' pi.pfam > pfam_hits.txt

# Find proteins that are in the PFAM list but not in the BLASTP list
grep -Fvxf no_blastp_hits.txt pfam_hits.txt > proteins_with_pfam_but_no_blastp.txt

# Count the number of such proteins
wc -l proteins_with_pfam_but_no_blastp.txt

# Extract protein IDs without a BLASTP hit
grep -E 'Query=|No hits found' pi.blastp | grep 'No hits found' | cut -d ";" -f 6 > no_blastp_hits.txt
grep -E 'Query=|No hits found' pi.blastp > no_blastp_hits1.txt # -E: extended regular expression
# Extract protein IDs with a PFAM hit
awk '$13 <= 1e-10 {print $1}' pi.pfam > pfam_hits.txt
# Find proteins that are in the PFAM list but not in the BLASTP list
grep -Fvxf no_blastp_hits.txt pfam_hits.txt > proteins_with_pfam_but_no_blastp.txt
grep -Fvxf no_blastp_hits1.txt pfam_hits.txt > proteins_with_pfam_but_no_blastp1.txt 
# -F: fixed string, -v: invert match, -x: match whole line, -f: file
# Count the number of such proteins
wc -l proteins_with_pfam_but_no_blastp.txt
wc -l proteins_with_pfam_but_no_blastp1.txt
echo "$( echo 8624-4698 | bc -l )" # 3926 (total number of proteins with a PFAM hit but no BLASTP hit)
# convert a multiple line fasta file to one-line fasta files
cat paxillus.faa | tr "\n" @ | sed "s/>/>anyString/g" | tr ">" "\n" | sed "s/anyString/>/" | sed "s/@/\n/" | tr -d @ | grep . > pi.oneline.faa
grep -B 5 "No hits found" pi.blastp
grep -B 5 "No hits found" pi.blastp | grep -c "Query="
grep -B 5 "No hits found" pi.pfam 
# 4.2 SignalP
# extract the first 143 sequence from the paxillus.faa
head -1000 paxillus.faa > 143.faa # 143.faa is a subset of paxillus.faa
# run signalP on the 143 proteins
signalp -fasta 143.faa
# output file: 143_summary.signalp5
# https://services.healthtech.dtu.dk/services/SignalP-5.0/
# Column:
# 1. Your protein ID
# 2. The protein prediction, one of: - SP(Sec/SPI); - LIPO(Sec/SPII); -TAT(Tat/SPI); or - OTHER
# 3. The probability of column 2’s prediction
# 4. The probability of another prediction
# 5. Cleavage Site (CS) position and associated probability.
# 1. How many proteins are predicted to be secreted (SPI or SPII)?
cat 143_summary.signalp5 | grep -v ^# | grep SP | cut -f 1 > paxillus143.signal.ids # -v: invert match, ^#: start with # -f: field
# 4.3 TargetP: use targetp to predict subcellular localization
targetp -fasta 143.faa
# The output is as follows in 143_summary.targetp5
# https://services.healthtech.dtu.dk/service.php?TargetP-2.0
# less 143_summary.targetp2
# Column:
# 1. Your protein ID
# 2. The protein prediction, one of: - SP (secreted); - mTP (mitochondrial);
# - plant specific: cTP (chloroplastic); - plant specific: luTP (thylakoid
# luminal transit peptide); - noTP (no prediction)
# 3. The probability of noTP prediction
# 4. The probability of SP
# 5. The probability of mTP prediction
# 6. Cleavage Site (CS) position and associated probability.
# 4.3.1 SignalP/TargetP Exercise
# 2. From the targetp output, how many proteins do we have for secreted proteins (SP), Mitochondrial (mTP), Chloroplastic?
cat 143_summary.targetp2 | grep -v ^# | grep -c SP | cut -f 1 # 7
cat 143_summary.targetp2 | grep -v ^# | grep -c mTP | cut -f 1 # 7
cat 143_summary.targetp2 | grep -v ^# | grep -c cTP | cut -f 1 # 0
cat 143_summary.signalp5 | grep -v ^# | grep -c SP | cut -f 1 
# 3. Compare the confidence scores between those proteins predicted to be
# secreted by signalP (‘SPI’ or ‘SPII’) with their score with targetP. Are there disagreements? Similarities?
cat 143_summary.signalp5 | grep -v ^# | grep SP | cut -f 1 # 5
cat 143_summary.targetp2 | grep -v ^# | grep SP | cut -f 1 # 7
# sort the common proteins
comm -12 <(cat 143_summary.signalp5 | grep -v ^# | grep SP | cut -f 1 | sort) <(cat 143_summary.targetp2 | grep -v ^# | grep SP | cut -f 1 | sort)
# Paxin1_151859
# Paxin1_151862
# Paxin1_164887
# Paxin1_164896
# Paxin1_64778
# sort the unique proteins
cat <(cat 143_summary.signalp5 | grep -v ^# | grep SP | cut -f 1) <(cat 143_summary.targetp2 | grep -v ^# | grep SP | cut -f 1) | sort | uniq
# Paxin1_151859
# Paxin1_151862
# Paxin1_159431
# Paxin1_164887
# Paxin1_164896
# Paxin1_48533
# Paxin1_64778
# print the number of unique proteins in targetp but not in signalp
comm -23 <(cat 143_summary.targetp2 | grep -v ^# | grep SP | cut -f 1 | sort) <(cat 143_summary.signalp5 | grep -v ^# | grep SP | cut -f 1 | sort)
# Paxin1_159431
# Paxin1_48533
# The eggNOG database & emapper
# http://eggnog5.embl.de/#/app/methods
# https://academic.oup.com/nar/article/47/D1/D309/5173662
# https://en.wikipedia.org/wiki/Sequence_homology#Orthology
# https://github.com/eggnogdb/eggnog-mapper/wiki
# make a directory
mkdir emapper_tmp
# conda new environment
conda create -n eggnog --file \
/resources/binp28/Programs/eggNOG/newest/requirements.txt \
diamond -y
# don't forget to activate the eggnog env
conda activate eggnog
# then launch emapper
nohup emapper.py -m diamond -i paxillus.faa -o paxillus.emap \
--temp_dir emapper_tmp --cpu 2 2>&1 &
# work with ready made output
cp -r /resources/binp28/Data/emapper_out/ ./
# paxillus.emap.annotation is pretty big
less -S paxillus.emap.emapper.annotations
# Information about the columns
# https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.6#Output_format
head -n 30 paxillus.emap.emapper.annotations | cut -f 1-8,12 # -f: field
# Columns information:
# • query_name
# • seed_eggNOG_ortholog: the best hit against the eggNOG db (taxonID.
# accession number)
# • seed_ortholog_evalue: the evalue of that match
# • seed_ortholog_score: the score of that match
# • eggNOG OGs: a comma-separated, clade depth-sorted (broadest to narrowest),
# list of Orthologous Groups (OGs) identified for this query. Note that
# each OG is represented in the following format: OG@tax_id|tax_name;
# for any OG, you can take the OG name (before the @) and search for it
# on http://eggnog5.embl.de/#/app/home
# • narr_og_name: OG@tax_id|tax_name for the narrowest OG found for
# this query.
# • narr_og_cat: COG category corresponding to narr_og_name
# • narr_og_desc: Description corresponding to narr_og_name
# • Preferred_name: preferred name of the protein
# look at one particular entry Paxin1_122960 - remember we’ll just need to look at the first 8 columns
# 1. What is the tax_name of the most refined OG for this sequence?
cat paxillus.emap.emapper.annotations | grep Paxin1_122960 | cut -f 1-8,12
# Paxin1_122960   85982.XP_007323734.1    2.62e-306       846.0   COG0688@1|root,KOG2420@2759|Eukaryota,38E13@33154|Opisthokonta,3NZD1@4751|
# Fungi,3UYYI@5204|Basidiomycota,226V2@155619|Agaricomycetes      4751|Fungi      I       Catalyzes the formation of phosphatidylethanolamine (PtdEtn) from phosphatidylserine (PtdSer). 
# Plays a central role in phospholipid metabolism and in the interorganelle trafficking of phosphatidylserine        ko:K01613
# Most refined OG: 226V2@155619|Agaricomycetes
# 2. What is the predicted function of this protein? What is the preferred protein name?
# Plays a central role in phospholipid metabolism and in the interorganelle trafficking of phosphatidylserine 
# 3. Using your results from a previous exercise, where would you expect this
# protein to function in the cell? Does this align with what is known about
# this protein (try and find information about it on wikipedia or uniprot)?
cat paxillus.emap.emapper.annotations | grep Paxin1_122960 | cut -f 1-8,12
# https://www.genome.jp/entry/K01613
awk '$13 <= 1e-10 {print $1, $6, $13, $12}' pi.pfam | grep Paxin1_122960
# It aligned with the previous search on interpro. https://www.ebi.ac.uk/interpro/entry/pfam/PF02666/domain_architecture/
# PS_Dcarbxylase

# Quiz
# The most frequent Repeat found in Paxillus pfam annotations is WD40
sed 's/ \+ /\t/g' pi.pfam > tab.pfam # replace multiple spaces with tab
cat tab.pfam | grep -v \#  | cut -f 6 | sort | uniq -c | sort -nr | head # sort by column 6, count unique, sort by count, show the first 10 lines


