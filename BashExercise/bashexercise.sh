# Ctrl + R history
# ^ + a beginning of the sentence
# ^ + e the end of the sentence
# Command line patterns and rename the .fna files to .fasta
for line in *.fna ; do basename=$(echo "$line" | cut -d. -f 1); mv -i "$line" "$basename".fasta ; done
# ctrl c to quit
##################
animal=dog
if [[ $animal == "dog" ]]; then echo "It's a dog"; fi
animal=cat
if [[ $animal != "dog" ]]; then echo "It's not a dog";
fi
######################
touch 1.fastq
if [[ -e 1.fastq ]]; then echo " File already exists "; fi    
######################
if [[ $((2+3)) -eq 4 ]]; then echo " This is true "; else
echo " Not true "; fil
##################
chmod +x try.sh # to run a script
#################
find ~/bin -name "*.sh"
find ~/bin -maxdepth 1 -name "*.sh"
find ~/bin -maxdepth 1 -name "*.sh" -mtime -1
find ~/bin -maxdepth 1 -name "*.sh" -mtime +1
##################
# find all large files (>=100Mb) 
find ~ -not -name " *.gz " -size +100M
# compress the files
find ~ -not -name " *.gz " -size +100M -exec gzip {} \;
#
find ~/bin -exec echo -n " md5sum: " \; -exec md5sum {} \;
# run script average script
./average.sh 10 20 40
##############
# calculate GC content
# 4.8:
grep -v '^>' paxillus.fna | while read -r seq ; do gc=$(tr -dc 'GC' <<< "$seq" | wc -m) ; printf '100*%s/%s\n' "$gc" "${#seq}" | bc -l ; done
#######################33
echo $(($(echo -n Höör | wc -c)+ $(echo -n Hoor | wc -c)))
Investigate why. Compare it with:
echo $(($(echo -n Höör | wc -m)+ $(echo -n Hoor | wc -m)))
# 10 and 8 respectively, because the ö character is 2 bytes (in Swedish UTF-8)
#########################
awk ' BEGIN { print " Hello world " } '
awk ' BEGIN { print 4+5} '
awk ' BEGIN { print sin (3.14/6) } '
awk ' BEGIN { print sin ( atan2 (0 , -1) /6) } ' # arctangent of 0/ -1
awk ' BEGIN { a =2+3; b =4/8; print a , b } '
# OFS stands for output field separator
awk ' BEGIN { OFS = " xxxxx " ; a =2+3; b =4/8; print a , b } '
awk ' BEGIN { OFS = " \ n " ; a =2+3; b =4/8; print a , b } '
#############
awk ' /^>/ {print} ' regions.fna
awk ' /^>/ {print $0} ' regions.fna
awk ' /^>/ {print} ' regions.fna
awk ' /^>/ {print $1} ' regions.fna
# count the number of sequences from regions 1 (FPV3KKJ01) and 2 (FPV3KKJ02)
awk ' /^>FPV3KKJ01/ {countRegion1++} /^>FPV3KKJ02/
{countRegion2++} END {print countRegion1 ,
countRegion2} ' regions.fna
# count the number of sequences from regions 1 (FPV3KKJ01)
awk ' /^>FPV3KKJ01/ {countRegion1++} {print countRegion1} ' regions.fna
# To make a more readable output
awk '/^>FPV3KKJ01/ {countRegion1++} /^>FPV3KKJ02/ {countRegion2++} END {print "Region 1:" countRegion1 "\n" "Region 2:" countRegion2}' regions.fna
#  output a list of ids (without the larger than sign) and the lengths
awk ' /^>/ {sub( ">" ," " , $1) ; print $1} ' regions.fna
# The gsub (global substitution) command replaces every occurrence
awk ' /^[^>]/ {gsub(/[ AGC]/, "", $1); print $1}' regions.fna
# calculate the length of a string
awk ' /^[^>]/ {print length($0)} ' regions.fna
# make the sequence in lowercase
awk ' /^>/ {print} /^[^>]/ {print tolower($0)} ' regions.fna
# look for bashrc_profile
find ~ -name " .bashrc_profile "
cd
velvetg
velveth

