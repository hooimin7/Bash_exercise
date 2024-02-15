# original article: https://www.nature.com/articles/s41559-017-0437-7#Sec2
# another article: https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03438-7#Sec14
# VCF format: https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
# Vcftools manuals: https://vcftools.sourceforge.net/man_latest.html
# bcftools manuals: https://samtools.github.io/bcftools/bcftools.html#view
# plink manuals: https://www.cog-genomics.org/plink/1.9/general_usage
# vcfR https://knausb.github.io/vcfR_documentation/visualization_1.html
# github tutorial: https://speciationgenomics.github.io/ADMIXTURE/
# ADMIXTURE tutorial: https://owensgl.github.io/biol525D/Topic_8-9/plotting_structure.html
# FST tutorial: https://www.youtube.com/watch?v=g9ftx6taogA
# plink     ~/miniconda3/bin/plink
# Admixutre    ~/miniconda3/bin/admixture
ssh UniServer
mkdir Research_project
cd Research_project
cp /resources/binp28/Data/ProjTaxa.vcf.gz .
# unzip the file
gunzip ProjTaxa.vcf.gz
# To analyse the file
vcftools --vcf ProjTaxa.vcf

less ProjTaxa.vcf
cat ProjTaxa.vcf | grep "^#CHROM" | tr "\t" "\n" | cat -n
#      1 CHROM
#      2  POS
#      3  ID
#      4  REF
#      5  ALT
#      6  QUAL
#      7  FILTER
#      8  INFO
#      9  FORMAT
#     10  8N05240
#     11  8N05890
#     12  8N06612
#     13  8N73248
#     14  8N73604
#     15  K006
#     16  K010
#     17  K011
#     18  K015
#     19  K019
#     20  Lesina_280
#     21  Lesina_281
#     22  Lesina_282
#     23  Lesina_285
#     24  Lesina_286
#     25  Naxos2

# check the chromosome
grep -v "^#" ProjTaxa.vcf | cut -f 1 | uniq 
# Identify chromosomes
awk '/^#/ {next} {print $1}' ProjTaxa.vcf | uniq

# print the site mean depth
grep -v "^#" ProjTaxa.vcf | cut -f 10 | cut -d ":" -f 3 | head -n 10

# Count SNPs
grep -v "^#" ProjTaxa.vcf | wc -l
# To count the number of SNPs in a VCF file
vcftools --vcf ProjTaxa.vcf --freq --out output
echo $(($(wc -l < output.frq) - 1))

# Count the number of individuals
vcftools --vcf ProjTaxa.vcf --depth
# Estimating mean coverage across all sites and individuals 
vcftools --vcf ProjTaxa.vcf --site-mean-depth 
rsync UniServer:/home/inf-51-2023/Research_project/out.ldepth.mean .

# remove the individuals Naxos2
vcftools --vcf ProjTaxa.vcf --remove-indv Naxos2 --recode --recode-INFO-all --out ProjTaxa_no_Naxos2
vcftools --vcf ProjTaxa_no_Naxos2.recode.vcf --site-mean-depth 
rsync UniServer:/home/inf-51-2023/Research_project/out.ldepth.mean .

# remove the minimum depth < 3
vcftools --vcf ProjTaxa_no_Naxos2.recode.vcf --min-meanDP 3 --recode --recode-INFO-all --out ProjTaxa_no_Naxos2_minDP3

# remove the maximum depth > 15
vcftools --vcf ProjTaxa_no_Naxos2_minDP3.recode.vcf --max-meanDP 15 --recode --recode-INFO-all --out ProjTaxa_no_Naxos2_minDP3_maxDP15

# remove the missing data
vcftools --vcf ProjTaxa_no_Naxos2_minDP3_maxDP15.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out ProjTaxa_no_Naxos2_minDP3_maxDP15_maxMiss0.5

# combine the commands
vcftools --vcf ProjTaxa.vcf --remove-indv Naxos2 --min-meanDP 3 --max-meanDP 15 --max-missing 0.8 --recode --recode-INFO-all --out ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss

# To count the number of SNPs in a VCF file after filtering
vcftools --vcf ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf --freq --out filtered_output
echo $(($(wc -l < filtered_output.frq) - 1))

# create a new environment
conda create -n research_project
conda activate research_project
conda install bioconda::plink

# print the first 10 lines of the file
grep -v "^#" ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf | head -10 
# SNP density per chromosome
vcftools --vcf ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf --SNPdensity 1000000

cat out.snpden | grep -E "^CHROM|chr5" > out.snpden_chr5 
cat out.snpden | grep -v "^CHROM" | grep "chr5" | wc -l
cat out.snpden_chr5 | wc -l
cat out.snpden | grep -E "^CHROM|chrZ" > out.snpden_chrZ
cat out.snpden | grep -v "^CHROM" | grep "chrZ" | wc -l
cat out.snpden_chrZ | wc -l

rsync UniServer:/home/inf-51-2023/Research_project/out.snpden_chr5 .
rsync UniServer:/home/inf-51-2023/Research_project/out.snpden_chrZ .

# make a plink directory
mkdir plink
# move into it
cd plink
VCF=~/Research_project/ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf
# perform linkage pruning - i.e. identify prune sites
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out cichlids
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out cichlids

# Perform PCA analysis
# prune and create pca
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract cichlids.prune.in \
--make-bed --pca --out cichlids
rsync UniServer:/home/inf-51-2023/Research_project/plink/cichlids.eigenvec .
rsync UniServer:/home/inf-51-2023/Research_project/plink/cichlids.eigenval .

# Install admixture
conda install bioconda::admixture

# Make a directory in the home directory
mkdir ADMIXTURE
cd ADMIXTURE

# Generate the input file in plink format
plink --vcf /home/inf-51-2023/Research_project/ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf --make-bed --out ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.bim > ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.bim.tmp
mv ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.bim.tmp ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.bim

# remove SNPs with more than 10% missing data
plink --bfile ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf --geno 0.1 --make-bed --out ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered 
# run Admixture with cross-validation (the default is 5-fold CV, for higher, choose e.g. cv=10) and K=2.
admixture --cv ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.bed 2 > log2.out
# run Admixture with cross-validation (the default is 5-fold CV, for higher, choose e.g. cv=10) and K=3-5.
for i in {3..5}
do
 admixture --cv ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.bed $i > log${i}.out
done
# ADMIXTURE produced 2 files: .Q which contains cluster assignments for each individual and .P which contains for each SNP the population allele frequencies.

# To check the CV error
grep -i 'CV error' log2.out
# extract the number of K and the CV error for each corresponding K
# Identify the best value of k clusters which is the value with lowest cross-validation error
# ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf = PnNmmmrv
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > PnNmmmrv_filtered.cv.error

# Extract the species name from the individual name
awk '{split($1,name,"."); print $1,name[2]}' ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.nosex > PnNmmmrv_filtered.list
awk '{print $2}' ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.nosex > PnNmmmrv_filtered.list
# add Lesina in front of *^2
sed -i '/^28/ s/^/Lesina_/' PnNmmmrv_filtered.list
# add species name to the list
awk 'BEGIN{OFS="\t"} /^8/{print $0, "species1"} /^K/{print $0, "species2"} /^L/{print $0, "species3"}' PnNmmmrv_filtered.list > PnNmmmrv_filtered_species.list
# move file to local machine
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.2.Q .
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.3.Q .
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.4.Q .
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/PnNmmmrv_filtered.cv.error .
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/PnNmmmrv_filtered_species.list .
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/log2.out .
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/log3.out .
rsync UniServer:/home/inf-51-2023/Research_project/ADMIXTURE/log4.out .


####################################################

# calculate a mean genome-wide FST value for our variant set
# make a genome scan directory
mkdir genome_scan
# make an environmental variable for the input vcf
VCF=~/Research_project/ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf
# population information in the sample names
# move into the genome scan directory
cd ~/genome_scan
# run bcftools query
bcftools query -l $VCF 
# grep all starting with 8
bcftools query -l $VCF | grep "^8" > 8N
# grep all starting with K
bcftools query -l $VCF | grep "^K" > K0
# grep all starting with L
bcftools query -l $VCF | grep "^L" > Lesina

# path to the VCF file
VCF=~/Research_project/ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf 
# calculated Weir & Cockerham’s FST for each SNP
# Population 8N and K0 that we can use for our FST analyses
vcftools --vcf ${VCF} --weir-fst-pop 8N --weir-fst-pop K0 --out ./8N_K0
# Population 8N and Lesina_280 that we can use for our FST analyses
vcftools --vcf ${VCF} --weir-fst-pop 8N --weir-fst-pop Lesina --out ./8N_Lesina
# Population K0 and Lesina that we can use for our FST analyses
vcftools --vcf ${VCF} --weir-fst-pop K0 --weir-fst-pop Lesina --out ./K0_Lesina

rsync UniServer:/home/inf-51-2023/Research_project/genome_scan/8N_K0.weir.fst .
rsync UniServer:/home/inf-51-2023/Research_project/genome_scan/8N_Lesina.weir.fst .
rsync UniServer:/home/inf-51-2023/Research_project/genome_scan/K0_Lesina.weir.fst .
