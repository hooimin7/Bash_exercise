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

################## in rstudio (creates histogram of mean depth)
# Load necessary library
library(ggplot2)

# Data
data <- read.csv("out.ldepth.mean", sep = "\t", header = TRUE)

# Create histograms
mean_depth_hist <- ggplot(data, aes(x = MEAN_DEPTH)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Mean Depth", x = "Mean Depth", y = "Frequency") +
  theme_minimal() + xlim(0, 50)


# Display histograms
mean_depth_hist
#####################

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

####################### in rstudio (create barplot for SNP density)
# Load the ggplot2 package
library(ggplot2)

# Read the data into a data frame
df <- read.csv("out.snpden_chr5", sep = "\t", header = TRUE)

# Ensure the BIN_START and SNP_COUNT columns are numeric
df$BIN_START <- as.numeric(df$BIN_START)
df$SNP_COUNT <- as.numeric(df$SNP_COUNT)

# Create a bar plot of SNP count
ggplot(df, aes(x = BIN_START, y = SNP_COUNT)) +
  geom_bar(stat = "identity") +
  labs(x = "Bin Start Position", y = "SNP Count", title = "chr5 SNP Density")

# Read the data into a data frame
df <- read.csv("out.snpden_chrZ", sep = "\t", header = TRUE)

# Ensure the BIN_START and SNP_COUNT columns are numeric
df$BIN_START <- as.numeric(df$BIN_START)
df$SNP_COUNT <- as.numeric(df$SNP_COUNT)

# Create a bar plot of SNP count
ggplot(df, aes(x = BIN_START, y = SNP_COUNT)) +
  geom_bar(stat = "identity") +
  labs(x = "Bin Start Position", y = "SNP Count", title = "chrZ SNP Density")

###########################

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

################################ (in rstudio create PCA)
# load tidyverse package
library(tidyverse)

# read in data
pca <- read_table2("cichlids.eigenvec", col_names = FALSE) # read in eigenvec
eigenval <- scan("cichlids.eigenval") # read in eigenvalues

# sort out the pca data
# remove nuisance column
pca <- pca[,-1] # remove the first column

# set names
names(pca)[1] <- "ind" # set the first column name to "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1)) # set the rest of the
# column names to PC1, PC2, etc.

# sort out the individual species and pops
# species
spp <- rep(NA, length(pca$ind))
spp[grep("^8", pca$ind)] <- "8N"
spp[grep("^K", pca$ind)] <- "K0"
spp[grep("^L", pca$ind)] <- "Lesina"

spp_loc <- spp # copy spp to spp_loc
# remake data.frame

pca <- as_tibble(data.frame(pca, spp, spp_loc))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100) 

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light() + xlab("Principal component") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_text(aes(label = round(pve, 2)), vjust = -0.5, size = 3) + 
  theme(plot.title = element_text(hjust = 0.5))

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp)) + geom_point(size = 3) 
b <- b + scale_colour_manual(values = c("yellow", "blue", "red"), name = "Taxa") # set the colours
b <- b + coord_equal() + theme_light() # set the aspect ratio

# add labels
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) 
b <- b + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b
####################################

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

######################################### (in rstudio create Admixture proportions plot)

# Load necessary library
library(ggplot2)
library(dplyr)
library(tidyverse)

samplelist <- read_tsv("PnNmmmrv_filtered_species.list", 
                       col_names = c("sample", "species"),
                       show_col_types = FALSE) 

read_delim("ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.2.Q",
           col_names = paste0("Q",seq(1:2)),
           delim=" ") 

# create an empty tibble to store the data
all_data <- tibble(sample=character(), 
                   k=numeric(),
                   Q=character(),
                   value=numeric()) 

for (k in 2:4){ # loop through the Q files
  data <- read_delim(paste0("ProjTaxa_no_Naxos2_minDP_maxDP_maxMiss.recode.vcf.filtered.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)), # create column names
                     delim=" ") # specify the delimiter
  data$sample <- samplelist$sample # add sample column
  data$k <- k # add k column
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data # print the data


# plot the data for k = 2
all_data %>%
  filter(k == 2) %>% # filter for k = 2
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  ggplot(aes(x=sample, y=value, fill=Q)) + 
  geom_bar(stat="identity", position="stack") +
  scale_x_discrete(name = "Sample") +
  ylab("Admixture proportions") +
  labs(fill = "K", subtitle = "K = 2") + # change legend title to "K"
  scale_fill_discrete(labels = c("1", "2")) + # change legend labels to "1", "2"
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))


# plot the data for k = 2, 3, 4
all_data %>% 
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Admixture proportions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=c("1", "2", "3", "4")) + # change legend labels to "1", "2", "3", "4"
  facet_wrap(~k,ncol=1)

##########################################

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

############################### (in rstudio create FST manhattan plot)
8N vs K0
######
# Load necessary library
library(ggplot2)
library(dplyr)
library(tidyverse)
library(qqman)

# read in the file
fst <- read_table("8N_K0.weir.fst")

# add snp column
# this data does not include snp name, I will add a column with snp name
# assigning a number to each snp
length_ = dim(fst)[1] # number of rows
fst$SNP = paste("SNP", 1:length_) # add snp name
head(fst)

# remove 'chr' prefix and convert to numeric
fst$CHROM <- as.numeric(gsub("chr", "", fst$CHROM))

# replace NA values with 1
fst$CHROM[is.na(fst$CHROM)] <- 1

# replace -nan with NA
fst$WEIR_AND_COCKERHAM_FST <- gsub("-nan", NA, fst$WEIR_AND_COCKERHAM_FST)

# convert to numeric
fst$WEIR_AND_COCKERHAM_FST <- as.numeric(fst$WEIR_AND_COCKERHAM_FST)

# remove rows with NA in WEIR_AND_COCKERHAM_FST
fst <- fst[!is.na(fst$WEIR_AND_COCKERHAM_FST),]

# Adjust margins
par(mar = c(5, 4, 4, 2) + 0.1)

# generate a manhattan plot
manhattan(fst, chr = 'CHROM', bp = 'POS', snp = 'SNP', p = 'WEIR_AND_COCKERHAM_FST', 
          col = c("blue", "red"), logp = FALSE, ylab = 'WEIR AND COCKERHAM FST',
          xlab = 'CHR')

# add note to the plot
mtext("Note: 8N vs K0", side = 1, line = 4)

# Add legend outside of plot
legend(x = 0.5, y = -0.1, legend = c("1 = chrZ", "5 = chr5"), fill = c("blue", "red"), xpd = TRUE)

############
8N vs Lesina
##########
# Load necessary library
library(ggplot2)
library(dplyr)
library(tidyverse)
library(qqman)

# read in the file
fst <- read_table("8N_Lesina.weir.fst")

# add snp column
# this data does not include snp name, I will add a column with snp name
# assigning a number to each snp
length_ = dim(fst)[1] # number of rows
fst$SNP = paste("SNP", 1:length_) # add snp name
head(fst)

# remove 'chr' prefix and convert to numeric
fst$CHROM <- as.numeric(gsub("chr", "", fst$CHROM))

# replace NA values with 1
fst$CHROM[is.na(fst$CHROM)] <- 1

# replace -nan with NA
fst$WEIR_AND_COCKERHAM_FST <- gsub("-nan", NA, fst$WEIR_AND_COCKERHAM_FST)

# convert to numeric
fst$WEIR_AND_COCKERHAM_FST <- as.numeric(fst$WEIR_AND_COCKERHAM_FST)

# remove rows with NA in WEIR_AND_COCKERHAM_FST
fst <- fst[!is.na(fst$WEIR_AND_COCKERHAM_FST),]

# Adjust margins
par(mar = c(5, 4, 4, 2) + 0.1)

# generate a manhattan plot
manhattan(fst, chr = 'CHROM', bp = 'POS', snp = 'SNP', p = 'WEIR_AND_COCKERHAM_FST', 
          col = c("blue", "red"), logp = FALSE, ylab = 'WEIR AND COCKERHAM FST',
          xlab = 'CHR')

# add note to the plot
mtext("Note: 8N vs Lesina", side = 1, line = 4)

# Add legend outside of plot
legend(x = 0.5, y = -0.1, legend = c("1 = chrZ", "5 = chr5"), fill = c("blue", "red"), xpd = TRUE)
####################

K0 vs Lesina
############
# Load necessary library
library(ggplot2)
library(dplyr)
library(tidyverse)
library(qqman)

# read in the file
fst <- read_table("K0_Lesina.weir.fst")

# add snp column
# this data does not include snp name, I will add a column with snp name
# assigning a number to each snp
length_ = dim(fst)[1] # number of rows
fst$SNP = paste("SNP", 1:length_) # add snp name
head(fst)

# remove 'chr' prefix and convert to numeric
fst$CHROM <- as.numeric(gsub("chr", "", fst$CHROM))

# replace NA values with 1
fst$CHROM[is.na(fst$CHROM)] <- 1

# replace -nan with NA
fst$WEIR_AND_COCKERHAM_FST <- gsub("-nan", NA, fst$WEIR_AND_COCKERHAM_FST)

# convert to numeric
fst$WEIR_AND_COCKERHAM_FST <- as.numeric(fst$WEIR_AND_COCKERHAM_FST)

# remove rows with NA in WEIR_AND_COCKERHAM_FST
fst <- fst[!is.na(fst$WEIR_AND_COCKERHAM_FST),]

# Adjust margins
par(mar = c(5, 4, 4, 2) + 0.1)

# generate a manhattan plot
manhattan(fst, chr = 'CHROM', bp = 'POS', snp = 'SNP', p = 'WEIR_AND_COCKERHAM_FST', 
          col = c("blue", "red"), logp = FALSE, ylab = 'WEIR AND COCKERHAM FST',
          xlab = 'CHR')

# add note to the plot
mtext("Note: K0 vs Lesina", side = 1, line = 4)

# Add legend outside of plot
legend(x = 0.5, y = -0.1, legend = c("1 = chrZ", "5 = chr5"), fill = c("blue", "red"), xpd = TRUE)
#######################
