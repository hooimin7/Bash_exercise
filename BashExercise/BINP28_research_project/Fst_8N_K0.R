rm(list = ls())
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

