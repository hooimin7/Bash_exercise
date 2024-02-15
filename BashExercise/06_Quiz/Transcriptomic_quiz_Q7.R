#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) # Get user input

# myScript.R infile_forMyScript.txt outfile_forMyScript.txt
infile <- args[1]
outfile <- args[2]

# load the suggested packages
library(DESeq2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(magrittr)

# Get user input for the host plant to compare with Ribes
other_host <- readline(prompt="Enter the other host plant to compare with Ribes: ")

# Define the sample names based on user input
sampleNames <- c(other_host_id, Ribes_id)

# Define the sample conditions based on user input
sampleConditions <- c(other_host_id, Ribes_id)

# Create a sample table from the sample conditions
sampleTable <- data.frame(condition = as.factor(sampleConditions))
# Set the first column as row names
row.names(sampleTable) <- sampleNames 
sampleTable # Check
# compare matrix column names with sample names
colnames(countMatrix) %in% rownames(sampleTable)

# Specify the design veriable as condition
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = sampleTable,
                              design = ~ condition)

dds <- estimateSizeFactors(dds) # Estimate the size factors
sizeFactors(dds) # Check the size factors
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE)) # N
# Define the contrast based on user input

rld <- rlog(dds, blind=TRUE) # Transform the data using rlog
DESeq2::plotPCA(rld, intgroup=c("condition"), ntop = "USERINPUT") # Plot PCA

dds <- DESeq(dds) # Run DESeq
DESeq2::resultsNames(dds) # Check the results names
contrast <- c("condition", other_host, "Ribes") 

res_table <- results(dds, contrast=contrast) # Extract the results
res_table <- res_table[order(res_table$padj),] # Order the results by padj

# Subset the results for genes that are significant differentiated on Ribes and the other host
resSig <- subset(res_table , padj < 0.01 & log2FoldChange >= 1) 

# identify genes that are "upregulated" in Ribes and the other host
resSig_upreg <- rownames(resSig)[resSig$log2FoldChange > 0]
length(resSig_upreg)

# identify genes that are "downregulated" in Ribes and the other host
resSig_downreg <- rownames(resSig)[resSig$log2FoldChange < 0]
length(resSig_downreg)

top50 <-res_table$other_host_id[res_table$padj < 0.05 & abs(res_table$log2FoldChange) > 1][1:50] # Get the top 50 genes

# Print the total number of genes that are upregulated on Ribes and the other host
print(paste("Total number of genes upregulated on Ribes and", other_host, ":", length(resSig_upreg)))

# Write the results to a file
write.table(top50, file = paste("diffExpr_", other_host, "_Ribes.tab", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE) 