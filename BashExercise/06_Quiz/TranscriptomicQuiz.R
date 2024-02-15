# rm(list = ls())
# load the packages
library(DESeq2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
getwd()
# read the count matrix
countMatrix <- read.csv("transcriptomicsQuizData.txt",
                        sep = "\t", row.names = "gene_id") %>% 
  as.matrix()
# check the column names of the count matrix
colnames(countMatrix)
# define the sample names
sampleNames <- c("Pca_Ribes_rep1", "Pca_Ribes_rep2", "Pca_Ribes_rep3",
                 "Pca_Ribes_rep4", "Pca_Ribes_rep5", "Pca_Ribes_rep6",
                 "Pca_Ribes_rep7", "Pca_Ribes_rep8", "Pca_Ribes_rep9",
                 "Pca_Salix_rep1", "Pca_Salix_rep2", "Pca_Salix_rep3",
                 "Pca_Salix_rep4", "Pca_Salix_rep5", "Pca_Salix_rep6",
                 "Pca_Salix_rep7", "Pca_Salix_rep8", "Pca_Salix_rep9",
                 "Pca_Urtica_rep1", "Pca_Urtica_rep2", "Pca_Urtica_rep3",
                 "Pca_Urtica_rep4", "Pca_Urtica_rep5", "Pca_Urtica_rep6",
                 "Pca_Urtica_rep7", "Pca_Urtica_rep8", "Pca_Urtica_rep9")
sampleNames # Check
# Define the sample conditions
sampleConditions <- c("Ribes", "Ribes", "Ribes",
                      "Ribes", "Ribes", "Ribes",
                      "Ribes", "Ribes", "Ribes",
                      "Salix", "Salix", "Salix",
                      "Salix", "Salix", "Salix",
                      "Salix", "Salix", "Salix",
                      "Urtica", "Urtica", "Urtica",
                      "Urtica", "Urtica", "Urtica",
                      "Urtica", "Urtica", "Urtica")
sampleConditions # check
# Create a sample table from the sample conditions
sampleTable <- data.frame(condition = as.factor(sampleConditions))
row.names(sampleTable) <- sampleNames # Set the first column as row names
sampleTable # Check
# compare matrix column names with sample names
colnames(countMatrix) %in% rownames(sampleTable)
# Create the deseq object using the count matrix and sample table.
# Specify the design veriable as condition
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = sampleTable,
                              design = ~ condition)
dds #Check
dds <- estimateSizeFactors(dds) # Estimate the size factors
sizeFactors(dds) # Check the size factors
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE)) # Normalize the counts
head(normalized_counts) # Check
summary(normalized_counts) # Check
rld <- rlog(dds, blind=TRUE) # Transform the data using rlog
DESeq2::plotPCA(rld, intgroup=c("condition"), ntop = 1000)
# Using a regularized log transformation, make a PCA plot of the top 1000 most variable genes. Do the
# samples cluster more by condition on the first or second PC axes, both or neither
# Answer: PC1
pdf("fig1.pca.pdf") # open the pdf device
DESeq2::plotPCA(rld, intgroup=c("condition"), ntop = 1000)
dev.off() 
# 3. Run DESeq on the dds. Extract the results for Urtica vs. Ribes. Identify the most differentially
# expressed gene using the adjusted P-value (padj)
dds <- DESeq(dds) # Run DESeq
DESeq2::resultsNames(dds) # Check the results names
# [1] "Intercept"                 "condition_Salix_vs_Ribes"  "condition_Urtica_vs_Ribes"
contrast_UR <- c("condition", "Urtica", "Ribes") # Define the contrast
res_table <- results(dds, contrast=contrast_UR) # Extract the results
res_table <- res_table[order(res_table$padj),] # Order the results by padj
head(res_table) # Check
#                baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   Polcal_g14648   151.939       -7.88699  0.641400  -12.2965 9.45626e-35 9.72482e-31
# 4. What is the expression of the gene you identified in Question 3 in Ribes relative to Urtica (rounded
# to 3 significant digits)?
# Answer: -7.88699 downregulated, 
result <- 2^7.88699
print(result)
# 5. Using an adjusted p-value threshold of 0.01 and an absolute log2FoldChange threshold of 1, how many
# genes were differentially expressed between larvae fed on Ribes and larvae fed on Urtica?
resSig <- subset(res_table , padj < 0.01 & abs(log2FoldChange) >= 1) # Subset the results
write.table(resSig, file = "diffExpr.sig06.tab", sep = "\t", quote = FALSE,
            row.names = FALSE) # Write the results to a file
resSig <- read.table("diffExpr.sig06.tab", header = TRUE, sep = "\t") # Read the results
# Print the number of significantly differentially expressed genes
print(nrow(resSig))
