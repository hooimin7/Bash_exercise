rm(list = ls())
# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("gplots")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("magrittr")

# load the packages
library(DESeq2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
#check your wd
getwd()
list.files(recursive = T)
# Load the count matrix into a new matrix object named count_mat 
count_mat <- as.matrix(read.csv("04_Counting/final.matrix.txt", sep = "\t",
                                row.names = "gene_id"))
head(count_mat) # Check
# fh1 fh2 fh3 ref1 ref2 ref3
# CE100114_10434  17  10   5   23   10    6
# CE100464_580     0   1   0    1    2    1
# CE100520_1543    7   0   2    7    4    0
# CE100533_337     1   2   0    2    0    0
# CE100719_908     1   0   2    3    1    0
# CE101042_51      1   0   0    1    0    0

# new character vector called sampleNames, 
# specify the names of each file in the same order they appear in the count matrix
sampleNames <- c("fh1", "fh2", "fh3", "ref1", "ref2", "ref3")
sampleNames # Check
# specify the conditions for each file
# forest hot is a rich media and the control replicates are grown in poor media.
sampleConditions <- c("rich", "rich", "rich", "poor", "poor", "poor")
sampleConditions # Check
# Create a metadata table
sampleTable <- data.frame(condition = as.factor(sampleConditions))
row.names(sampleTable) <- sampleNames # Set the first column as row names
sampleTable # Check
# 6.2.1 Create the DESeq object
# DESeqDataSet object, we will need the count matrix (count_mat) and the metadata table
# (sampleTable) as input. We will also need to specify a design formula
# Design formula indicates how the column(s) in the metadata table should be used in the analysis
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = sampleTable,
                              design = ~ condition)
# Take a look:
dds
#DESeq2 import functions with help("DESeqDataSet-class")
# DESeq2 has a function to import HTSeq count files directly (DESeqDataSetFromHTSeqCount) 
# 6.3 Normalization
# systematic biases such as differences in total sequencing depth
# main factors considered during normalization are sequencing depth and RNA composition
# Median of Ratios normalization (MRN) method
# 1. A pseudo-reference sample is created for each gene by taking the row-wise geometric mean.
# 2. The ratio is calculated of each sample in relation to this pseudo-reference.
# 3. A normalization factor (aka size factor) is calculated from the median value of all ratios for each
# column.
# 4. Normalized count values are calculated using the normalization factor.
# 6.3.1 Normalize the counts within the DESeq object
# Estimate the size factors
dds <- estimateSizeFactors(dds)
# Take a look at the size factors
sizeFactors(dds) # Check the size factors
# 16. Based on the size factor calculated for fh1, were read counts in this sample systematically lower or
# higher than expected? What other characteristic of fh1 (and ref1) would hint at this answer?
# Answer: the size factor for fh1 is higher than others, which we have seen this 
# on the coverage reads (multiQc hmtl after trimmed)
# fh1       fh2       fh3      ref1      ref2      ref3 
# 2.1100351 0.7052493 0.7021197 2.1686521 0.7170036 0.7116748 
# normalized counts
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
head(normalized_counts)
# fh1       fh2      fh3       ref1      ref2     ref3
# CE100114_10434 8.0567381 14.179382 7.121293 10.6056658 13.946931 8.430817
# CE100464_580   0.0000000  1.417938 0.000000  0.4611159  2.789386 1.405136
# CE100520_1543  3.3174804  0.000000 2.848517  3.2278113  5.578772 0.000000
# CE100533_337   0.4739258  2.835876 0.000000  0.9222318  0.000000 0.000000
# CE100719_908   0.4739258  0.000000 2.848517  1.3833477  1.394693 0.000000
# CE101042_51    0.4739258  0.000000 0.000000  0.4611159  0.000000 0.000000
# alternatively
normalized_counts %>% View()
# the %>% symbolizes a pipe in R using the magrittr package, similar to | in bash
# basic statistics for each sample
summary(normalized_counts)
# fh1                fh2                fh3                ref1               ref2               ref3         
# Min.   :   0.000   Min.   :   0.000   Min.   :   0.000   Min.   :   0.000   Min.   :   0.000   Min.   :   0.000  
# 1st Qu.:   0.000   1st Qu.:   0.000   1st Qu.:   0.000   1st Qu.:   0.000   1st Qu.:   0.000   1st Qu.:   0.000  
# Median :   1.422   Median :   1.418   Median :   1.424   Median :   1.383   Median :   1.395   Median :   1.405  
# Mean   :   6.891   Mean   :   6.809   Mean   :   6.846   Mean   :   6.687   Mean   :   6.736   Mean   :   6.799  
# 3rd Qu.:   5.687   3rd Qu.:   5.672   3rd Qu.:   5.697   3rd Qu.:   5.995   3rd Qu.:   5.579   3rd Qu.:   5.620  
# Max.   :1252.112   Max.   :1213.755   Max.   :1231.984   Max.   :1138.956   Max.   :1030.678   Max.   :1131.135  
# 17. Are the median and mean read counts similar? Why or why not?
# Answer: No, the mean is higher than the median, which means the data is skewed to the right
# 6.4 Quality control and preliminary visualization
# Principal Component Analysis (PCA) and hierarchical clustering
# PCA is a technique used to emphasize variation and bring out strong patterns in the data.
# PCA:biological replicates belonging to the same experimental group to cluster together
# Hierarchical clustering also groups samples by expression similarity.
# A hierarchical tree indicates which samples are more similar to each other based on the normalized gene expression values
# A heatmap is used to show the correlation of gene expression for all pairwise combinations of samples in the data
# 6.4.1 Transformation
# unsupervised clustering methods such as PCA and hierarchical clustering are sensitive to the
# magnitude of the data ,transformation of the normalized counts improves the estimated correlations 
# or distances between samples or genes for visualization.
# Transform the data using regularized log transformation
rld <- rlog(dds, blind=TRUE) # blind=TRUE to avoid using the condition information
# if conditions are expected to be VERY different in their total counts, it is a good idea to change this to blind == FALSE.
# 6.4.2 Principal component analysis
# Create a PCA using the DESeq2 function plotPCA()
# 500 most variable genes, which is the default value for the ntop option within the plotPCA() function
# PCA plot
# make directory for figures
dir.create("05_DiffExpr")
DESeq2::plotPCA(rld, intgroup=c("condition"), ntop = 500) # intgroup is the grouping variable
# Save to a pdf using the pdf() device.
pdf("05_DiffExpr/fig1.pca.pdf") # open the pdf device
DESeq2::plotPCA(rld, intgroup=c("condition"), ntop = 500) # remember the default for ntop is 500,
# so the results would be the same if you left this out
# remember the default for ntop is 500,
# so the results would be the same if you left this out
dev.off() # close the pdf device
# 6.4.3 Hierarchical clustering
# visualize which samples cluster together when using all of the filtered, normalized and
# transformed expression data
# creating a distance matrix among samples using stats::dist().
distsRL <- dist(t(assay(rld))) # distsRL is a distance matrix
mat <- as.matrix(distsRL) # convert to matrix
# rename the columns and rows of the matrix with the column names from the dds
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, sampleNames, sep=" : "))
# clustering analysis on the samples using the function stats::hclust()
hc <- hclust(distsRL)
# Visualize the distance between samples using a heatmap
# show the clustering of samples using a dendrogram
install.packages("ComplexHeatmap") # no available for my R version
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # create a color palette
# Plot the heatmap
heatmap.2(mat,
          Rowv=as.dendrogram(hc),
          symm=TRUE,
          trace="none",
          col = rev(hmcol),
          margin=c(13, 13))
# save the heatmap with the pdf() device
pdf("05_DiffExpr/fig2.heatmap.pdf")
heatmap.2(mat,
          Rowv=as.dendrogram(hc),
          symm=TRUE,
          trace="none",
          col = rev(hmcol),
          margin=c(13, 13))
dev.off()
# 18. Look at your PCA and heatmap. Do the samples cluster as you would expect from the experimental
# design? Are there any major outliers or concerning clusters
# Answer: The PCA plot shows that the samples cluster by condition, which is what we would expect.
# The heatmap also shows that the samples cluster by condition.
# 6.5 Differential expression
# we have checked the quality of our samples to see if they fit the expectation of our experimental design
# have checked our samples for outliers and explored patterns using dimension reduction and clustering
# conduct differential expression analysis
# run the actual analysis, we use the function DESeq().
# It performs the normalization of the counts using the size factors, estimates dispersion, and fits the linear
# model
# Use the assign arrow(<-) to replace the old dds object with the newly fitted object
# Run the analysis
# re-assigning the results of the function back to the same variable name (dds), we can fill in the slots of
# our DESeqDataSet object with the new values calculated by DESeq2
dds <- DESeq(dds)
# using pre-existing size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# final dispersion estimates
# fitting model and testing
# 6.5.1 Identifying DE genes
# With DESeq2, the Wald test is commonly used for hypothesis testing when comparing two groups
# To indicate to DESEq2 the two groups we want to compare, we provide contrasts, or 
# the specific pair of groups we want to test
# we have to check what contrasts are possible given our data and design formula
# check what contrasts are even possible
DESeq2::resultsNames(dds)
# [1] "Intercept"              "condition_rich_vs_poor"
# Extract the results for Poor vs. Rich
# specify the contrast by naming the factor you want to test and
# the levels of the factor you want to compare.
contrast_pr <- c("condition", "poor", "rich")
# extract the results for your specified contrast
res_table <- results(dds, contrast=contrast_pr)
# res_table_pm <- results(dds, contrast=contrast_pm)
# alternatively you could use:
# res_table <- results(dds, name = "condition_rich_vs_poor")
# Wald test statistic is computed per-gene, and evaluates the probability that a test statistic at least as
# extreme as the observed value would occur at at random
# We reject the null hypothesis (i.e. there is no difference in expression between groups) and conclude 
# a gene is differentially expressed when the calculated p-value is below our significance threshold 
# (usually 0.05).
# sort results based on significance
res_table <- res_table[order(res_table$padj),]
# You will now have the smallest padj values at the top of your table
head(res_table)
# log2 fold change (MLE): condition poor vs rich 
# Wald test p-value: condition poor vs rich 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   fgenesh1_kg.12__183__Locus432v1rpkm387.16   588.740        3.86140 0.1208459   31.9531 4.88813e-224 3.28629e-220
# estExt_Genemark1.C_130215                   386.051        5.25702 0.1950922   26.9463 6.30090e-160 2.11805e-156
# fgenesh1_kg.7__381__Locus862v1rpkm207.81    222.891        6.19791 0.3318836   18.6750  7.91533e-78  1.77383e-74
# fgenesh1_kg.23__105__Locus2331v1rpkm69.56   143.510       -3.02321 0.1832531  -16.4975  3.82638e-61  6.43119e-58
# fgenesh1_kg.21__184__Locus30v1rpkm2293.91   609.877       -1.45415 0.0890015  -16.3385  5.25582e-60  7.06697e-57
# e_gw1.29.217.1                              147.326        2.23810 0.1605483   13.9404  3.60119e-44  4.03513e-41
# The column headers are:
#   • baseMean: mean of normalized counts for all samples
# • log2FoldChange: log2 fold change
# • lfcSE: standard error of log2 fold change
# • stat: the Wald test statistic
# • pval: Wald test p-value
# • padj: Benjamini-Hochberg (BH) adjusted p-values
# write your results into a separate file and save them to your computer
# first convert the rownames to the first column of the dataframe so that
# when we write it out to a table, we have our gene names as well.
gene_id = rownames(res_table) # extract the gene names
res_table <- cbind(gene_id , data.frame(res_table, row.names=NULL)) # add the gene names as a column
# write the result table in a file.
write.table(res_table , file = "05_DiffExpr/diffExpr.tab", sep = "\t", quote = FALSE, 
            row.names = FALSE)
# Adjusting the p-values to account for testing multiple genes helps us to have more confidence in our conclusions.
# can use log2 fold change (LFC), to identify DE genes
# common threshold is an absolute LFC of 1, which means expression of that gene is twice as much in one of the conditions compared
# to the other
# If you just want a list of the genes that have an adjusted p-value cutoff below a certain threshold and above
# an absolute log fold change threshold, we can easily subset the results table. Here we will use thresholds of
# padj < 0.05 and log2FoldChange < -1 or log2FoldChange > 1.
resSig <- subset(res_table , padj < 0.05 & abs(log2FoldChange) > 1) # subset the results table
write.table(resSig, file = "05_DiffExpr/diffExpr.sig05.tab", sep = "\t", quote = FALSE,
            row.names = FALSE)
# significantly differentially expressed genes had a negative log2FoldChange 
# and how many had a positive log2FoldChange.
# call genes with a positive log2FoldChange as “upregulated” and genes with 
# a negative log2FoldChange as downregulated
# identify genes that are "upregulated" in group 1 relative to group 2
resSig_upreg <- rownames(resSig)[resSig$log2FoldChange > 0]
length(resSig_upreg)
# [1] 116
# identify genes that are "downregulated" in group 1 relative to group 2
resSig_downreg <- rownames(resSig)[resSig$log2FoldChange < 0]
length(resSig_downreg)
# [1] 156
# What does a positive or negative log2FoldChange actually mean for this comparison?
# Positive log2FoldChange: Upregulated in "poor" compared to "rich."
# Negative log2FoldChange: Downregulated in "poor" compared to "rich."
# The log2FoldChange is the log2 of the fold change, which is the ratio of the normalized counts in the two groups.
# A positive log2FoldChange means that the normalized counts in the first group are higher than the normalized counts in the second group.
# A negative log2FoldChange means that the normalized counts in the second group are higher than the normalized counts in the first group.
# 6.5.2 Visualizing DE Genes
# 6.5.2.1 Individual genes
# the gene that is most differentially expressed (i.e., lowest padj) in the sorted results sort using head().
head(resSig, 1) # 1 = the first row
tail(resSig, 1) # 1 = the last row)

# gene_id baseMean log2FoldChange     lfcSE     stat        pvalue          padj
# 1 fgenesh1_kg.12__183__Locus432v1rpkm387.16 588.7405       3.861403 0.1208459 31.95313 4.888133e-224 3.286292e-220
# Is it up or down-regulated in P. involutus on poor nutrient media relative to rich nutrient media?
plotCounts(dds, gene="fgenesh1_kg.12__183__Locus432v1rpkm387.16", intgroup="condition")
# why the plot keep changing?
# the x-axis keep changing due to the different number of samples in each group (plotCounts), but the y-axis is the same
# The plotCounts() function plots the normalized counts for each sample for a given gene.
# 6.5.3 All genes: MA plot
# MA plot shows the mean of the normalized counts (x-axis) versus the LFC (y-axis) for all genes.
# genes that are significantly differentially expressed are colored to be easily identified
# The default significance threshold is 0.1, but today let’s use p < 0.05
# positive LFC values mean the genes are more expressed in poor nutrients 
# MA plot needs the dds object
plotMA(dds, ylim=c(-2,2), main="DESeq2", alpha = 0.05)
# save the MA plot using the pdf() device
pdf("05_DiffExpr/fig3.MA.pdf")
plotMA(dds, ylim=c(-2,2), main="DESeq2", alpha = 0.05)
dev.off()
# This plot allows us to evaluate the magnitude of fold changes and how they are distributed relative to mean expression
# When within-group variation is high, LFC values can be a little misleading. Just look at all the nonsignificant
# genes with very high or low LFC in our MA and volcano plots. In this example data set, you have lots of
# genes with very small count numbers (between 0 to 10) that also have large log2fold changes (>1). Using
# DESeq2, there is a way to shrink the log2 fold change for genes with low counts so that you can have more
# confidence in the fold change estimates. You can have a look on documentation and read about log2 fold
# change shrinkage.
# 6.5.4 Subset of genes: heatmap
# Make a heatmap of the 50 most differentially expressed (lowest p-value) genes:
#   • subset the resuls sort to only include significantly differentially expressed genes using an absolute log2
# fold change cutoff of 1 (abs(log2FoldChange) > 1) and an adjusted p-value (FDR=FalseDiscoveryRate) cutoff of 0.01.
# • create a list of the 50 most differentially expressed gene IDs from your results table
# • subset the VST transformed expression matrix to include only the rows for these 50 genes.
# • make a heatmap of this subset using ‘pheatmap.
# • save your result using the pdf() device
# Identify the top 50 differentially expressed genes from the sorted results table.
top50 <-res_table$gene_id[res_table$padj < 0.05 & abs(res_table$log2FoldChange) > 1][1:50]
# Extract the expression levels of the most expressed genes from the vst-transformed dds object
counts_top50 <- counts(dds, normalized=TRUE)[top50,] # subset the counts matrix to include only the rows for the top 50 genes
hmcol <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100) # create a color palette
heatmap.2(counts_top50, 
          col = hmcol ,
          Rowv = FALSE ,
          Colv = FALSE,
          scale = "none",
          dendrogram = "none",
          trace = "none")
# Save with pdf()
pdf("05_DiffExpr/fig4.heatmap_sig_diff_exp.pdf")
heatmap.2(counts_top50,
          col = hmcol ,
          Rowv = FALSE ,
          Colv = FALSE,
          scale = "none",
          dendrogram = "none",
          trace = "none",
          margin = c(10,6))
dev.off()
# 7 Additional Exercises
# 19. How many genes show a significant difference in expression on the FDR level of 0.01? Compare with
# the above result
# Answer: 272
resSig <- subset(res_table , padj < 0.01 & abs(log2FoldChange) > 1) # abs(log2FoldChange) > 1 means the absolute value of log2FoldChange is greater than 1
# How many genes show a significant difference 
print(resSig)
write.table(resSig, file = "05_DiffExpr/diffExpr.sig05.tab", sep = "\t", quote = FALSE,
            row.names = FALSE)
# 20. Produce a heatmap of the 30 most variable genes, not the 30 most differentially expressed genes
# Answer: 

