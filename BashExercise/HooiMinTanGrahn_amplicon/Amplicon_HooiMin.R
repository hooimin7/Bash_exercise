rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyverse)

# read the count matrix

Bacterial <- read.csv("new_file1", sep = "\t", header = FALSE)
Otus <- read.csv("otus.tsv", sep = "\t")

row.names(Otus) <- Otus[,1] # set the row names to the first column
row.names(Bacterial) <- Bacterial$V1 # set the row names to the first column

Otus$Phylum <- Bacterial[rownames(Otus),2] # add the phylum column to the otus 
# data frame
Otus <- na.omit(Otus) # remove the rows with NA values
Otus <- as.data.frame(Otus) # convert the otus data set to a data frame
# make aggregate data frame
otu_agg <- aggregate(. ~ Phylum, data = Otus[,-1], sum) # summing the values of 
# all columns (except the first one) in the `Otus` data frame, for each unique 
# value of `Phylum`

dim(otu_agg) # check the dimensions of the data frame

otu_agg$Phylum <- factor(otu_agg$Phylum, levels = otu_agg$Phylum) # convert the 
# Phylum column to a factor

sample_names <- colnames(otu_agg)[-1] # get the sample names

filtered <- otu_agg[rowSums(otu_agg[,sample_names]) > 6000,] # filter the data


# long format
otu_long <- pivot_longer(data = filtered, cols = -Phylum, values_to = "abundance", 
                         names_to = "Sample")
dim(otu_long)

# plot
ggplot(data = otu_long, aes(x = Sample, y = abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  labs(title = "Phylum abundance", x = "Phylum", y = "Abundance") +
  scale_fill_brewer(palette = "Set3")

