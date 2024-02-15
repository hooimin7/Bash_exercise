rm(list = ls())
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
