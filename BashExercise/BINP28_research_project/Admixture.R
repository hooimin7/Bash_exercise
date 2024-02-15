rm(list = ls())
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

