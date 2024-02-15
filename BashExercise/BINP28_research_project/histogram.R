rm(list = ls())
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

