#!/usr/bin/env Rscript

# load appropriate libraries
library(ggplot2)
library(dplyr)

# get args args[1] = telo_stats (file), args[2] = plot_telo_length (boolean), args[3] = plot_vrr_length (boolean), args[4] = strand_comparison (boolean)
args = commandArgs(trailingOnly=TRUE)
# open telo stats file into dataframe
telo_stats <- read.table(args[1], header=TRUE)

plt <- ggplot(data=telo_stats) +
  geom_histogram(mapping=aes(x=telo_length), binwidth=10) +
  theme_minimal() + xlab("Telomere Length (bp)") + ylab("Number of Telomeres") +
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=9))

ggsave("telomere_length.pdf", width=12, height=8, units="cm", plot=plt)


plt <- ggplot(data=telo_stats) +
  geom_histogram(mapping=aes(x=telo_length, color=strand), binwidth=10) +
  theme_minimal() + xlab("Telomere Length (bp)") + ylab("Number of Telomeres") +
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=9))

ggsave("telomere_length_by_strand.pdf", width=12, height=8, units="cm", plot=plt)

print(paste("Number of Telomeres: ", length(telo_stats$telo_length)))
print(paste("Mean Telomere Length: ", mean(telo_stats$telo_length)))
print(paste("Median Telomere Length: ", median(telo_stats$telo_length)))