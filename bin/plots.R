#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)


# script funtions to generate summary stats of entire sequencing run/experiment
args = commandArgs(trailingOnly=TRUE)
data_table <- read.table(args[1], sep="\t", header=TRUE)

if (args[2] == "individual"){
      telo_length_histogram <- ggplot(data=data_table) +
      geom_histogram(mapping=aes(x=telo_length), fill="#D81B60", binwidth=10) +
      theme_minimal() +
      xlab("Telomere Length (bp)") + ylab("Number of Telomeres") +
      theme(axis.title = element_text(size=20),
            axis.text=element_text(size=15))
      ggsave("telomere_length_histogram.pdf", device="pdf", width=12, height=10, plot=telo_length_histogram)

      stranded_telo_histo <- ggplot(data=data_table) +
      geom_histogram(mapping=aes(x=telo_length, fill=strand), binwidth=10) +
      facet_wrap(~strand, ncol=1, scales="free_y") +
      theme_minimal() +
      theme(strip.background=element_blank(),
            strip.text= element_blank(),
            axis.text=element_text(size=15),
            axis.title= element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20)) +
      xlab("Telomere Length (bp)") + ylab("Number of Telomers") +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
      ggsave("telomere_length_by_strand.pdf", device=pdf, width=12, height=10, plot=stranded_telo_histo)

      stranded_telo_box <- ggplot(data=data_table) +
      geom_boxplot(mapping=aes(y=telo_length, fill=strand)) +
      theme_minimal() +
      theme(axis.text=element_text(size=15),
            axis.title= element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20),
            axis.text.x=element_blank()) +
      ylab("Telomere Length (bp)") +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
      ggsave("telomere_length_by_strand_boxplot.pdf", device=pdf, width=6, height=10, plot=stranded_telo_box)


      cluster_telo_hist <- ggplot(data=data_table) +
      geom_histogram(mapping=aes(x=telo_length, fill=Cluster), binwidth=10) +
      facet_wrap(~Cluster, ncol=1, scales="free_y") +
      theme_minimal() +
      theme(strip.background=element_blank(),
            strip.text= element_blank(),
            axis.text=element_text(size=15),
            axis.title= element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20),
            axis.text.y=element_text(size=10)) +
      xlab("Telomere Length (bp)") + ylab("Number of Telomers")
      cluster_telo_hist
      ggsave("telomere_length_by_cluster.pdf", device=pdf, width=12, height=10, plot=cluster_telo_hist)


      cluster_telo_box <- ggplot(data=data_table) +
      geom_boxplot(mapping=aes(y=telo_length, fill=Cluster)) +
      theme_minimal() +
      theme(strip.background=element_blank(),
            strip.text= element_blank(),
            axis.text=element_text(size=15),
            axis.title= element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20),
            axis.text.x = element_blank()) +
      ylab("Telomere Length (bp)")
      ggsave("telomere_length_by_cluster_boxplot.pdf", device=pdf, width=12, height=10, plot=cluster_telo_box)


      stranded_telo_quality <- ggplot(data=data_table) +
      geom_boxplot(mapping=aes(y=telo_quality, fill=strand)) +
      theme_minimal() +
      theme(axis.text=element_text(size=15),
            axis.title= element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20),
            axis.text.x=element_blank()) +
      ylab("Telomere Quality") +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
      stranded_telo_quality
      ggsave("telomere_quality.pdf", device=pdf, width=6, height=10, plot=stranded_telo_quality)

      stranded_read_quality <- ggplot(data=data_table) +
      geom_boxplot(mapping=aes(y=read_quality, fill=strand)) +
      theme_minimal() +
      theme(axis.text=element_text(size=15),
            axis.title= element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20),
            axis.text.x=element_blank()) +
      ylab("Read Quality") +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
      ggsave("read_quality.pdf", device=pdf, width=6, height=10, plot=stranded_read_quality)

      read_len_vs_telo_len <- ggplot(data=data_table) +
      geom_point(mapping=aes(x=read_length, y=telo_length, color=strand)) +
      theme_minimal() +
      xlab("Read Length (bp)") + ylab("Telomere Length (bp)") +
      theme(axis.title = element_text(size=20),
            axis.text=element_text(size=15)) +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")

      ggsave("read_length_by_telo_length.pdf", device="pdf", width=12, height=10, plot=read_len_vs_telo_len)

      read_qual_vs_telo_qual <- ggplot(data=data_table) +
      geom_point(mapping=aes(x=read_quality, y=telo_quality, color=strand)) +
      theme_minimal() +
      xlab("Read Quality") + ylab("Telomere Quality") +
      theme(axis.title = element_text(size=20),
            axis.text=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20)) +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
      ggsave("read_quality_by_telo_quality.pdf", device="pdf", width=12, height=10, plot=read_qual_vs_telo_qual)

      read_length_vs_qual <- ggplot(data=data_table) +
      geom_point(mapping=aes(x=read_length, y=read_quality, color=strand)) +
      theme_minimal() +
      xlab("Read Length") + ylab("Read Quality") +
      theme(axis.title = element_text(size=20),
            axis.text=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20)) +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
      ggsave("read_quality_by_length.pdf", device="pdf", width=12, height=10, plot=read_length_vs_qual)

      telo_length_vs_qual <- ggplot(data=data_table) +
      geom_point(mapping=aes(x=telo_length, y=telo_quality, color=strand)) +
      theme_minimal() +
      xlab("Telomere Length (bp)") + ylab("Telomere Quality") +
      theme(axis.title = element_text(size=20),
            axis.text=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20)) +
      scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
      ggsave("telo_quality_by_length.pdf", device="pdf", width=12, height=10, plot=telo_length_vs_qual)
      telo_length_vs_qual


      telo_lengths_for_binning <- c(525, 475, 425, 375, 325, 275, 225, 175, 125, 75, 25)

      # 0-50, 50-100, 150-200, 200-250, 250-300, 30-350, 400-450, 450-500, 500-550, 550-600, 600+
      data_table$bin_telo_length <- unlist(lapply(data_table$telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
      data_table$bin_telo_length <- factor(data_table$bin_telo_length, levels=c(525, 475, 425, 375, 325, 275, 225, 200, 175, 125, 75, 25))

      telo_bar_hist <- ggplot(data=data_table) +
      geom_bar(mapping=aes(x=1, fill=bin_telo_length), position="fill") +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text = element_text(size=15),
            axis.title=element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20)) +
      ylab("Proportion of Telomere Reads") +
      guides(fill=guide_legend(title="Telomere Length (BP)")) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(breaks=c(525, 475, 425, 375, 325, 275, 225, 175, 125, 75, 25),
                        labels=c("500+", "450-499","400-449", "350-399", "300-349", "250-299", "200-249", "150-199", "100-149", "50-99", "0-49"),
                        values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))
      ggsave("telo_length_bar_plot.pdf", plot = telo_bar_hist, device="pdf", width=6, height=10)

      telo_bar_hist_by_strand <- ggplot(data=data_table) +
      geom_bar(mapping=aes(x=strand, fill=bin_telo_length), position="fill") +
      theme_minimal() +
      theme(axis.text = element_text(size=15),
            axis.title=element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20)) +
      ylab("Proportion of Telomere Reads") + xlab("Strand") +
      guides(fill=guide_legend(title="Telomere Length (BP)")) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(breaks=c(525, 475, 425, 375, 325, 275, 225, 175, 125, 75, 25),
                        labels=c("500+", "450-499","400-449", "350-399", "300-349", "250-299", "200-249", "150-199", "100-149", "50-99", "0-49"),
                        values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))
      ggsave("telo_length_bar_plot_by_strand.pdf", plot = telo_bar_hist_by_strand, device="pdf", width=6, height=10)

      telo_bar_hist_by_cluster <- ggplot(data=data_table) +
      geom_bar(mapping=aes(x=Cluster, fill=bin_telo_length), position="fill") +
      theme_minimal() +
      theme(axis.text = element_text(size=15),
            axis.title=element_text(size=20),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20)) +
      ylab("Proportion of Telomere Reads") + xlab("Strand") +
      guides(fill=guide_legend(title="Telomere Length (BP)")) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(breaks=c(525, 475, 425, 375, 325, 275, 225, 175, 125, 75, 25),
                        labels=c("500+", "450-499","400-449", "350-399", "300-349", "250-299", "200-249", "150-199", "100-149", "50-99", "0-49"),
                        values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))
      telo_bar_hist_by_cluster
      ggsave("telo_length_bar_plot_by_cluster.pdf", plot = telo_bar_hist_by_cluster, device="pdf", width=10, height=10)


      #D81B60
      #1E88E5
      #FFC107
      #004D40
} else if (args[2] == "comparison"){

  data_table <- read.table("combined.stats.txt", header=TRUE, sep="\t")
  
  final_telo_counts <- data_table %>% group_by(Sample, strand) %>%
        summarize(n=n())
  final_count_comparison <- ggplot(data=final_telo_counts) +
    geom_bar(mapping=aes(x=Sample, y=n, fill=strand), stat="identity") +
    theme_minimal() +
    xlab("Sample") + ylab("Telomere Read Count") +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_text(size=20),
          legend.text=element_text(size=15),
          axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    scale_fill_manual(breaks=c("C", "G"), values=c("#1E88E5", "#FFC107"), name="Strand")
  ggsave("telomere_count_by_sample.pdf", device="pdf", width=12, height=10, plot=final_count_comparison)
  
  telo_length_boxplot_comparison <- ggplot(data=data_table) +
    geom_boxplot(mapping=aes(x=Sample, y=telo_length, fill=Sample)) +
    theme_minimal() +
    xlab("Sample") + ylab("Telomere Lenght (bp)") +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_text(size=20),
          legend.text=element_text(size=15),
          axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  ggsave("telomere_length_boxplot_by_sample.pdf", device="pdf", width=12, height=10, plot=telo_length_boxplot_comparison)
  
  telo_lengths_for_binning <- c(525, 475, 425, 375, 325, 275, 225, 175, 125, 75, 25)
  
  # 0-50, 50-100, 150-200, 200-250, 250-300, 30-350, 400-450, 450-500, 500-550, 550-600, 600+
  data_table$bin_telo_length <- unlist(lapply(data_table$telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
  data_table$bin_telo_length <- factor(data_table$bin_telo_length, levels=c(525, 475, 425, 375, 325, 275, 225, 200, 175, 125, 75, 25))
  
  telo_bar_hist <- ggplot(data=data_table) +
    geom_bar(mapping=aes(x=Sample, fill=bin_telo_length), position="fill") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text(size=15),
          axis.title=element_text(size=20),
          legend.text=element_text(size=15),
          legend.title=element_text(size=20)) +
    ylab("Proportion of Telomere Reads") +
    guides(fill=guide_legend(title="Telomere Length (BP)")) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(breaks=c(525, 475, 425, 375, 325, 275, 225, 175, 125, 75, 25),
                      labels=c("500+", "450-499","400-449", "350-399", "300-349", "250-299", "200-249", "150-199", "100-149", "50-99", "0-49"),
                      values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))

  ggsave("telo_length_bar_plot_by_sample.pdf", plot = telo_bar_hist, device="pdf", width=6, height=10)
  
}