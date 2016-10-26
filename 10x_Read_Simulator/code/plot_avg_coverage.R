#!/usr/bin/env Rscript
# R 3.3

# USAGE: plot_avg_coverage.R data/depth.averages.txt /path/to/outdir

# DESCRIPTION:
# This script will read in a text file generated with the 'average_coverages.py' script
# and output a plot showing the average coverage per chromosome per genome (column)

library("reshape2")
library("ggplot2")


# get commands passed to the script
args <- commandArgs(TRUE)

cov_file <- args[1]
outdir <- args[2]
# save.image(file=file.path(outdir, "plot_avg_args.Rdata"),compress = TRUE)
# load("../test_output/plot_avg_args.Rdata")

# read in the file
coverage_df <- read.table(cov_file)

# fix colnames
colnames(coverage_df)[1] <- "chrom"
colnames(coverage_df)[-1] <- paste("genome_", seq_along(colnames(coverage_df)[-1]), sep="")

# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = "chrom", value.name = "coverage", variable.name = "sample")


# peel off the coverage stats column and turn into a grouping factor
stat_strings <- strsplit(as.character(coverage_df$coverage), ',')

stats_df <- data.frame(matrix(as.numeric(unlist(stat_strings)), nrow=length(stat_strings), byrow=T))
colnames(stats_df) <- c("average", "std_dev")

coverage_df <- cbind(coverage_df[c("chrom", "sample")], stats_df)
# colnames(stats_df) <- c("chrom", "sample", "average", "std_dev")

# melt again
coverage_df <- reshape2::melt(coverage_df, id.vars = c("chrom","sample"), value.name = "coverage", variable.name = "statistic")

# fix chrom order for plot
coverage_df <- coverage_df[with(coverage_df, order(chrom)), ]

# coverage_df
# make horizontal stacked grouped barplot
# plot by chrom
# pdf(file = file.path(outdir, "avg_cov_byChrom.pdf"), height = 8, width = 8)
# ggplot(coverage_df, aes(x = chrom, y = coverage, fill = factor(sample))) +
#   geom_bar(stat="identity", position="dodge") + # remove 'position' for stacked plot
#     coord_flip() + 
#     labs(title="Average Coverage Per Chromosome\nPer Samples", x="Chromosome", y = "Average Coverage")
# dev.off()

# plot by genome
coverage_df_avg <- subset(coverage_df, statistic == "average")
pdf(file = file.path(outdir, "avg_cov_byGenome.pdf"), height = 8, width = 8)
ggplot(coverage_df_avg, aes(x = sample, y = coverage, fill = factor(chrom))) +
  geom_bar(stat="identity", position="dodge") + 
    coord_flip() + 
    # scale_y_continuous(breaks = waiver()) + 
    # coord_cartesian(ylim=c(0,max(numeric(coverage_df_avg[["coverage"]])))) +
    labs(title="Average Coverage Per Chromosome\nPer Samples", x="Sample", y = "Average Coverage", fill="Chromosome")
    # need to fix axis scale markers !
dev.off()
save.image(file=file.path(outdir, "plot_avg.Rdata"),compress = TRUE)
# load("plot_avg.Rdata")
