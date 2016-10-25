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

# read in the file
coverage_df <- read.table(cov_file)

# fix colnames
colnames(coverage_df)[1] <- "chrom"
colnames(coverage_df)[-1] <- paste("genome_", seq_along(colnames(coverage_df)[-1]), sep="")

# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = "chrom", value.name = "avg_coverage", variable.name = "sample")

# fix chrom order for plot
coverage_df <- coverage_df[with(coverage_df, order(chrom)), ]


# make horizontal stacked grouped barplot
# plot by chrom
# pdf(file = file.path(outdir, "avg_cov_byChrom.pdf"), height = 8, width = 8)
# ggplot(coverage_df, aes(x = chrom, y = avg_coverage, fill = factor(sample))) +
#   geom_bar(stat="identity", position="dodge") + # remove 'position' for stacked plot
#     coord_flip() + 
#     labs(title="Average Coverage Per Chromosome\nPer Samples", x="Chromosome", y = "Average Coverage")
# dev.off()

# plot by genome
pdf(file = file.path(outdir, "avg_cov_byGenome.pdf"), height = 8, width = 8)
ggplot(coverage_df, aes(x = sample, y = avg_coverage, fill = factor(chrom))) +
  geom_bar(stat="identity", position="dodge") + 
    coord_flip() + 
    labs(title="Average Coverage Per Chromosome\nPer Samples", x="Sample", y = "Average Coverage")
dev.off()
