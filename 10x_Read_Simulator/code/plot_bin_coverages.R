#!/usr/bin/env Rscript
# R 3.3

# USAGE: plot_bin_coverages.R data/depth.bin_coverages.txt /path/to/outdir

# DESCRIPTION:
# This script will plot the 

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
colnames(coverage_df)[2] <- "start"
colnames(coverage_df)[3] <- "stop"
colnames(coverage_df)[-(1:3)] <- paste("genome_", seq_along(colnames(coverage_df)[-(1:3)]), sep="")

head(coverage_df)

# melt into long format
# coverage_df <- reshape2::melt(coverage_df, id.vars = "chrom", value.name = "avg_coverage", variable.name = "sample")

# fix chrom order for plot
# coverage_df <- coverage_df[with(coverage_df, order(chrom)), ]



# # # # # # 
# library("reshape2")
# interval_summ_df <- dcast(interval_summ_df_sub, Target + sample ~ method, value.var = "average_coverage")
# interval_summ_df[["Difference"]] <- interval_summ_df[["KP"]] - interval_summ_df[["HP"]]

# ggplot(interval_summ_df, aes(x = Difference)) + geom_density(alpha=.5, fill="#FF6666") +
#     # geom_histogram(alpha=.5, position="identity") + 
#     labs(title="Difference in Avgerage Coverage KP vs. HP", x="Difference (KP - HP)", y = "Density")
