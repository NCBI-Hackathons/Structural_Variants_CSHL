#!/usr/bin/env Rscript
# R 3.3

# USAGE: plot_bin_coverages.R data/depth.bin_coverages.txt /path/to/outdir

# DESCRIPTION:
# This script will plot the binned coverage data

# ! ! update with ggplot2 ribbon plot http://docs.ggplot2.org/current/geom_ribbon.html
# plot avg +/- sd

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

coverage_df['region'] <- paste(coverage_df[[2]], coverage_df[[3]], sep = '-')

coverage_df <- coverage_df[-(2:3)]


# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = c("chrom", "region"), value.name = "total_coverage", variable.name = "sample")

# fix chrom order for plot
coverage_df <- coverage_df[with(coverage_df, order(chrom, region, sample)), ]

# head(coverage_df)

pdf(file = file.path(outdir, "total_cov_byRegion.pdf"), height = 8, width = 8, onefile = TRUE)
for(i in seq_along(levels(coverage_df[['chrom']]))){
    ichrom <- levels(coverage_df[['chrom']])[i]
    cov_subset <- subset(coverage_df, chrom == ichrom)
    # print(head(cov_subset))
    myplot <- ggplot(cov_subset, aes(x = total_coverage))
    myplot <- myplot + geom_density(alpha=.5, fill="#FF6666") # geom_histogram(alpha=.5, position="identity") + 
    myplot <- myplot + labs(title=paste0("Chromosome: ", ichrom, "\nTotal Coverage"), x="Total Coverage", y = "Density")
    print(myplot)

}
dev.off()

pdf(file = file.path(outdir, "total_cov_byRegion_hist.pdf"), height = 8, width = 8, onefile = TRUE)
for(i in seq_along(levels(coverage_df[['chrom']]))){
    ichrom <- levels(coverage_df[['chrom']])[i]
    cov_subset <- subset(coverage_df, chrom == ichrom)
    # print(head(cov_subset))
    myplot <- ggplot(cov_subset, aes(x = total_coverage, fill = factor(sample)))
    myplot <- myplot + geom_histogram(alpha=.5, position="identity")
    myplot <- myplot + labs(title=paste0("Chromosome: ", ichrom, "\nTotal Coverage"), x="Total Coverage", y = "Count")
    print(myplot)

}
dev.off()

# pdf(file = file.path(outdir, "total_cov_byRegion.pdf"), height = 8, width = 8)
# dev.off()