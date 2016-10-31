#!/usr/bin/env Rscript
# R 3.3

# USAGE: plot_avg_coverage.R data/depth.averages.txt /path/to/outdir

# DESCRIPTION:
# This script will read in a text file generated with the 'average_coverages.py' script
# and output a plot showing the average coverage per chromosome per genome (column)

library("reshape2")
library("ggplot2")
library("plotly")

# get commands passed to the script
args <- commandArgs(TRUE)

cov_file <- args[1]
outdir <- args[2]
# save.image(file=file.path(outdir, "plot_avg_args.Rdata"),compress = TRUE)
# load("../test_output/plot_avg_args.Rdata")

# read in the file
coverage_df <- read.table(cov_file)
save.image(file=file.path(outdir, "plot_avg_start.Rdata"),compress = TRUE)

# fix colnames
colnames(coverage_df)[1] <- "chrom"
colnames(coverage_df)[-1] <- paste("genome_", seq_along(colnames(coverage_df)[-1]), sep="")

# melt into long format
coverage_df <- reshape2::melt(coverage_df, id.vars = "chrom", value.name = "coverage", variable.name = "sample")


# peel off the coverage stats column and turn into a grouping factor
stat_strings <- strsplit(as.character(coverage_df$coverage), ',')

stats_df <- data.frame(matrix(as.numeric(unlist(stat_strings)), nrow=length(stat_strings), byrow=T))
colnames(stats_df) <- c("average", "std_dev", "count")

coverage_df <- cbind(coverage_df[c("chrom", "sample")], stats_df)
# colnames(stats_df) <- c("chrom", "sample", "average", "std_dev")

# melt again
# coverage_df <- reshape2::melt(coverage_df, id.vars = c("chrom","sample"), value.name = "coverage", variable.name = "statistic")


# fix chrom order for plot
# coverage_df <- coverage_df[with(coverage_df, order(chrom)), ]


# plot by genome
# horizontal barplot per genome
# coverage_df_avg <- subset(coverage_df, statistic == "average")
# pdf(file = file.path(outdir, "avg_cov_byGenome.pdf"), height = 8, width = 8)
chrom_plot <- ggplot(coverage_df, aes(x = sample, y = average, fill = factor(chrom)))
chrom_plot <-chrom_plot + geom_bar(stat="identity", position="dodge")
chrom_plot <-chrom_plot + coord_flip()
chrom_plot <-chrom_plot + labs(title="Average Coverage Per Chromosome\nPer Samples", x="Sample", y = "Average Coverage", fill="Chromosome")
print(chrom_plot)
# dev.off()
quit()




# ribbon plot
# coverage_df_avg <- subset(coverage_df, statistic == "average" | statistic == "std_dev")
# coverage_sample_df <- subset(coverage_df, sample == unique(as.character(coverage_df[["sample"]]))[1] & statistic == "average")
chrom_ribbon <- ggplot(coverage_df, x=chrom, x=average, colour=sample)
chrom_ribbon <- chrom_ribbon + geom_line()
chrom_ribbon <- chrom_ribbon + geom_ribbon(aes(x=chrom, ymin=average - std_dev, ymax=average + std_dev), alpha = 0.3)

print(chrom_ribbon)


# save.image(file=file.path(outdir, "plot_avg.Rdata"),compress = TRUE)
# load("plot_avg.Rdata")
print(coverage_sample_df)
quit()






# test
# coverage_df_sb <- subset(coverage_df, chrom == "chr1" | chrom == "chr2" | chrom == "chr3" | chrom == "chr4") 
# coverage_df_sb <- subset(coverage_df_sb, sample == "genome_1" | sample == "genome_2" )#
# coverage_df_sb <- droplevels(coverage_df_sb)
# levels(coverage_df_sb$sample) <- c("sample_1", "sample_2")
# coverage_df_sb <- coverage_df_sb[c("chrom", "sample", "average", "std_dev")]
# colnames(coverage_df_sb)
# dput(coverage_df_sb)

sample_data <- structure(list(chrom = structure(c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L), 
                                                .Label = c("chr1", "chr2", "chr3", "chr4"), 
                                                class = "factor"), 
                              sample = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L), 
                                                 .Label = c("sample_1", "sample_2"), 
                                                 class = "factor"),
                              average = c(358.017, 34.452, 409.7959, 117.0805, 345.6717, 34.3544, 362.3519, 110.7264), 
                              std_dev = c(1484.33280699, 97.332895241, 1460.24099656, 
                                          519.299214731, 1439.86318396, 114.04659662, 1340.67100158, 
                                          499.901605662)), 
                         .Names = c("chrom", "sample", "average", "std_dev"), 
                         row.names = c(1L, 2L, 3L, 4L, 25L, 26L, 27L, 28L), class = "data.frame")

sample_ribbon <- ggplot(sample_data, x=chrom, y=average, colour=sample)
sample_ribbon <- sample_ribbon + geom_line()
sample_ribbon <- sample_ribbon + geom_ribbon(aes(x=chrom, ymin=average - std_dev, ymax=average + std_dev), alpha = 0.3)
# Error in order(data$PANEL, data$group, data$x) : 
# argument 3 is not a vector
print(sample_ribbon)









set.seed(1)
y <- sin(seq(1, 2*pi, length.out = 100))
x <- 1:100
plotdata <- data.frame(x=x, y=y, lower = (y+runif(100, -1, -0.5)), upper = (y+runif(100, 0.5, 1)))




# plotly
chrom_plotly <- ggplotly(chrom_plot)
htmlwidgets::saveWidget(as.widget(chrom_plotly), file.path(outdir, "avg_cov_byGenome.html"))


# coverage_df
# make horizontal stacked grouped barplot
# plot by chrom
# pdf(file = file.path(outdir, "avg_cov_byChrom.pdf"), height = 8, width = 8)
# ggplot(coverage_df, aes(x = chrom, y = coverage, fill = factor(sample))) +
#   geom_bar(stat="identity", position="dodge") + # remove 'position' for stacked plot
#     coord_flip() + 
#     labs(title="Average Coverage Per Chromosome\nPer Samples", x="Chromosome", y = "Average Coverage")
# dev.off()
