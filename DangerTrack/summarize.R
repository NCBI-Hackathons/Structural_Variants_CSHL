#!/usr/bin/env Rscript


# R commands for summarizing and plotting the data


####################


# load relevant libraries
library(GGally, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(reshape2, quietly = TRUE)


####################


# find all summary files
summary_files = list.files(path = ".", pattern = "summary.*.5000.txt")
summary_files

# import and merge all summary files (619,150 rows for 5kb bins)
imported_files = lapply(summary_files, read.delim, check.names = FALSE, stringsAsFactors = FALSE)
combined_summary = Reduce(function(x, y) merge(x, y, all = FALSE, by = "#BIN"), imported_files, accumulate = FALSE)
head(combined_summary)
dim(combined_summary)


####################


# plot distributions of values

png("dist.Mapability50mer.png", res = 200, width = 8, height = 5, units = "in")
ggplot(combined_summary, aes(Mapability50mer)) + geom_freqpoly(binwidth=0.01, size=2)
dev.off()

png("dist.gc.png", res = 200, width = 8, height = 5, units = "in")
ggplot(combined_summary, aes(GC_fraction)) + geom_freqpoly(binwidth=0.01, size=2)
dev.off()

png("dist.n.png", res = 200, width = 8, height = 5, units = "in")
ggplot(combined_summary, aes(N_fraction)) + geom_freqpoly(binwidth=0.01, size=2)
dev.off()

# boxplot of mapability versus number break points
png("map-breaks.box.png", res = 200, width = 8, height = 5, units = "in")
ggplot(subset(combined_summary, GIAB_H002_events < 21),
       aes(x = GIAB_H002_events, y = Mapability50mer, group = GIAB_H002_events)) +
geom_boxplot()
dev.off()

# violin plot of mapability versus number break points
png("map-breaks.violin.png", res = 200, width = 8, height = 5, units = "in")
ggplot(subset(combined_summary, GIAB_H002_events < 21),
       aes(x = GIAB_H002_events, y = Mapability50mer, group = GIAB_H002_events)) +
geom_violin(fill = "black")
dev.off()


####################


# a matrix of correlation plots with regression lines

# calculate correlation for a random subset of entries to save time (remove bin name column)
correlation_subset = combined_summary[sample(nrow(combined_summary), 10000), 2:ncol(combined_summary)]
correlation_subset = as.matrix(correlation_subset)
head(correlation_subset)
dim(correlation_subset)

# regression lines
smoothed_mean = function(data, mapping, ...) {
  p = ggplot(data = data, mapping = mapping) +
    geom_point() +
    geom_smooth(method="loess", fill="red", color="red", ...) +
    geom_smooth(method="gam", fill="blue", color="blue", ...)
  p
}

# plot correlations (full subset)
png("correlations.full.png", res = 100, width = 15, height = 10, units = "in")
ggpairs(correlation_subset, lower = list(continuous = smoothed_mean))
dev.off()

# filter for alignable bins to remove non-informative stats
correlation_subset_filtered = correlation_subset[correlation_subset[,"N_fraction"] < 0.5,]
head(correlation_subset_filtered)
dim(correlation_subset_filtered)

# plot correlations (filtered subset)
png("correlations.filtered.png", res = 100, width = 15, height = 10, units = "in")
ggpairs(correlation_subset_filtered, lower = list(continuous = smoothed_mean))
dev.off()


####################


# calculate score

scores = combined_summary
scores[,"chr"] = sub("(.*):.*-.*", "\\1",  scores[,"#BIN"])
scores[,"start_pos"] = as.integer(sub(".*:(.*)-.*", "\\1",  scores[,"#BIN"]))
scores[,"end_pos"] = as.integer(sub(".*:.*-(.*)", "\\1",  scores[,"#BIN"]))

scores[,"map_score"] = 1 - scores[,"Mapability50mer"]
scores[,"breaks_score"] = scores[,"GIAB_H002_events"] / max(scores[,"GIAB_H002_events"])

scores[,"dangerscore"] = scores[,"map_score"] + scores[,"breaks_score"]
scores[,"dangerscore"] = scores[,"dangerscore"] / 2

# chr3 dangerscore plot
png("dangerscore.chr3.png", res = 100, width = 15, height = 5, units = "in")
ggplot(subset(scores, chr == "chr3"), aes(y = dangerscore, x = start_pos)) +
geom_point(shape = 20, size = 5, alpha = 0.1) +
scale_x_continuous(labels = scales::comma)
dev.off()

# chr5 dangerscore plot
png("dangerscore.chr5.png", res = 100, width = 15, height = 5, units = "in")
ggplot(subset(scores, chr == "chr5"), aes(y = dangerscore, x = start_pos)) +
geom_point(shape = 20, size = 5, alpha = 0.1) +
scale_x_continuous(labels = scales::comma)
dev.off()

# export BED file
scores_bed = scores[, c("chr", "start_pos", "end_pos", "#BIN", "dangerscore")]
scores_bed[,"dangerscore"] = round(scores_bed[,"dangerscore"] * 1000, digits = 0)
write.table(scores_bed, file = "scores.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# export bedGraph file
scores_bedgraph = scores_bed[, c("chr", "start_pos", "end_pos", "dangerscore")]
write.table(scores_bedgraph, file = "scores.bedgraph", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


####################



# end
