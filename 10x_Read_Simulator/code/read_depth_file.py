#!/usr/bin/env python
# python 2.7

import sys
import csv
import collections

input_file = sys.argv[1]
print input_file
# sys.exit()


def get_num_cols(infile):
    # get the number of columns in the file
    with open(infile) as file:
        reader = csv.reader(file, delimiter='\t')
        first_row = next(reader)
        num_cols = len(first_row)
    return num_cols

# number of genomes in depth file 
num_genomes = get_num_cols(input_file) - 2
print "Number of genomes", num_genomes


def get_avg_coverage(infile):
    # get the average coverage per genome in the file
    # num_genomes = get_num_cols(infile) - 2
    total_coverage = 0 
    line_count = 0
    # calculate the average coverage for a column in the file
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            line_count += 1
            coverage = line[2]
            total_coverage = total_coverage + int(coverage)
    avg_coverage = total_coverage / line_count
    return avg_coverage

avg_coverage = get_avg_coverage(input_file)
print 'Average Coverage', avg_coverage


def get_binned_stats(infile, avg_coverage):
    # separate the input into binned regions
    bin_size = 100 # size of the bins
    bin_count = 0 # bin iterator
    bin_coverage = 0 # starting coverage vale
    postion_stats = [] # placeholder for bin stats
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            if (bin_count <= bin_size):
                # get values from the file line
                chrom = line[0]
                position = line[1]
                bin_position = int(position) / int(bin_size)
                coverage = line[2]
                bin_coverage = int(coverage) + bin_coverage
                if (bin_count == 0):
                    # start of the position
                    position_values = []
                    position_values.append(chrom)
                    # position_values.append(position)
                    position_values.append(bin_position)
                if (bin_count == bin_size):
                    # end of the position
                    # bin_avg_coverage = bin_coverage / bin_count
                    # position_values.append(position)
                    position_values.append(bin_position)
                    # coverage for the binned region / avg cov whole file
                    position_values.append(bin_coverage / avg_coverage) 
                    postion_stats.append(position_values)
                    # RESET
                    bin_count = 0
                    bin_coverage = 0
                    continue
                bin_count += 1
    return postion_stats

postion_stats = get_binned_stats(infile = input_file, avg_coverage = avg_coverage)

# for line in postion_stats:
#     print('\t'.join(map(str,line)))




def genome_avg_coverages(infile):
    # get average coverage per chromosome per genome
    # EXAMPLE:
    # {genome1: {chr1: 10000, chr2:5000}, genome2:etc.}
    # ~~~~~ # 
    # get the number of genomes in the file
    num_genomes = get_num_cols(infile) - 2
    # dict to hold total coverage for each genome
    genome_coverages = collections.defaultdict(dict)
    # dict to hold total number entries for each genome
    genome_counts = collections.defaultdict(dict)
    # dict to hold average coverage for each genome
    genome_average_coverages = collections.defaultdict(dict)
    for i in range(1, num_genomes + 1):
        genome_coverages[str(i)] = collections.defaultdict(int)
        genome_counts [str(i)] = collections.defaultdict(int)
        genome_average_coverages [str(i)] = collections.defaultdict(int)
    # calculate the total coverages
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            print line
            chrom = line.pop(0)
            position = line.pop(0)
            for i in range(1, len(line) + 1):
                print i - 1 
                print line[i - 1]
                genome_counts[str(i)][chrom] += 1
                genome_coverages[str(i)][chrom] += int(line[i - 1])
    # print genome_counts.keys()
    # print genome_coverages.keys()
    # calculate the average coverages
    for genome in genome_counts.iterkeys():
        for chrom in genome_counts[genome].iterkeys():
            genome_average_coverages[genome][chrom] = genome_coverages[genome][chrom] / genome_counts[genome][chrom]
    # print the results
    for genome in genome_average_coverages.iterkeys():
        for chrom in genome_average_coverages[genome].iterkeys():
            print chrom, '\t', genome_average_coverages[genome][chrom]
genome_avg_coverages(input_file)
