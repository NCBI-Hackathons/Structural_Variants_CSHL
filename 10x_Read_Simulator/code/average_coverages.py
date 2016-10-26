#!/usr/bin/env python
# python 2.7

# USAGE: average_coverages.py data/depth_file.txt > output_file

# DESCRIPTION:
# This script will read in a coverage file and output the average coverage
# for chromosome per genome (column) in the file

# INPUT FILE FORMAT:
# chr1	10003	6656	7998	8002	7996	7997	7998	7421	8000	7999	5620	7999	6887	7998
# chr1	10004	7997	7999	8000	7997	7996	7997	7996	8000	8000	7786	8000	7998	7998
# chr1	10005	8000	8000	8005	7998	7999	8002	7999	8002	8001	8000	8001	8001	8000
# chr1	10006	8001	8004	8005	8003	8003	8005	8002	8003	8002	8000	8003	8002	8003
# chr1	10007	8001	8004	8006	8003	8000	8005	8002	8001	7994	7998	8000	7992	8000
# chr1	10008	7996	8000	8002	7994	7990	7999	7992	7985	7991	7995	7997	7991	7989
# chr1	10009	8005	8007	8009	8006	8006	8007	8005	8006	8004	8003	8006	8005	8005

# EXAMPLE OUTPUT:
# chr4	120	110	213	213	165	275	157	155	147	159	281	252	152
# chr3	688	616	908	916	729	972	642	743	661	754	979	888	742
# chr2	36	40	83	105	46	117	47	48	40	49	103	71	47
# chr1	599	587	680	683	619	689	597	615	588	631	678	701	634

from __future__ import division
import sys
import csv
import collections
import math

input_file = sys.argv[1]

def get_num_cols(infile):
    # get the number of columns in the file
    with open(infile) as file:
        reader = csv.reader(file, delimiter='\t')
        first_row = next(reader)
        num_cols = len(first_row)
    return num_cols

def genome_avg_coverages(infile):
    # get average coverage per chromosome per genome
    # EXAMPLE:
    # {genome1: {chr1: 10000, chr2:5000}, genome2:etc.}
    # ~~~~~ # 
    # ~~ SETUP ~~~ # 
    # get the number of genomes in the file
    num_genomes = get_num_cols(infile) - 2
    # dict to hold total number entries for each genome
    # n, s0
    genome_counts = collections.defaultdict(dict)
    # dict to hold total coverage for each genome
    # s1
    genome_coverages = collections.defaultdict(dict)
    # dict to hold the sum of squares of coverage for each genome
    # s2
    genome_SS_coverages = collections.defaultdict(dict)
    # dict to hold average coverage for each genome
    genome_average_coverages = collections.defaultdict(dict)
    # dict to hold the std dev of coverage for each genome
    genome_std_coverages = collections.defaultdict(dict)
    # initialize defaultdict for every genome in the file
    for i in range(1, num_genomes + 1):
        genome_coverages[i] = collections.defaultdict(int)
        genome_counts[i] = collections.defaultdict(int)
        genome_average_coverages[i] = collections.defaultdict(int)
        genome_SS_coverages[i] = collections.defaultdict(int)
    # ~~~ READ FILE ~~ # 
    # calculate the total & SS coverages 
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            chrom = line.pop(0)
            position = line.pop(0)
            for i in range(1, len(line) + 1):
                genome_counts[i][chrom] += 1
                genome_coverages[i][chrom] += int(line[i - 1])
                genome_SS_coverages[i][chrom] += int(line[i - 1]) * int(line[i - 1])
    # ~~~ STATS ~~ # 
    # calculate the averages & std dev
    for genome in genome_counts.iterkeys():
        for chrom in genome_counts[genome].iterkeys():
            s0 = genome_counts[genome][chrom]
            s1 = genome_coverages[genome][chrom]
            s2 = genome_SS_coverages[genome][chrom]
            # average
            genome_average_coverages[genome][chrom] = s1 / s0
            # std dev; watch for count of 1 = divide by zero
            if (s0 <= 1):
                genome_std_coverages[genome][chrom] = 0
            else :
                # sample std
                # genome_std_coverages[genome][chrom] = math.sqrt((s0 * s2 - s1 * s1)/(s0 * (s0 - 1)))
                # population std
                genome_std_coverages[genome][chrom] = math.sqrt((s0 * s2 - s1 * s1)/(s0 * (s0)))
    # ~~~ REFORMAT ~~ # 
    # format for printing to stdout; need to keep columns & entries in order !!
    chrom_list = [genome_average_coverages, genome_std_coverages]
    chrom_output = collections.defaultdict(list)
    # for k in d1.iterkeys():
    # d[k] = tuple(d[k] for d in ds)
    chrom_output_avg = collections.defaultdict(list)
    chrom_output_sd = collections.defaultdict(list)
    for genome in sorted(genome_average_coverages.keys()):
        for chrom in sorted(genome_average_coverages[genome].keys()):
            stats_tup = (genome_average_coverages[genome][chrom], genome_std_coverages[genome][chrom])
            chrom_output[chrom].append(stats_tup)
            # for av_value in chrom:
            #     value = value
            # chrom_output_avg[chrom].append(str(genome_average_coverages[genome][chrom]))
            #  + ','.join(str(genome_std_coverages[genome][chrom])
    for genome in sorted(genome_std_coverages.keys()):
        for chrom in sorted(genome_std_coverages[genome].keys()):
            chrom_output_sd[chrom].append(genome_std_coverages[genome][chrom])
    return chrom_output
    # return chrom_output_avg, chrom_output_sd

if __name__ == '__main__':
    chrom_output_avg = genome_avg_coverages(input_file)
    # chrom_output_avg, chrom_output_sd = genome_avg_coverages(input_file)
    # print chrom_output_avg.items()
    # print '' 
    # print chrom_output_sd.items()
    for chrom in sorted(chrom_output_avg.keys()):
        chrom_stats = chrom_output_avg[chrom]
        # print chrom + '\t' + '\t'.join(map(str,chrom_stats))
        # ["%s=%s" % (k, v) for k, v in params.items()]
        print chrom + '\t' + '\t'.join(map(str,["%s,%s" % (av, sd) for av, sd in chrom_stats]))

