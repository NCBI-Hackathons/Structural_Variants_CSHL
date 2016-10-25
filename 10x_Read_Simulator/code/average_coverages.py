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
    # get the number of genomes in the file
    num_genomes = get_num_cols(infile) - 2
    # dict to hold total coverage for each genome
    genome_coverages = collections.defaultdict(dict)
    # dict to hold total number entries for each genome
    genome_counts = collections.defaultdict(dict)
    # dict to hold average coverage for each genome
    genome_average_coverages = collections.defaultdict(dict)
    for i in range(1, num_genomes + 1):
        genome_coverages[i] = collections.defaultdict(int)
        genome_counts[i] = collections.defaultdict(int)
        genome_average_coverages[i] = collections.defaultdict(int)
    # calculate the total coverages
    with open(infile) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            chrom = line.pop(0)
            position = line.pop(0)
            for i in range(1, len(line) + 1):
                genome_counts[i][chrom] += 1
                genome_coverages[i][chrom] += int(line[i - 1])
    # calculate the averages
    for genome in genome_counts.iterkeys():
        for chrom in genome_counts[genome].iterkeys():
            genome_average_coverages[genome][chrom] = genome_coverages[genome][chrom] / genome_counts[genome][chrom]
    # format for printing to stdout; need to keep columns & entries in order !!
    chrom_output = collections.defaultdict(list)
    for genome in sorted(genome_average_coverages.keys()):
        for chrom in sorted(genome_average_coverages[genome].keys()):
            chrom_output[chrom].append(genome_average_coverages[genome][chrom])
    for chrom in sorted(chrom_output.keys()):
        print chrom + '\t' + '\t'.join(map(str,chrom_output[chrom]))
genome_avg_coverages(input_file)
