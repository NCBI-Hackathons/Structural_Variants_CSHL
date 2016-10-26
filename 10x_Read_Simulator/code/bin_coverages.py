#!/usr/bin/env python
# python 2.7

# USAGE: bin_coverages.py data/depth_file.txt

# DESCRIPTION:
# This script will parse a depth of coverage file into binned regions
# and calculate the average coverage per region

import sys
import csv
import collections
import average_coverages as av


def test_bins(value, start_bin, end_bin):
    if (int(value) < int(start_bin)) == True:
        print 'ERROR: value {} smaller than Start bin {}'.format(value, start_bin)
        sys.exit()
    if (int(value) > int(end_bin)) == True:
        print 'ERROR: value {} greater than end bin {}'.format(value, end_bin)
        sys.exit()



input_file = sys.argv[1]
chrom_sizes_file = "code/hg38_chrom_sizes.txt"
bin_size = 1000
chrom_coverage_counts = collections.defaultdict(int)
# chrom_sizes = collections.defaultdict(int)
chrom_bins = collections.defaultdict(lambda : collections.defaultdict(list))
# keep track of bins while parsing file
chrom_bin_counts = collections.defaultdict(lambda : collections.defaultdict(int))


if __name__ == '__main__':
    # # create bins from chom sizes # don't use this
    # with open(chrom_sizes_file) as chrom_file:
    #     chrom_file = csv.reader(chrom_file, delimiter='\t')
    #     for line in chrom_file:
    #         print line
    #         chrom = line.pop(0)
    #         size = line.pop(0)
    #         chrom_bins[chrom] = range(1, int(size))
    # print chrom_bins
    with open(input_file) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for line in tsvin:
            chrom = line.pop(0)
            position = line.pop(0)
            # convert remaining values to int
            line = map(int, line)
            # set position bins
            position_start_bin = int(position) - (int(position) % bin_size)
            position_end_bin = position_start_bin + bin_size
            position_bin = str(position_start_bin) + '-' + str(position_end_bin)
            # sanity check - make sure we binned correctly
            test_bins(position, position_start_bin, position_end_bin)
            # sanity check - make sure we don't add too many more entries to the bin
            chrom_bin_counts[chrom][position_bin] += 1 
            if chrom_bin_counts[chrom][position_bin] <= bin_size:
                if not chrom_bins[chrom][position_bin]:
                    chrom_bins[chrom][position_bin] = line
                else :
                    chrom_bins[chrom][position_bin] = [x + y for x, y in zip(chrom_bins[chrom][position_bin], line)]
            else :
                print 'ERROR: bin exceeded!'
                print chrom, position, position_bin, chrom_bin_counts[chrom][position_bin]
                sys.exit()
                # chrom_bins[chrom][position_bin] = [sum(x) for x in zip(*chrom_bins[chrom][position_bin])]
            # ~~~~~~ # # ~~~~~~ # # ~~~~~~ #
            # include shifted bin window
            # position_starthalf_bin = int(position) - (int(position) % bin_size) - (bin_size / 2)
            # position_endhalf_bin = position_starthalf_bin + bin_size
            # position_halfbin = str(position_starthalf_bin) + '-' + str(position_endhalf_bin)
            # test_bins(position, position_starthalf_bin, position_endhalf_bin)
            # ~~~~~~ # # ~~~~~~ # # ~~~~~~ # 
            print chrom, position, position_bin #, position_halfbin
    #
    for chrom, statsdict in chrom_bins.iteritems():
        for statstat in statsdict.iteritems():
            print chrom, statstat
    #
    # chrom_average_coverages = av.genome_avg_coverages(input_file)
    # print, for stdout piping
    # for chrom in sorted(chrom_average_coverages.keys()):
    #     print chrom + '\t' + '\t'.join(map(str,chrom_average_coverages[chrom]))
