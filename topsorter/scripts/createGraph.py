#!/usr/bin/env python

## ----------------------
# A python module for creating a graph from parsed VCF files
## ----------------------
# Han Fang (hanfang.cshl@gmail.com)
# Cold Spring Harbor Laboratory
## ----------------------

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import open
import os
import sys
import argparse
import vcf
import collections
import networkx as nx
import uuid


# the class for the data structure
class TreeNode(object):
    def __init__(self, chr):
        self.chr = chr
        self.start = None
        self.end = None
        self.next = []

# the class for filtering scalpel vcf file
class parseVcf:
    'filter based on parental coverage, chi2 score, and alternative allele coverage'
    def __init__(self, in_vcf):# , mc, ac, zyg, cr, chi2, fp, out_vcf):
        self.in_vcf  = in_vcf
        # self.sub_vcf = []
        self.large_variants = []
        self.chr_sizes = {}

    # a read function for reading and filtering
    def readVcf(self):
        # read a vcf file
        self.vcf_reader = vcf.Reader(open(self.in_vcf, 'r'))

        # [temp] filter out variants without END
        for variant in self.vcf_reader:
            if 'END' in variant.INFO:
                if abs(variant.INFO['END'] - variant.POS) > 10000:
                    self.large_variants.append(variant)

        # create a dict of chromosome: size
        for k,v in self.vcf_reader.contigs.items():
            self.chr_sizes[k] = v[1]

    def exportVcfBed(self):
        # export the subset of large SVs to a vcf
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.in_vcf + ".large_sv.vcf", 'w'), self.vcf_reader)
        for record in self.large_variants:
            vcf_writer.write_record(record)

        # export the flanking regions as a bed file
        bed_str = ""
        for var in self.large_variants:
            chr = var.CHROM
            before = chr + "\t" + str(int(var.POS)-1001) + "\t" + str(int(var.POS)-1) + "\t" + var.ID + ",before\n"
            left = chr + "\t" + str(int(var.POS)-1) + "\t" + str(int(var.POS)+999) + "\t" + var.ID + ",left\n"
            right = chr + "\t" + str(int(var.INFO['END'])-1001) + "\t" + str(int(var.INFO['END'])-1) + "\t" + var.ID + ",right\n"
            after = chr + "\t" + str(int(var.INFO['END'])-1) + "\t" + str(int(var.INFO['END'])+999) + "\t" + var.ID + ",after\n"
            bed_str += before + left + right + after

        # write to local
        bed_out = open(self.in_vcf + ".bed", "w")
        bed_out.write(bed_str)

    def createNodes(self):
        # sort the variants by start coordinates
        self.large_variants = sorted(self.large_variants, key=lambda x: x.POS)

        # scan the chromosomes and split by variants
        self.graph = collections.defaultdict(list)
        for chr in self.chr_sizes.keys():
            i = 0
            for variant in self.large_variants:
                if variant.CHROM == chr:
                    # nodes for left of SV
                    prev_val = (variant.CHROM, i, variant.POS-2, "prev")
                    prev_node = TreeNode(prev_val)
                    prev_node.start, prev_node.end = i, variant.POS-2
                    if i != 0:
                        curr_node.next.append(prev_node)
                    # current SV node
                    curr_val = (variant.CHROM, variant.POS-1, variant.INFO['END']-1, "curr")
                    curr_node = TreeNode(curr_val)
                    curr_node.start, curr_node.end = variant.POS-1, variant.INFO['END']-1
                    # link the prev to curr node
                    prev_node.next.append(curr_node)
                    self.graph[chr].extend([prev_node, curr_node])
                    i = variant.INFO['END']
            if self.graph[chr] != []:
                self.graph[chr][-1].next = [None]

        #for chr in graph.keys():
        #    for i in graph[chr]:
        #        if i.next is not None:
        #            print (i.chr, i.start, i.end, i.next[0])

    def weightedDag(self):
        DG=nx.DiGraph()
        example = self.graph["chr20"]
        for i in range(len(example)-1):
            DG.add_edge(example[i], example[i+1])
        # DG.add_weighted_edges_from([("A","B", 10000), ("B","C",0.75)])
        # topological sorting
        print ("[execute]\tPerforming topological sorting", flush=True)
        order = nx.topological_sort(DG)
        print ("[status]\tOrder of the nodes after topological sorting: ", flush=True)
        print (order)

        # build a tree
        start = order[0]
        nodes = [order[0]] # start with first node in topological order
        labels = {}
        print ("[status]\tEdges:")
        tree = nx.Graph()
        while nodes:
            source = nodes.pop()
            labels[source] = source
            for target in DG.neighbors(source):
                if target in tree:
                    t = uuid.uuid1() # new unique id
                else:
                    t = target
                labels[t] = target
                tree.add_edge(source,t)
                print (source, target, source, t)
                nodes.append(target)
        plt.figure()
        nx.draw(tree)# ,labels=labels)
        plt.gcf()
        plt.savefig("test.pdf")

    # write function for exporting vcf files
    def parse(self):
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.out_vcf, 'w'), self.vcf_reader)
        elif int(sys.version_info.major) == 2 :
            vcf_writer = vcf.Writer(open(self.out_vcf, 'wb'), self.vcf_reader)

        for record in self.sub_vcf:
            vcf_writer.write_record(record)

# the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='createGraph.py - a script for constructing a graph from a vcf file')
    parser.add_argument("-i", help="input vcf file [REQUIRED]",required=True)

    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        print ("[help] please use [-h] for details")
        sys.exit(1)
    else:
        args = parser.parse_args()

    # process the file if the input files exist
    if (args.i!=None):
        worker = parseVcf(args.i)
        print("[status]\tReading the vcf file: " + str(args.i), flush=True)
        worker.readVcf()
        print("[execute]\tExporting the large SVs to vcf and flanking regions to a bed file", flush=True)
        worker.exportVcfBed()
        print("[execute]\tCreating the nodes", flush=True)
        worker.createNodes()
        print("[execute]\tConstructing the graph", flush=True)
        worker.weightedDag()

    # print usage message if any argument is missing
    else:
        print ("[error]\tmissing argument")
        parser.print_usage()