#!/usr/bin/env python

## ----------------------
# Topsorter
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


class topSorter:
    'parse & filter vcf files, construct & traverse DAGs'
    def __init__(self, in_vcf):
        self.in_vcf  = in_vcf
        # self.sub_vcf = []
        self.large_variants = []
        self.chr_sizes = {}

    # function for reading and filtering vcf files
    def readVcf(self):
        # read a vcf file
        self.vcf_reader = vcf.Reader(open(self.in_vcf, 'r'))

        # [temp] filter out variants without END
        hashset = set()
        for variant in self.vcf_reader:
            if 'END' in variant.INFO and 'PAIR_COUNT' in variant.INFO \
                and variant.INFO['SVTYPE'] in ("DEL", "DUP", "INV"):
                if abs(variant.INFO['END'] - variant.POS) > 10000 and variant.INFO['PAIR_COUNT'] >= 10:
                    id = (variant.CHROM, variant.POS, variant.INFO['END'])
                    if id not in hashset:
                        self.large_variants.append(variant)
                        hashset.add(id)

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
            before = chr +"\t"+ str(var.POS-1001) +"\t"+ str(var.POS-1) +"\t"+ var.ID + ",before\n"
            left   = chr +"\t"+ str(var.POS-1) +"\t"+ str(var.POS+999) +"\t"+ var.ID + ",left\n"
            right = chr +"\t"+ str(var.INFO['END']-1001) +"\t"+ str(var.INFO['END']-1) +"\t"+ var.ID + ",right\n"
            after = chr +"\t"+ str(var.INFO['END']-1) +"\t"+ str(var.INFO['END']+999) +"\t"+ var.ID + ",after\n"
            bed_str += before + left + right + after

        # write bed file to local
        bed_out = open(self.in_vcf + ".bed", "w")
        bed_out.write(bed_str)

    def createNodes(self):
        # sort the variants by start coordinates
        self.large_variants = sorted(self.large_variants, key=lambda x: x.POS)

        # scan the chromosomes and split by variants
        self.graph = collections.defaultdict(list)
        last_pos = collections.defaultdict(int)
        for chr in self.chr_sizes.keys():
            i = 0
            for variant in self.large_variants:
                if variant.CHROM == chr:
                    # add ref node
                    ref_node = (variant.CHROM, i, variant.POS-1, "REF")
                    self.graph[chr].append(ref_node)
                    # add variant node
                    var_node = (variant.CHROM, variant.POS-1, variant.INFO['END']-1, variant.INFO['SVTYPE'])
                    self.graph[chr].append(var_node)
                    i = variant.INFO['END']-1
                    last_pos[chr] = i # keep track of the last pos
        # add the last node
        for chr in self.chr_sizes.keys():
            last_node = (chr, last_pos[chr], self.chr_sizes[chr]-1, "REF")
            self.graph[chr].append(last_node)

        '''
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
        '''

    def weightedDAG(self):
        DAG = nx.DiGraph()
        chroms = self.graph.keys()
        # add nodes
        n = len(self.graph["chr20"])
        # print (n)
        example = self.graph["chr20"]
        # print(example)
        #for i in example:
        #    DAG.add_node(i)
        #DAG.add_nodes_from(example)
        #print(DAG.nodes())
        # initialize the graph with equal weight
        for i in range(len(example)-1):
            DAG.add_node(example[i]) #
            DAG.add_node(example[i+1]) #
            DAG.add_edge(example[i], example[i+1], weight=1)
            if example[i][3] == "DEL":
                DAG.add_edge(example[i-1], example[i+1], weight=30)
            elif example[i][3] == "DUP":
                nodeCopy = (example[i][0], example[i][1], example[i][2], "DUP_COPY")
                DAG.add_node(nodeCopy)
                DAG.add_edge(example[i], nodeCopy, weight=1)
                DAG.add_edge(nodeCopy, example[i+1], weight=1)
            elif example[i][3] == "INV":
                nodeInv = (example[i][0], example[i][2], example[i][1], "INV_FLIP")
                DAG.add_edge(example[i-1], nodeInv, weight=1)
                DAG.add_edge(nodeInv, example[i+1], weight=1)
        # print(nx.find_cycle(DAG))# .edges())
        print(DAG.nodes())
        print(DAG.edges())

        # topological sorting
        print ("[execute]\tPerforming topological sorting", flush=True)
        order = nx.topological_sort(DAG)
        print ("[status]\tOrder of the nodes after topological sorting: ", flush=True)
        print ("[results]\t", order, flush=True)
        longest = nx.dag_longest_path(DAG)
        print ("[results]\t Longest\n", longest, flush=True)
        longestDis = nx.dag_longest_path_length(DAG)
        print ("[results]\t Longest length: ", longestDis, flush=True)

        '''
        # build a tree
        start = order[0]
        nodes = [order[0]] # start with first node in topological order
        labels = {}
        print ("[status]\tEdges:")
        tree = nx.Graph()
        while nodes:
            source = nodes.pop()
            labels[source] = source
            for target in DAG.neighbors(source):
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
        '''
        plt.figure()
        nx.draw(DAG)# ,labels=labels)
        plt.gcf()
        plt.savefig("test_new.pdf")

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
        worker = topSorter(args.i)
        print("[status]\tReading the vcf file: " + str(args.i), flush=True)
        worker.readVcf()
        print("[execute]\tExporting the large SVs to vcf and flanking regions to a bed file", flush=True)
        worker.exportVcfBed()
        print("[execute]\tCreating the nodes", flush=True)
        worker.createNodes()
        print("[execute]\tConstructing the graph", flush=True)
        worker.weightedDAG()

    # print usage message if any argument is missing
    else:
        print ("[error]\tmissing argument")
        parser.print_usage()