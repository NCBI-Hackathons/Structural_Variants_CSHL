# 10x genomics read simulator
This method simulates 10x genomics reads in the following steps:

1. It simulates haplotypes introducing SNPs and SVs
2. It generates the fragments from the haplotypes and assigns barcodes
3. It uses a modified version of dwgsim to simulate illumina reads per fragment
