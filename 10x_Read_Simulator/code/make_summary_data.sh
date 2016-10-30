#!/bin/bash

input_file="$1" # depth coverage file
chrom_list="$2" # list of chrom's to use

cat $chrom_list | while read i; do
chrom="${i}"
echo $chrom
done