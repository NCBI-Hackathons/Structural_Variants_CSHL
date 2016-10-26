## README
####  

1) topsorter.py

A script that creates a weighted directed acyclic graph (dag) for each chromosome in python using networkx. 
  1) It represents a consecutive region on the reference genome as a node 
  2) After the initial construction of the graph, the script updates the weights of the edges according to obtained barcodes from the script barcode_profiles.py and quality values.
  3) Finally the script sorts the graph to obtain the longest path in the weighted graph using a topological sorting algorithm.
  
2) barcode_profiling.py

 Given the list of all barcodes present in a regions left and right of each breakpoints for each SVs that are labeled as before, left, right and after , we generated a bed file that including all the 4 regions per SVs each of which 10kb long. 
 1) Next we used samtools to extract the reads overlapping these regions. 
 2) For each of these regions we collect the information of which barcodes occur and how many reads are supporting that barcode. 
 3) In the end a new bed file is generated holding the so obtained informations.

   


