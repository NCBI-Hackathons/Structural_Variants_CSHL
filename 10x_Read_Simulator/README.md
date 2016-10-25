

```
# make the averages file
code/average_coverages.py data
/depth13Genome.depth.summary2.txt > sample_output/depth13Genome.depth.summary2.averages.txt

# make the diagnostic plot of the averages
code/plot_avg_coverage.R sample_output/depth13Genome.depth.summary2.averages.txt sample_output/
```

Example output

!()[sample_output/avg_cov_byGenome.png]

# Dependencies

- Python 2.7

- R 3.3

- - ggplot2

- - reshape2