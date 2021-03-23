# Changes Completed
Snakefile changes from original (in Master) are now in minimal_snakefile.py, callable with the -s flag. Changes to the snakefile are:

1. Removed rule.NanoFilt. NanoFilt is hanging, and also isn't strictly necessary: between cutadapt and isonclust we can achieve the same filtering. Cutadapt filters by length. isONclust filters by quality.
2. Skipped purging small clusters. Couldn't get this work either with bash or python. It also isn't necessary - normal downstream processing of the OTU table will accomplish the same goal.
3. Removed minimap_from_filter. Without NanoFilt, this doesn't make sense and it also doesn't produce useful results.
4. Removed related count_reads and histograms related to minimap_from_filter and NanoFilt.
5. Remove the box-and-whisker plot, because it kept throwing an error within python. If fixable, this would be a helpful graph.
6. Explicitly parallelize cutadapt, isonclust, SPOA, minimap2. These programs all take a parameter specifying resource availability, with defaults. Snakemake now tells them exactly how many threads are available. Also added an argument increasing memory usage by minimap2.

7. Added a rule to create a minimap index from the reference database as a separate step. Best practice for large indices so can recover with less computation.


# Changes Necessary

1. Drop basecalling? I'm not sure how to get snakemake, singularity, and CUDA to communicate. When explicitly opening singularity in a bash, need to add a CUDA flag to make it available, as in:
```singularity --nv shell /path/to/container```
Then assume that basecalling and barcoding have already been done before executing the pipeline. 
2. *Alternatively* if there is a good way to get this communication, then great.

2. Include GCC as a SPOA dependency to the conda environment. On scinet, need to explicitly call `module load gcc`. Probably should have a SPOA container.

3. Better than conda is to build a full pipeline image using `snakemake --containerize > Dockerfile`?
