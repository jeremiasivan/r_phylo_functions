## Window Tree Annotation
<p>The following R script can be used to generate concatenated alignment from multiple fasta files, build windows trees, and annotate based on the tree topology and chromosomal position. In addition to R libraries, <a href="https://bioinf.shenwei.me/seqkit/">Seqkit</a> and <a href="http://www.iqtree.org">IQTree2</a> are also used for specific analyses, and users should provide the executable themselves.</p>

### Usage (based on *example*)
```
$ cd example/
$ Rscript ../annotate_windows.R -a 1 --chr chr21 --fasta fastas/ --topology input_topology.txt --seqkit seqkit.exe --iqtree2 iqtree2.exe -o output
```

For further details, including input files for each analysis, you can run:
```
$ Rscript annotate_windows.R --help
```

### Known Error(s)
- By default, chromosomal annotation uses `Set1` color pallete from `RColorBrewer`, which will return a warning if there are more than 8 tree topologies.
- Construction of window trees (i.e. `-a 3`) might return an error if there is a lot of window alignments. Please use *IQTree2* from command-line if that happens.

### Notes
- Annotation only works for a single chromosome
- The code is developed and tested in macOS Monterey v12.6 and R v4.1.3 (installed in Conda environment)

### Dataset
<p>Edelman, N. B., Frandsen, P. B., Miyagi, M., Clavijo, B., Davey, J., Dikow, R. B., . . . Mallet, J. (2019). Genomic architecture and introgression shape a butterfly radiation. <i>Science, 366</i>(6465), 594-599. doi:<a href="https://doi.org/10.1126/science.aaw2090">10.1126/science.aaw2090</a></p>