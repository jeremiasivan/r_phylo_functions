## Tree Visualization
---
The following R script can be used to read and plot multitrees in Newick format.

<br>

### Usage (based on *example*)
```
$ cd example/
$ Rscript ../visualize_trees.R -f input.treefile -o output
```

For further details, you can run:
```
$ Rscript visualize_trees.R --help
```
<br>

### Notes
- If the treefile has no edge length, the labels will show number of children nodes by default
- The code is developed and tested in macOS Monterey v12.6 and R v4.1.3 (installed in Conda environment)

<br>

### Dataset
<p>Edelman, N. B., Frandsen, P. B., Miyagi, M., Clavijo, B., Davey, J., Dikow, R. B., . . . Mallet, J. (2019). Genomic architecture and introgression shape a butterfly radiation. <i>Science, 366</i>(6465), 594-599. doi:<a href="https://doi.org/10.1126/science.aaw2090">10.1126/science.aaw2090</a></p>