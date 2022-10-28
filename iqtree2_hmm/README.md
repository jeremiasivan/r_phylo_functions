## Comparison of MAST+HMM and Sliding Windows
The following R script can be used to run HMM and compare the results of MAST and sliding-windows techniques

### HMM Installation via Conda
```
$ conda create hmm
$ conda activate hmm
$ conda install -c conda-forge r-essentials r-devtools r-testit r-kmer
```

Then, install the rest of the dependencies via R
```
$ install.packages("aphid", dependencies = F)
$ library("devtools")
$ install_github("roblanf/MixtureModelHMM", dependencies = F)
```

### Notes
- The code is developed and tested in macOS Monterey v12.6 and R v4.1.3 (installed in Conda environment)

### Dataset
<p>Edelman, N. B., Frandsen, P. B., Miyagi, M., Clavijo, B., Davey, J., Dikow, R. B., . . . Mallet, J. (2019). Genomic architecture and introgression shape a butterfly radiation. <i>Science, 366</i>(6465), 594-599. doi:<a href="https://doi.org/10.1126/science.aaw2090">10.1126/science.aaw2090</a></p>