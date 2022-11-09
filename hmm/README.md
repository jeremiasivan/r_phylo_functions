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

### Known Error(s)
- `run_hmm.R` might not work due to an issue in `save_report` function in <a href="https://github.com/roblanf/MixtureModelHMM/issues/20">`HMMMixtureModel`</a>. As it only happens in Linux OS, you might want to use another OS or comment out the code if not necessary.

### Notes
- By default, sites with <10% likelihood difference between topologies are considered to support multiple topologies
- The code is developed and tested in macOS Monterey v12.6 and R v4.1.3 (installed in Conda environment)