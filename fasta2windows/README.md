## Generate Window Alignments
The following R script can be used to construct fixed-length window alignments from concatenated FASTA.

### Notes
- During development, the concatenated FASTA comes from multiple MAF blocks which is converted using PHAST (*<a href="http://compgen.cshl.edu/phast/">msa_view</a>*).
- Default filters: (i) DNA regions must present in all species; (ii) fully-aligned regions (i.e. no gaps or *N*) must be >10% of the window size
- The code is developed and tested in macOS Monterey v12.6 and R v4.1.3 (installed in Conda environment)