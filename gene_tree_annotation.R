# this code is used to generate concatenated alignment, gene/window trees, and statistics
# required packages: optparse, logger, ape, hash, seqinr, dplyr, data.table, ggplot2, RColorBrewer
# required software:
# - concatenation: seqkit (https://github.com/shenwei356/seqkit)
# - gene trees: iqtree2 (https://github.com/iqtree/iqtree2)
# - statistics: n/a (but it requires tree topologies, which you can get using PAML)

library(optparse)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-a", "--analysis"), type="integer", default=1,
              help="type of analysis: 1=all (default); 2=concatenation; 3=gene trees; 4=annotation 5=visualization", metavar="integer"),
  make_option("--chr", type="character", default="chr", 
              help="chromosome name [default=%default]", metavar="character"),
  make_option("--fasta", type="character", default=NULL, 
              help="directory of fasta files (required for analysis 2,3)", metavar="character"),
  make_option("--topology", type="character", default=NULL, 
              help="directory of topology file (required for analysis 4)", metavar="character"),
  make_option("--seqkit", type="character", default=NULL, 
              help="seqkit executable (required for analysis 2)", metavar="character"),
  make_option("--iqtree2", type="character", default=NULL, 
              help="iqtree2 executable (required for analysis 3)", metavar="character"),
  make_option(c("-t", "--thread"), type="integer", default=4, 
              help="number of thread [default=%default]", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (file.exists(opt$out)==F) {
  dir.create(opt$out)
}

# import functions
source("gene_tree_annotation_functions.R")

# switch according to the type of analysis
switch(opt$analysis,
       allSteps(opt$fasta, opt$out, opt$seqkit, opt$iqtree2, opt$topology, opt$thread, opt$chr),
       concatFasta(opt$fasta, opt$out, opt$seqkit, opt$thread, opt$chr),
       geneTree(opt$fasta, opt$out, opt$iqtree2, opt$thread, opt$chr),
       geneAnnotate(opt$out, opt$topology, opt$chr),
       geneChr(opt$out, opt$chr)
)
