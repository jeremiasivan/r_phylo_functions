# this code is used to generate concatenated alignment, window trees, and annotation
# required packages: optparse, tidyverse, logger, gtools, ape, hash, seqinr, dplyr, data.table, ggplot2, RColorBrewer
# required software:
# - concatenation: seqkit (https://github.com/shenwei356/seqkit)
# - window trees: iqtree2 (https://github.com/iqtree/iqtree2)

library(optparse)
library(tidyverse)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-a", "--analysis"), type="integer", default=1,
              help="type of analysis: 1=all (default); 2=concatenation; 3=gene trees; 4=annotation", metavar="integer"),
  make_option("--chr", type="character", default="chr", 
              help="chromosome name [default=%default]", metavar="character"),
  make_option("--fasta", type="character", default=NULL, 
              help="directory of fasta files (required for analysis 2,3)", metavar="character"),
  make_option("--outgroup", type="character", default=NULL, 
              help="outgroup species (optional for analysis 3)", metavar="character"),
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

# check current file location (source: https://stackoverflow.com/a/55322344)
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
    if (length(this_file)==0)
    {
      this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

# import functions
source(paste(getCurrentFileLocation(),"/annotate_windows_f.R",sep=""))

# switch according to the type of analysis
switch(opt$analysis,
       allSteps(opt$fasta, opt$out, opt$seqkit, opt$iqtree2, opt$topology, opt$outgroup, opt$thread, opt$chr),
       concatFasta(opt$fasta, opt$out, opt$seqkit, opt$thread, opt$chr),
       windowTree(opt$fasta, opt$out, opt$iqtree2, opt$outgroup, opt$thread, opt$chr),
       windowAnnotate(opt$out, opt$topology, opt$chr)
)
