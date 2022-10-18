# this code is used to visualize phylogenetic trees with the corresponding branch lengths
# image resolutions and fonts can be changed in line 42-44
# required packages: ape, optparse, logger

library(logger)
library(optparse)
library(ape)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input treefile", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input treefile)", call.=FALSE)
}

if (file.exists(opt$out)==F) {
  dir.create(opt$out)
}

# read the treefile
tre <- read.tree(opt$file)
log_info(paste("Number of trees: ",as.character(length(tre))),sep="")

# plot and save each tree topology
count <- 1
for (tree in 1:length(tre)){
  tiff(filename=paste(opt$out,"/T",as.character(count),".tiff",sep=""), units="px", width=2880, height=1800)
  plot.phylo(tre[[count]], cex=2.4, font=3)
  edgelabels(round(tre[[count]]$edge.length,7), frame="none", cex=2, adj=c(0.5,-0.5))
  dev.off()
  
  count <- count + 1
}

log_info("Done!")
