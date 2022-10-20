# this code is used to plot sliding tree and MAST tree topologies
# required packages: logger, optparse

library(logger)
library(optparse)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option("--siteprob", type="character", default=NULL, 
              help=".siteprob iqtree2 file", metavar="character"),
  make_option("--annotation", type="character", default=NULL,
              help="annotation file from window_annotation.R", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$siteprob) || is.null(opt$annotation)){
  print_help(opt_parser)
  stop("Missing argument: .siteprob or annotation files (see --help)", call.=FALSE)
}

if (file.exists(opt$out)==F) {
  dir.create(opt$out)
}

# read input files
log_info("Reading files...")
annotation <- read.csv(opt$annotation, sep = "\t")
siteprob <- read.csv(opt$siteprob, sep = "\t")
siteprob$Site <- NULL

# create empty dataframe
tree_header <- paste("T", 1:ncol(siteprob), sep="")
empty_matrix <- as.data.frame(matrix(rep(0, 2+length(tree_header)), nrow=ncol(siteprob)+1, ncol=ncol(siteprob)+2))
names(empty_matrix) <- c("NT",tree_header,"match (%)")

output <- data.frame(topology = c("NT",tree_header))
output <- cbind(output, empty_matrix)

# calculate the proportion of sites per window
log_info("Calculating sites-supported windows...")
for (i in 1:nrow(annotation)){
  subset <- siteprob[annotation$start[i]:annotation$stop[i],]
  idx_topology <- which(output$topology == annotation$topology[i])
  
  for (j in 1:nrow(subset)) {
    if (max(subset[j,])-min(subset[j,]) < 0.01){
      output$NT[idx_topology] = output$NT[idx_topology] + 1
    
    } else {
      idx_max <- which(subset[j,] == max(subset[j,]))
      output[idx_topology,idx_max+2] = output[idx_topology,idx_max+2] + 1
    }
  }
}

# calculate the proportion of matching sites per topology
for (i in 1:nrow(output)) {
  output$`match (%)`[i] <- round(output[i,i+1]/sum(output[i,3:ncol(output)-1])*100,5)
}

# visualization

write.csv(paste(opt$out,"/mast_to_window.csv",sep=""), sep = "\t", quote = F)

log_info("Done!")
