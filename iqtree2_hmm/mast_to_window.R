# this code is used to plot sliding tree and MAST tree topologies
# required packages: logger, optparse, utils, ggplot2

library(logger)
library(optparse)
library(utils)
library(ggplot2)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option("--siteprob", type="character", default=NULL, 
              help=".siteprob iqtree2 file", metavar="character"),
  make_option("--annotation", type="character", default=NULL,
              help="annotation file from window_annotation.R", metavar="character"),
  make_option("--alninfo", type="character", default=NULL, 
              help=".alninfo iqtree2 file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$siteprob) || is.null(opt$annotation) || is.null(opt$alninfo)){
  print_help(opt_parser)
  stop("Missing argument: .siteprob, .alninfo, or annotation files (see --help)", call.=FALSE)
}

if (file.exists(opt$out)==F) {
  dir.create(opt$out)
}

# read input files
log_info("Reading files...")
alninfo <- read.table(opt$alninfo, header = T)
annotation <- read.csv(opt$annotation, sep = "\t")
siteprob <- read.csv(opt$siteprob, sep = "\t")
siteprob$Site <- NULL

# create empty dataframe
tree_header <- paste("T", 1:ncol(siteprob), sep="")
empty_matrix <- as.data.frame(matrix(rep(0, 2+length(tree_header)), nrow=ncol(siteprob)+1, ncol=ncol(siteprob)+2))
names(empty_matrix) <- c("NT",tree_header,"match")

output <- data.frame(topology = c("NT",tree_header))
output <- cbind(output, empty_matrix)

# calculate the proportion of sites per window
log_info("Calculating sites-supported windows...")
pb = txtProgressBar(min = 0, max = nrow(annotation), initial = 0, style = 3) 

for (i in 1:nrow(annotation)){
  subset <- siteprob[annotation$start[i]:annotation$stop[i],]
  inf_sites <- alninfo[annotation$start[i]:annotation$stop[i],]
  idx_topology <- which(output$topology == annotation$topology[i])
  
  for (j in 1:nrow(subset)) {
    if (inf_sites$Stat[j] != "I"){
      output$NT[idx_topology] <- output$NT[idx_topology] + 1
    } else {
      idx_max <- which(subset[j,] == max(subset[j,]))
      if (length(idx_max) != 1){
        cat("\n")
        log_info(paste("Mutiple best topology for window",i,"site",j,sep=" "))
      }
      output[idx_topology,idx_max+2] <- output[idx_topology,idx_max+2] + 1
    }
  }
  
  setTxtProgressBar(pb,i)
}

close(pb)

# calculate the proportion of matching sites per topology
log_info("Calculating matching sites proportion...")
for (i in 1:nrow(output)) {
  output$match[i] <- round((output[i,i+1]/sum(output[i,3:ncol(output)-1])*100),5)
}

write.table(output, file=paste(opt$out,"/mast_to_window.tsv",sep=""), sep = "\t", quote = F, row.names = F)

# visualization
sub_output <- data.frame()
for (i in 1:nrow(output)) {
  for (j in 4:ncol(output)-1) {
    sub_output <- rbind(sub_output, c(output$topology[i],colnames(output[j]),output[i,j]))
  }
}
colnames(sub_output) <- c("topology","sites","frequency")
sub_output$sites <- factor(sub_output$sites,levels = unique(sub_output$sites),ordered = T)

tiff(filename=paste(opt$out,"/sites_ditribution.tiff",sep=""), units="px", width=2880, height=1800)
print(ggplot(sub_output, aes(fill=sites, y=as.numeric(frequency), x=topology)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Distribution of Informative Sites in Sliding Windows Topology") + 
  xlab("sliding windows topology") +
  ylab("frequency of sites") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_text(size=30),
    axis.text.x=element_text(size=30),
    legend.title=element_blank(),
    legend.text=element_text(size=30),
    legend.key.size=unit(2,"cm")
  ))
dev.off()

log_info("Done!")
