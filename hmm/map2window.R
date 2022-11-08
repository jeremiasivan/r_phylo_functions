# this code is used to plot sliding tree and MAST tree topologies
# required packages: logger, optparse, stringr, utils, ggplot2

library(logger)
library(optparse)
library(stringr)
library(utils)
library(ggplot2)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-a","--analysis"), type="integer", default=1, 
              help="type of analysis: 1=all (default); 2=siteprob; 3=hmm", metavar="integer"),
  make_option("--siteprob", type="character", default=NULL, 
              help=".siteprob iqtree2 file (required for analysis 2)", metavar="character"),
  make_option("--annotation", type="character", default=NULL,
              help="annotation file from window_annotation.R", metavar="character"),
  make_option("--alninfo", type="character", default=NULL, 
              help=".alninfo iqtree2 file (required for analysis 2)", metavar="character"),
  make_option("--hmm", type="character", default=NULL, 
              help="hmm classification output (required for analysis 3)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (file.exists(opt$out)==F) {
  dir.create(opt$out)
}

plotSiteProb <- function(siteprobf, alninfof, annotationf, outdir) {
  if (is.null(siteprobf) || is.null(annotationf) || is.null(alninfof)){
    print_help(opt_parser)
    stop("Missing argument: .siteprob, .alninfo, or annotation files (see --help)", call.=FALSE)
  }
  
  # read input files
  log_info("Reading files...")
  alninfo <- read.table(alninfof, header = T)
  annotation <- read.csv(annotationf, sep = "\t")
  siteprob <- read.csv(siteprobf, sep = "\t")
  siteprob$Site <- NULL
  
  # create empty dataframe
  tree_header <- str_sort(colnames(siteprob))
  tree_topology <- str_sort(unique(annotation$topology))
  
  empty_matrix <- as.data.frame(matrix(rep(0, length(tree_header)), nrow=length(tree_topology), ncol=length(tree_header)+1))
  names(empty_matrix) <- c(tree_header,"NT")
  
  output <- data.frame(topology = tree_topology)
  output <- cbind(output, empty_matrix)
  
  # calculate the proportion of sites per window
  log_info("Calculating sites-supported windows...")
  pb = txtProgressBar(min = 0, max = nrow(annotation), initial = 0, style = 3) 
  
  for (i in 1:nrow(annotation)){
    subset <- siteprob[annotation$start[i]:annotation$stop[i],]
    inf_sites <- alninfo[annotation$start[i]:annotation$stop[i],]
    idx_topology <- which(output$topology == annotation$topology[i])
    
    for (j in 1:nrow(subset)) {
      order_lnl <- sort(subset[j,])
      low_bound <- as.numeric(order_lnl[ncol(subset)] - abs(0.1*order_lnl[ncol(subset)]))
      
      # exclude non-informative sites
      if (inf_sites$Stat[j] != "I"){
        output$NT[idx_topology] <- output$NT[idx_topology] + 1
      } else {
        idx_max <- which(subset[j,] > low_bound)

        # exclude sites that prefer all topologies
        if (length(idx_max) == ncol(subset)){
          output$NT[idx_topology] <- output$NT[idx_topology] + 1
        } else {
          output[idx_topology,idx_max+1] <- output[idx_topology,idx_max+1] + 1
        }
      }
    }
    
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  write.table(output, file=paste(outdir,"/mast_to_window.tsv",sep=""), sep = "\t", quote = F, row.names = F)
  
  # visualization
  log_info("Visualizing...")
  sub_output <- data.frame()
  for (i in 1:nrow(output)) {
    for (j in 3:ncol(output)-1) {
      sub_output <- rbind(sub_output, c(output$topology[i],colnames(output[j]),output[i,j]))
    }
  }
  colnames(sub_output) <- c("topology","sites","frequency")
  sub_output$sites <- factor(sub_output$sites,levels = unique(sub_output$sites),ordered = T)
  
  tiff(filename=paste(outdir,"/sites_ditribution.tiff",sep=""), units="px", width=2880, height=1800)
  print(ggplot(sub_output, aes(fill=sites, y=as.numeric(frequency), x=topology)) + 
          geom_bar(position="dodge", stat="identity") +
          ggtitle("Distribution of Informative Sites in Sliding Windows Topology") + 
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
}

plotHMMClassification <- function(hmmf, annotationf, outdir){
  if (is.null(annotationf) || is.null(hmmf)){
    print_help(opt_parser)
    stop("Missing argument: hmm or annotation files (see --help)", call.=FALSE)
  }
  
  # read input files
  log_info("Reading files...")
  annotation <- read.csv(annotationf, sep = "\t")
  hmm <- read.table(hmmf)
  
  # create empty dataframe
  tree_header <- str_sort(unique(hmm[,1]))
  tree_topology <- str_sort(unique(annotation$topology))
  
  empty_matrix <- as.data.frame(matrix(rep(0, length(tree_header)), nrow=length(tree_topology), ncol=length(tree_header)))
  names(empty_matrix) <- tree_header
  
  output <- data.frame(topology = tree_topology)
  output <- cbind(output, empty_matrix)
  
  # calculate the proportion of sites per window
  log_info("Calculating sites-supported windows...")
  pb = txtProgressBar(min = 0, max = nrow(annotation), initial = 0, style = 3) 
  
  for (i in 1:nrow(annotation)){
    subset <- hmm[annotation$start[i]:annotation$stop[i],]
    idx_topology <- which(output$topology == annotation$topology[i])
    
    for (j in 1:length(subset)) {
      idx_max <- match(subset[j],tree_header)
      output[idx_topology,idx_max+1] <- output[idx_topology,idx_max+1] + 1
    }
    
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  write.table(output, file=paste(outdir,"/hmm_to_window.tsv",sep=""), sep = "\t", quote = F, row.names = F)
  
  # visualization
  log_info("Visualizing...")
  sub_output <- data.frame()
  for (i in 1:nrow(output)) {
    for (j in 2:ncol(output)) {
      sub_output <- rbind(sub_output, c(output$topology[i],colnames(output[j]),output[i,j]))
    }
  }
  colnames(sub_output) <- c("topology","sites","frequency")
  sub_output$sites <- factor(sub_output$sites,levels = unique(sub_output$sites),ordered = T)
  
  tiff(filename=paste(outdir,"/hmm_ditribution.tiff",sep=""), units="px", width=2880, height=1800)
  print(ggplot(sub_output, aes(fill=sites, y=as.numeric(frequency), x=topology)) + 
          geom_bar(position="dodge", stat="identity") +
          ggtitle("Distribution of HMM Sites in Sliding Windows Topology") + 
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
}

runAll <- function(siteprob, alninfo, annotation, hmm, outdir){
  if (is.null(siteprob) || is.null(annotation) || is.null(alninfo) || is.null(hmm)){
    print_help(opt_parser)
    stop("Missing argument: .siteprob, .alninfo, hmm, or annotation files (see --help)", call.=FALSE)
  }
  
  plotSiteProb(siteprob,alninfo,annotation,outdir)
  plotHMMClassification(hmm,annotation,outdir)
}

# switch according to the type of analysis
switch(opt$analysis,
       runAll(opt$siteprob, opt$alninfo, opt$annotation, opt$hmm, opt$out),
       plotSiteProb(opt$siteprob, opt$alninfo, opt$annotation, opt$out),
       plotHMMClassification(opt$hmm, opt$annotation, opt$out)
)

log_info("Done!")
