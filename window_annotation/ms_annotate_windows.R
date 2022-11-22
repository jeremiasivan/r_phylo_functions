# this code is used to summarize trees from ms (https://uchicago.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13) and compare the results with other techniques
# required packages: optparse, logger, dplyr, ape, data.table, ggplot2, forcats
library(optparse)
library(logger)
library(dplyr)
library(ape)
library(data.table)
library(ggplot2)
library(forcats)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-a", "--analysis"), type="integer", default=1,
              help="type of analysis: 1=ms trees; 2=ms summary; 3=ms vs sliding-window; 4=ms vs hmm", metavar="integer"),
  make_option("--ms", type="character", default=NULL,
              help="ms file", metavar="character"),
  make_option("--mssum", type="character", default=NULL,
              help="ms summary file", metavar="character"),
  make_option("--swsum", type="character", default=NULL,
              help="sliding-window summary file", metavar="character"),
  make_option("--hmm", type="character", default=NULL,
              help="hmm file", metavar="character"),
  make_option("--treemap", type="character", default=NULL,
              help="map of tab-separated rooted and unrooted trees", metavar="character"),            
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# set output directory
outdir <- paste(opt$out,"/ms/",sep="")
if (file.exists(outdir)==F) {
  dir.create(outdir, recursive = T)
}

# ms summary
msSummary <- function(ms, outdir) {
  if (is.null(ms)){
    stop("Missing argument: ms file (see --help)", call.=FALSE)
  }
  
  # summarize ms (tree) output
  log_info("Reading input files...")
  msfile <- read.delim(ms)
  mstrees <- subset(msfile, grepl("^\\[", msfile[,1]))
  
  trees <- data.frame(do.call('rbind', strsplit(as.character(mstrees[,1]),'\\[|\\]')))
  trees[,1] <- NULL
  colnames(trees) <- c("length","tree")
  trees$length <- as.numeric(trees$length)
  
  # summarize trees
  summary <- trees %>%
    group_by(tree,
             group_run = data.table::rleid(tree)) %>%
    summarise_all(sum) %>%
    arrange(group_run)
  summary$group_run <- NULL
  
  # extract start, stop, and tree topology
  log_info("Annotating...")
  start <- c(1)
  stop <- c(summary$length[1])
  
  first_tree <- read.tree(text=summary$tree[1])
  first_tree$edge.length <- NULL
  topology <- c(write.tree(first_tree))
  
  for (i in 2:nrow(summary)){
    start <- c(start,stop[length(stop)]+1)
    stop <- c(stop,stop[length(stop)]+summary$length[i])
    
    temp_tree <- read.tree(text=summary$tree[i])
    temp_tree$edge.length <- NULL
    topology <- c(topology,write.tree(temp_tree))
  }
  
  summary <- cbind(summary,start=start,stop=stop,topology=topology)
  write.table(summary,file=paste(outdir,"tree_summary.txt",sep=""), sep="\t", quote=F, row.names=F)
  
  # summarize topologies
  t_summary <- summary %>%
    group_by(topology,
             group_run = data.table::rleid(topology)) %>%
    summarise(length=sum(length)) %>%
    arrange(group_run)
  t_summary$group_run <- NULL
  
  # extract start, stop, and tree topology
  start <- c(1)
  stop <- c(t_summary$length[1])
  
  for (i in 2:nrow(t_summary)){
    start <- c(start,stop[length(stop)]+1)
    stop <- c(stop,stop[length(stop)]+t_summary$length[i])
  }
  
  t_summary <- cbind(t_summary,start=start,stop=stop)
  write.table(t_summary,file=paste(outdir,"topology_summary.txt",sep=""), sep="\t", quote=F, row.names=F)
  
  # summarize unique topology
  write.table(data.frame(unique(t_summary$topology)),file=paste(outdir,"unique_topology.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  log_info("Done!")
}

# extract treefile for Seq-Gen
msTrees <- function(ms, outdir){
  if (is.null(ms)){
    stop("Missing argument: ms file (see --help)", call.=FALSE)
  }
  
  # summarize ms (tree) output
  msfile <- read.delim(ms)
  mstrees <- subset(msfile, grepl("^\\[", msfile[,1]))
  
  write.table(mstrees,file=paste(outdir,"ms_trees.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  log_info("Done!")
}

# simulation vs sliding-window analysis
msSWComparison <- function(ms,sw,treeMap,outdir){
  if (is.null(ms) || is.null(sw) || is.null(treeMap)){
    stop("Missing argument: ms, sliding-window, or treemap file (see --help)", call.=FALSE)
  }

  log_info("Reading input files...")
  msfile <- read.table(ms)
  colnames(msfile) <- msfile[1,]
  msfile <- msfile[-1,]
  
  swfile <- read.table(sw)
  colnames(swfile) <- swfile[1,]
  swfile <- swfile[-1,]
  
  tmfile <- read.table(treeMap)
  for (i in 1:nrow(msfile)) {
    idx <- which(tmfile[,2]==msfile$topology[i])
    msfile$topology[i] <- tmfile[idx,1]
  }
  
  # extracting essential columns
  dat <- apply(swfile, 1,
               function(x) data.table(bp = x["start"]:x["stop"], newick = x["newick"], topology = x["topology"])) %>%
    rbindlist()
  
  # create empty dataframe
  tree_header <- tmfile[,1]
  tree_topology <- tree_header
  
  empty_matrix <- as.data.frame(matrix(rep(0, length(tree_header)+1), nrow=length(tree_topology), ncol=length(tree_header)+1))
  names(empty_matrix) <- c(tree_header,"NT")
  
  output <- data.frame(topology = tree_topology)
  output <- cbind(output, empty_matrix)
  
  # annotating
  log_info("Annotating...")
  for (i in 1:nrow(msfile)){
    sites <- dat[msfile$start[i]:msfile$stop[i],]
    topology_idx <- which(output$topology==msfile$topology[i])
    
    sites_sum <- sites %>%
      group_by(newick) %>%
      summarise(n=n())
    
    for (j in 1:nrow(sites_sum)) {
      window_idx <- which(colnames(output)==sites_sum$newick[j])
      
      if (length(window_idx) == 0) {
        output$NT[topology_idx] = output$NT[topology_idx] + sites_sum$n[j]
      } else {
        output[topology_idx,window_idx] = output[topology_idx,window_idx] + sites_sum$n[j]
      }
    }
  }
  
  write.table(output, file=paste(outdir,"/window_to_simulation.tsv",sep=""), sep = "\t", quote = F, row.names = F)
  
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
          aes(x = fct_inorder(topology)) +
          ggtitle("Distribution of Sites from Sliding-Window Analysis in Simulated Dataset") + 
          xlab("True topology from simulated data") +
          ylab("Number of sites from sliding-window analysis") +
          theme(
            plot.title = element_text(hjust = 0.5, size = 50),
            axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
            axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
            axis.text.y=element_text(size=30),
            axis.text.x=element_text(size=30),
            legend.title=element_blank(),
            legend.text=element_text(size=30),
            legend.key.size=unit(2,"cm")
          ))
  dev.off()
  log_info("Done!")
}

msMastComparison <- function(ms,hmm,treeMap,outdir){
  if (is.null(ms) || is.null(hmm) || is.null(treeMap)){
    stop("Missing argument: ms, hmm, or treemap file (see --help)", call.=FALSE)
  }

  log_info("Reading input files...")
  msfile <- read.table(ms)
  colnames(msfile) <- msfile[1,]
  msfile <- msfile[-1,]
  
  hmmfile <- readLines(hmm)
  hmmtrees <- subset(hmmfile, grepl("^\\[", hmmfile))
  
  trees <- data.frame(do.call('rbind', strsplit(as.character(hmmtrees),'\\[|,|\\]')))
  trees[,1] <- NULL
  colnames(trees) <- c("start","stop","tree")
  
  tmfile <- read.table(treeMap)
  for (i in 1:nrow(msfile)) {
    idx <- which(tmfile[,2]==msfile$topology[i])
    msfile$topology[i] <- tmfile[idx,1]
  }
  
  # extracting essential columns
  dat <- apply(trees, 1,
               function(x) data.table(bp = x["start"]:x["stop"], tree = x["tree"])) %>%
    rbindlist()
  
  # create empty dataframe
  tree_header <- tmfile[,1]
  tree_topology <- tree_header
  
  empty_matrix <- as.data.frame(matrix(rep(0, length(tree_header)+1), nrow=length(tree_topology), ncol=length(tree_header)+1))
  names(empty_matrix) <- c(tree_header,"NT")
  
  output <- data.frame(topology = tree_topology)
  output <- cbind(output, empty_matrix)
  
  # annotating
  log_info("Annotating...")
  for (i in 1:nrow(msfile)){
    sites <- dat[msfile$start[i]:msfile$stop[i],]
    topology_idx <- which(output$topology==msfile$topology[i])
    
    sites_sum <- sites %>%
      group_by(tree) %>%
      summarise(n=n())
    
    for (j in 1:nrow(sites_sum)) {
      window_idx <- as.numeric(sites_sum$tree[j]) + 1
      output[topology_idx,window_idx] = output[topology_idx,window_idx] + sites_sum$n[j]
    }
  }
  
  write.table(output, file=paste(outdir,"/hmm_to_simulation.tsv",sep=""), sep = "\t", quote = F, row.names = F)
  
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
  
  tiff(filename=paste(outdir,"/hmm_sites_ditribution.tiff",sep=""), units="px", width=2880, height=1800)
  print(ggplot(sub_output, aes(fill=sites, y=as.numeric(frequency), x=topology)) + 
          geom_bar(position="dodge", stat="identity") +
          aes(x = fct_inorder(topology)) +
          ggtitle("Distribution of Sites from MAST+HMM in Simulated Dataset") + 
          xlab("True topology from simulated data") +
          ylab("Number of sites from MAST+HMM") +
          theme(
            plot.title = element_text(hjust = 0.5, size = 50),
            axis.title.x=element_text(size = 40, margin = margin(t=20, r=0, b=0, l=0)),
            axis.title.y=element_text(size = 40, margin = margin(t=0, r=20, b=0, l=0)),
            axis.text.y=element_text(size=30),
            axis.text.x=element_text(size=30),
            legend.title=element_blank(),
            legend.text=element_text(size=30),
            legend.key.size=unit(2,"cm")
          ))
  dev.off()
  log_info("Done!")
}

# switch according to the type of analysis
switch(opt$analysis,
       msTrees(opt$ms,outdir),
       msSummary(opt$ms,outdir),
       msSWComparison(opt$mssum,opt$swsum,opt$treemap,outdir),
       msMastComparison(opt$mssum,opt$hmm,opt$treemap,outdir),
)

