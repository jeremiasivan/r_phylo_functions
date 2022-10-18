# this code contains list of functions for gene_tree_annotation.R
# required packages: logger, ape, hash, seqinr, dplyr, data.table, ggplot2, RColorBrewer
# required software:
# - concatenation: seqkit (https://github.com/shenwei356/seqkit)
# - gene trees: iqtree2 (https://github.com/iqtree/iqtree2)
# - statistics: n/a (but it requires tree topologies, which you can get using PAML)

library(logger)
library(ape)
library(hash)
library(seqinr)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)

# function to concatenate fasta sequences
concatFasta <- function(inputdir, outdir, seqkitdir, nthread, chr){
  if (is.null(inputdir)){
    stop("Missing argument: fasta files directory (see --help)", call.=FALSE)
  }
  
  if (is.null(seqkitdir)){
    stop("Missing argument: seqkit executable (see --help)", call.=FALSE)
  }
  
  outdir <- paste(outdir,"/concatenation",sep="")
  if (file.exists(outdir)==F) {
    dir.create(outdir)
  }
  
  log_info("Running concatenation steps...")
  listdir <- paste(outdir,"/list_fasta.txt", sep="")
  fastas <- setNames(data.frame(list.files(inputdir)), "file")
  fastas$file <- paste(inputdir,"/",fastas$file, sep="")
  write.table(fastas, file=listdir, row.names=F, col.names=F, quote=F)
  
  concatdir <- paste(outdir,"/",chr,"_concat.fa", sep="")
  system(paste(seqkitdir,"concat -f --quiet --infile-list",listdir,"-j",nthread,"-o",concatdir, sep=" "))
  
  log_info("Done!")
}

# function to generate gene/window trees
geneTree <- function(inputdir, outdir, iqtreedir, nthread, chr){
  if (is.null(inputdir)){
    stop("Missing argument: fasta files directory (see --help)", call.=FALSE)
  }
  
  if (is.null(iqtreedir)){
    stop("Missing argument: iqtree2 executable (see --help)", call.=FALSE)
  }
  
  log_info("Generating gene trees...")
  outdir <- paste(outdir,"/trees/",sep="")
  if (file.exists(outdir)==F) {
    dir.create(outdir)
  }
  
  prefix <- paste(outdir,chr, sep="")
  system(paste(iqtreedir,"-S",inputdir,"--quiet -redo -pre",prefix,"-T",nthread, sep=" "))
  
  log_info("Done!")
}

# function for annotation
geneAnnotate <- function(outdir, topologydir, chr) {
  if (is.null(topologydir)){
    stop("Missing argument: tree topology directory (see --help)", call.=FALSE)
  }
  
  log_info("Reading topology file...")
  t_hash <- hash()
  t_count <- 1
  topology_file <- file(description=topologydir, open="r", blocking = TRUE)
  repeat {
    tl <- readLines(topology_file, n=1, warn=F)
    if (identical(tl, character(0))) {
      close(topology_file)
      break
    }
    t_hash[tl] = paste("T",t_count,sep="")
    t_count <- t_count+1
  }
  
  seq_list <- file(description = paste(outdir,"/concatenation/list_fasta.txt",sep=""), open="r", blocking = TRUE)
  iqtree_file <- file(description = paste(outdir,"/trees/",chr,".treefile",sep=""), open="r", blocking = TRUE)
  
  seq_len <- data.frame()
  start <- 0
  
  log_info("Annotating...")
  repeat {
    pl <- readLines(seq_list, n=1)
    if (identical(pl, character(0))) {
      colnames(seq_len) <- c("source","chromosome","start","stop","length","newick","topology")
      closeAllConnections()
      break
    }
    
    ttop <- "NT"
    tl <- read.tree(text=readLines(iqtree_file, n=1))
    tl$edge.length <- NULL
    tree <- write.tree(tl)
    if (has.key(tree, t_hash) == TRUE) {
      ttop <- as.character(values(t_hash[tree]))
    } 
    
    s <- read.fasta(pl, as.string = T)
    seq_len <- rbind(seq_len, c(pl,chr,start+1,start+getLength(s)[1],getLength(s)[1],tree, ttop))
    
    start <- start+getLength(s)[1]
  }
  
  outdir <- paste(outdir,"/summary",sep="")
  if (file.exists(outdir)==F) {
    dir.create(outdir)
  }
  
  write.table(seq_len, paste(outdir,"/",chr,"_annotation.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
  log_info("Done!")
}

# function for topology-based chromosome visualization
geneChr <- function(outdir, chr) {
  outdir <- paste(outdir,"/summary/",sep="")
  seq_annot <- read.csv(paste(outdir,chr,"_annotation.txt",sep=""),sep="\t")
  
  log_info("Reading annotation data...")
  dat <- apply(seq_annot, 1,
               function(x) data.table(chromosome = x["chromosome"], bp = x["start"]:x["stop"], topology = x["topology"])) %>%
    rbindlist()
  
  log_info("Plotting...")
  tiff(filename=paste(outdir,chr,"_topology.tiff",sep=""),units="px",width=250,height=1800)
  print(ggplot() + 
    geom_rect(data=dat,aes(fill=topology,ymin=bp-1,ymax=bp,xmin=0,xmax=1)) +
    labs(x=chr) +
    # scale_fill_manual(values = c("black","firebrick2","dodgerblue3","darkseagreen1","darkgoldenrod1","rosybrown1")) +
    scale_fill_brewer(palette="Set1") +
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid")) +
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_text(size=20),
      legend.title=element_blank(),
      legend.text=element_text(size=20),
      legend.key.size=unit(2,"cm")
      ))
  dev.off()
  log_info("Done!")
}

# function to do all
allSteps <- function(inputdir, outdir, seqkitdir, iqtreedir, topologydir, nthread, chr) {
  if (is.null(inputdir)){
    stop("Missing argument: fasta files directory (see --help)", call.=FALSE)
  }
  
  if (is.null(seqkitdir)){
    stop("Missing argument: seqkit executable (see --help)", call.=FALSE)
  }
  
  if (is.null(iqtreedir)){
    stop("Missing argument: iqtree2 executable (see --help)", call.=FALSE)
  }
  
  if (is.null(topologydir)){
    stop("Missing argument: tree topology directory (see --help)", call.=FALSE)
  }
  
  concatFasta(inputdir, outdir, seqkitdir, nthread, chr)
  geneTree(inputdir, outdir, iqtreedir, nthread, chr)
  geneAnnotate(outdir, topologydir, chr)
  geneChr(outdir, chr)
}