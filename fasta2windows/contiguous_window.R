# this code is used to combine windows with the same topology as contiguous regions and build window trees afterwards

library(logger)
library(optparse)
library(seqinr)
library(dplyr)
library(data.table)
library(ggplot2)
library(forcats)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-f","--fasta"), type="character", default=NULL, 
              help="concatenated fasta alignment", metavar="character"),
  make_option("--swsum", type="character", default=NULL, 
              help="annotation of non-overlapping windows", metavar="character"),
  make_option("--iqtree2", type="character", default=NULL, 
              help="iqtree2 directory", metavar="character"),
  make_option("--outgroup", type="character", default=NULL, 
              help="outgroup for iqtree2 run (optional)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$fasta) || is.null(opt$swsum) || is.null(opt$iqtree2)){
  print_help(opt_parser)
  stop("Missing argument (see --help)", call.=FALSE)
}

outdir <- paste(opt$out,"/summary/contiguous/",sep="")
if (file.exists(outdir)==F) {
  dir.create(outdir, recursive = T)
}

# read files
log_info("Reading Non-overlapping Windows annotation file...")
swfile <- read.table(opt$swsum)
colnames(swfile) <- swfile[1,]
swfile <- swfile[-1,]

log_info("Counting contiguous blocks...")
count_contiguous <- swfile %>% 
  group_by(topology,
           group_run = data.table::rleid(topology)) %>% 
  summarise(count = n()) %>% 
  arrange(group_run)
count_contiguous$group_run <- NULL

write.table(count_contiguous, file=paste(outdir,"summary_count.txt",sep=""), quote=F, row.names=F)

cons_windows <- count_contiguous %>%
  group_by(topology, count) %>%
  summarise(total = n())

tiff(filename=paste(outdir,"cons_windows.tiff",sep=""), units="px", width=2880, height=1800)
print(ggplot(cons_windows, aes(fill=topology)) + 
        geom_bar(position="dodge", stat="identity", aes(y=total, x=count)) +
        aes(x = fct_inorder(topology)) +
        ggtitle("Distribution of Consecutive Windows Per Topology") + 
        # xlab("Length of consecutiveness") +
        # ylab("Number of windows") +
        facet_grid(~topology) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 50),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=30),
          axis.text.x=element_text(size=30),
          legend.title=element_blank(),
          legend.text=element_text(size=30),
          legend.key.size=unit(2,"cm")
        ))
dev.off()

log_info("Reading FASTA file...")
s <- read.fasta(opt$fasta, whole.header = T)
sdf <- do.call(rbind,s)

unique_topology <- sort(unique(swfile$topology))

log_info("Generating topology-based contiguous sequence...")
for (i in unique_topology) {
  sub_sw <- subset(swfile[swfile$topology == i,])
  align_idx <- unlist(Map(seq, sub_sw$start, sub_sw$stop))
    
  alignment <- sdf[,align_idx]
  subfasta <- setNames(split(alignment, seq(nrow(alignment))), rownames(alignment))
  
  write.fasta(sequences=subfasta, names=names(subfasta), file.out=paste(outdir,i,".fa",sep=""), nbchar = 80)
}

log_info("Generating contiguous window trees...")
treedir <- paste(outdir,"trees/",sep="")
if (file.exists(treedir)==F) {
  dir.create(treedir, recursive = T)
}

prefix <- paste(treedir,"contiguous", sep="")
if (is.null(opt$outgroup)){
  system(paste(opt$iqtree2,"-S",outdir,"--quiet -redo -pre",prefix,sep=" "))
} else {
  system(paste(opt$iqtree2,"-S",outdir,"--quiet -redo -pre",prefix,"-o",opt$outgroup,sep=" "))
}

log_info("Done!")
