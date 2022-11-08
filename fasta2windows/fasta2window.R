# this code is used to generate window alignments from concatenated fasta
# required packages: logger, optparse, seqinr

library(logger)
library(optparse)
library(seqinr)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-f","--fasta"), type="character", default=NULL, 
              help="concatenated fasta alignment", metavar="character"),
  make_option(c("-s","--size"), type="integer", default=NULL, 
              help="window size in bp", metavar="integer"),
  make_option(c("-n","--name"), type="character", default="chr", 
              help="chromosome name [current=%default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$fasta) || is.null(opt$size)){
  print_help(opt_parser)
  stop("Missing argument: fasta or window size (see --help)", call.=FALSE)
}

# set window size
wsize <- as.numeric(opt$size)
if (is.na(wsize) || wsize == 0) {
  stop("Invalid window size (see --help)", call.=FALSE)
}

# set output directory
outdir <- paste(opt$out,"/windows/",sep="")
if (file.exists(outdir)==F) {
  dir.create(outdir, recursive = T)
}

log_info("Reading FASTA file...")
s <- read.fasta(opt$fasta, whole.header = T)

log_info("Generating non-overlapping windows...")
start <- 1
pb = txtProgressBar(min = 0, max = floor(getLength(s)[1]/wsize), initial = 0, style = 3)

for (i in 1:floor(getLength(s)[1]/wsize)) {
  subfasta <- lapply(s, function(x) x[seq(from = start, to = as.numeric(i*wsize))])
  df <- do.call(rbind,subfasta)

  # remove sites that are not present in each individual
  df <- df[, !colSums(df == "*")]
  
  if (is.null(ncol(df))) {
    start <- start + wsize
    setTxtProgressBar(pb,i)
    next
  }
  
  # count number of fully-aligned sites
  numInfSites <- df[, !colSums(df == "n")]
  if (is.null(ncol(numInfSites))) {
    start <- start + wsize
    setTxtProgressBar(pb,i)
    next
  }
  
  numInfSites <- numInfSites[, !colSums(numInfSites == "-")]
  if (is.null(ncol(numInfSites)) || ncol(numInfSites) < as.numeric(0.1*wsize)) {
    start <- start + wsize
    setTxtProgressBar(pb,i)
    next
  }
  
  subfasta <- setNames(split(df, seq(nrow(df))), rownames(df))
  write.fasta(sequences=subfasta, names=names(subfasta),file.out=paste(outdir,"/",opt$name,"_",i,".fa",sep=""), nbchar = 80)
  
  start <- start + wsize
  setTxtProgressBar(pb,i)
}

close(pb)

log_info("Done!")