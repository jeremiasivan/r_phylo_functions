# this code is used to set up environment and run HMM using Conda
# required packages: logger, optparse, MixtureModelHMM (https://github.com/roblanf/MixtureModelHMM)

library(logger)
library(optparse)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="directory of iqtree2 files", metavar="character"),
  make_option("--conda", type="logical", default=F,
              help="set up conda environment; please make sure conda is installed", metavar="logical"),
  make_option(c("-o", "--out"), type="character", default=c_dir, 
              help=paste("output dir [current=",c_dir,"]", sep=""), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Missing argument: iqtree2 files directory (see --help)", call.=FALSE)
}

if (file.exists(opt$out)==F) {
  dir.create(opt$out)
}

# set-up conda environment
if (opt$conda) {
  system("conda install -c conda-forge r-essentials r-devtools r-testit r-kmer")
  install.packages("aphid", dependencies = F)
  library("devtools")
  install_github("roblanf/MixtureModelHMM", dependencies = F)
}

# run HMM
library("MixtureModelHMM")

log_info("Reading IQTree2 files...")
alninfo <- system(paste("ls | find ", opt$input, "/*.alninfo", sep=""))
siteprob <- system(paste("ls | find ", opt$input, "/*.siteprob", sep=""))
sitelh <- system(paste("ls | find ", opt$input, "/*.sitelh", sep=""))

log_info("Running HMM...")
hmm_result <- run_HMM(site_info = sitelh, aln_info = alninfo)
save_report(hmm_result = hmm_result, output_filename = paste(opt$out,"/hmm_report.Rmd",sep=""))

# visualization
log_info("Plotting...")
tiff(filename="siteprob.tiff")
plot_scatter(siteprob)
dev.off()

tiff(filename="alignment_plot.tiff")
hmm_result$alignment_plot
dev.off()

log_info("Done!")