# this code is used to set up environment and run HMM using Conda
# required packages: logger, optparse, data.table, MixtureModelHMM (https://github.com/roblanf/MixtureModelHMM)

library(logger)
library(optparse)
library(data.table)

# retrieve the input parameters
c_dir <- getwd()

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="directory of iqtree2 files", metavar="character"),
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

# run HMM
library("MixtureModelHMM")

log_info("Reading IQTree2 files...")
alninfo <- system(paste("ls | find ", opt$input, "/*.alninfo", sep=""), intern = T)
siteprob <- system(paste("ls | find ", opt$input, "/*.siteprob", sep=""), intern = T)
sitelh <- system(paste("ls | find ", opt$input, "/*.sitelh", sep=""), intern = T)

log_info("Running HMM...")
hmm_result <- run_HMM(site_info = sitelh, aln_info = alninfo)
save_report(hmm_result = hmm_result, output_filename = paste(opt$out,"/hmm_report.Rmd",sep=""))

# save files
write.table(hmm_result$hmm_transition_table, file = paste(opt$out,"/hmm_transition_table.txt",sep=""), quote=F, row.names=F)
fwrite(hmm_result$classification, file = paste(opt$out,"/hmm_classification.txt",sep=""))

# visualization
log_info("Plotting...")
tiff(filename=paste(opt$out,"/hmm_siteprob.tiff",sep=""), units="px", width=1000, height=1000)
plot_scatter(siteprob)
dev.off()

tiff(filename=paste(opt$out,"/hmm_alignment_plot.tiff",sep=""), units="px", width=1000, height=1000)
hmm_result$alignment_plot
dev.off()

log_info("Done!")