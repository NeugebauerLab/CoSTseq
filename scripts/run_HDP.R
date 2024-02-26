### PURPOSE OF THIS SCRIPT
## Run HDProbe using Rscript with command line options.

# Load dependencies ------------------------------------------------------------

library(tidyverse)
library(data.table)
library(devtools)
library(optparse)

# Parse command line arguments to retrieve input files

# Define the option parser
option_list <- list(
  make_option(c("-a", "--control"), type = "character", 
              help = "Name of control sample (must be same as in column header)."),
  make_option(c("-b", '--sample'), type = "character", 
              help = "Name of test sample (must be same as in column header)."),
  make_option(c("--hdploc"), type = "character", 
              help = "Location of the HDProbe package."),
  make_option(c("--mutcov"), type = "character", 
              help = "Input count table filename."),		    
  make_option(c("--outpath"), type = "character", 
              help = "Path to the folder where the output should be saved."),		    
  make_option(c("--cutoff"), type = "numeric",
              help = "Minimum coverage required for a nt to be included in the analysis.")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# load HDProbe
setwd(opt$hdploc)
devtools::load_all()

# Functions to make life easier
inv_logit <- function(x) exp(x)/(exp(x) + 1)
logit <- function(x) log(x/(1-x))


# HDProbe wrapper that generates useful reports
#' @param data Table in the form you have been recently sending me.
#' @param metadf Table with three columns; sample_name, E_ID, and R_ID. Contents
#' described in example code at bottom of this script
#' @param fit_name File name to save HDProbe Fit object to
#' @param results_path Path to where you want non-plot results saved
#' @param nreps Number of replicates in each experimental condition
#' @param Normalize Logical; if TRUE, mutation rates will be normalized
#' to the median mutation rate in all samples.
#' @param Homosked Logical; if TRUE, a homoskedastic variance model will be used,
#' rather than the default reads vs. variance heteroskedastic model. Nice if you
#' are analyzing a single RNA, or generally low-throughput data.
#' @param prior_slope Numeric; heteroskedastic slope to impute if not enough sites
#' make it past filtering. This was chosen to be equal to what was seen in
#' other datasets you have sent me, in cases where the data is too low-throughput
#' to confidently estimate the variance trend de novo
#' @param prior_int_het Numeric; heteroskedastic intercept to impute if not enough sites
#' make it past filtering. This was chosen to be equal to what was seen in
#' other datasets you have sent me, in cases where the data is too low-throughput
#' to confidently estimate the variance trend de novo
analyze_dms <- function(data, metadf,
                        results_name, results_path,
                        nreps = 3, Normalize = TRUE, Homosked = FALSE,
                        prior_slope = -0.6, prior_int_het = 0.7, cutoff = 700){

  message("Processing data")
  Muts_df <- data_prep(data, metadf)
  # Takes a minute or two

  Muts_df <- Muts_df %>%
    mutate(GF = paste(rname, gloc, nt, strand, sep = "&"))

  ### Analyze only As and Cs
  Muts_df <- Muts_df[Muts_df$nt %in% c("A", "C", "U", "T"),]


  message("Running HDProbe")
  Fit <- HDProbe(Muts_df, nreps = nreps,
                 homosked = Homosked, Normalize = Normalize,
                 prior_slope = prior_slope, prior_int_het = prior_int_het,
                 cutoff = cutoff)
  # This takes a bit

  message("Saving results")

  ### Create report
  GF_dict <- Muts_df %>%
    dplyr::select(rname, gloc, nt, strand, P_ID) %>%
    dplyr::distinct()

  Diffmut <- Fit$Diffmut

  # Add genomic info
  Diffmut <- left_join(Diffmut, GF_dict, by = "P_ID")


  # Make table of results
  new_df <- tibble(P_ID = Diffmut$P_ID,
                   chromosome = Diffmut$rname,
                   position = Diffmut$gloc,
                   strand = Diffmut$strand,
                   nt = Diffmut$nt,
                   tot_reads = Diffmut$reads,
                   avg_reads = Diffmut$reads/(nreps*2),
                   pval = Diffmut$lpval,
                   padj_BH = Diffmut$lpadj,
                   padj_BF = p.adjust(Diffmut$lpval, method = "bonferroni"),
                   difference = Diffmut$l_num,
                   uncertainty = Diffmut$lden)


  ## Add average mutation rate in wt and mutant
  Mutrates <- Fit$Mut_est

  mutrate_dict <- tibble(P_ID = Mutrates$P_ID[Mutrates$E_ID == 1],
                         wt_raw_mutrate = Mutrates$nmuts[Mutrates$E_ID == 1]/Mutrates$ntrials[Mutrates$E_ID == 1],
                         exp_raw_mutrate = Mutrates$nmuts[Mutrates$E_ID == 2]/Mutrates$ntrials[Mutrates$E_ID == 2],
                         wt_adj_mutrate = inv_logit(Mutrates$Avg_lp[Mutrates$E_ID == 1]),
                         exp_adj_mutrate = inv_logit(Mutrates$Avg_lp[Mutrates$E_ID == 2]))


  new_df <- inner_join(new_df, mutrate_dict, by = "P_ID")


  # Order by sig
  new_df <- new_df[order(new_df$pval),]

  # Save
  write_csv(new_df, file = results_path)

}

# Function that processes data for HDProbe
data_prep <- function(cU, metadf){

  # How many replicates? (have to be balanced)
  nreps <- max(metadf$R_ID)

  # add sample ID info to data
  Muts_df <- right_join(cU, metadf, by = "sample_name")

  # extract useful columns !!(this is hard coded and may have to be changed if column names are different)!!
  Muts_df <- Muts_df[,c("chromosome", "genomic_position", "nt", "strand", "n_mutations", "n_reads", "E_ID", "R_ID")]

  # Change column names to make HDProbe happy
  colnames(Muts_df) <- c("rname", "gloc", "nt", "strand", "nmuts", "ntrials", "E_ID", "R_ID")

  # Check that data exists for each GF in each sample
  Test_df <- Muts_df %>% mutate(check = 1) %>% group_by(rname, gloc, nt, strand) %>%
    summarise(check = sum(check)) %>% filter(check == nrow(metadf)) %>% select(rname, gloc, nt, strand)

  # Filter out sites lacking data in one or more replicates
  Muts_df <- inner_join(Muts_df, Test_df, by = c("rname", "gloc", "nt", "strand"))

  # Order to prep for position ID addition
  Muts_df <- Muts_df[order(Muts_df$nt, Muts_df$gloc, Muts_df$rname, Muts_df$strand, Muts_df$E_ID, Muts_df$R_ID),]

  # Each position will appear n times in a row, where n = 2*nreps
  Muts_df$P_ID <- rep(1:(nrow(Muts_df)/(2*nreps)), each = nreps*2)

  return(Muts_df)

}

# Analyze----------------------------------------------------------------

# load
data <- fread(opt$mutcov)

# filter for sample of interest
data_ctrl <- data[grepl(opt$control,data$sample_name)]
data_test <- data[grepl(opt$sample,data$sample_name)]

# Combine datasets
data <- bind_rows(data_ctrl, data_test)

# Metadata
metadf <- tibble(sample_name = unique(data$sample_name),
                 E_ID = c(1, 1, 1, 2, 2, 2),
                 R_ID = c(1, 2, 3, 1, 2, 3))

# Call the analyze_dms function with updated parameters
analyze_dms(data, metadf, cutoff = opt$cutoff,
            results_path = opt$outpath)
