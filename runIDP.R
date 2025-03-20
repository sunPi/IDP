#---- Handle Requirements ----
# Install necessary packages if not already installed and load libraries
# BiocManager::install("sva")
# BiocManager::install("genefilter")
# remotes::install_github("icbi-lab/immunedeconv")
pkgs <- c("docopt", "readr", "dplyr", "remotes", "immunedeconv", "mesocore", "here")
suppressMessages(mesocore::handleRequirements(pkgs))

#---- Load Functions ----
# Function to estimate the deconvolution score
estimate.score <- function(gene_expression_matrix){
  tryCatch({
    message("Calculating deconvolution estimate...")
    return(immunedeconv::deconvolute_estimate(gene_expression_matrix))
  }, error = function(e) {
    cat("Error in deconvolution: ", e$message, "\n")
    NULL
  })
}
# Function to run deconvolution and save results
run_deconvolution <- function(method, gene_expression_matrix, indications = NULL) {
  result <- tryCatch({
    if (is.null(indications)) {
      immunedeconv::deconvolute(gene_expression_matrix, method)
    } else {
      immunedeconv::deconvolute(gene_expression_matrix, method, indication = indications)
    }
  }, error = function(e) {
    cat(sprintf("Error in %s deconvolution: %s\n", method, e$message))
    NULL
  })
  
  # if (!is.null(result)) {
  #   output_file <- paste0(method, "_results.csv")
  #   write.csv(result, output_file, row.names = FALSE)
  #   cat(sprintf("%s deconvolution completed successfully and results saved to %s.\n", method, output_file))
  # } else {
  #   cat(sprintf("No results for %s deconvolution due to an error.\n", method))
  # }
  
  return(result)
}

#---- Header ----
# Please make sure there are no gene duplicates in your dataset !!!!!!!
"Immune Deconvolution Pipeline - Deconvolutes TPM-normalized RNA-seq data by 
  methods such as 
  - 'quantiseq'
  - 'epic'
  - 'mcp_counter' 
  - 'abis' 
  - 'xcell'
  - 'consensus_tme' 
  - 'timer'

Usage: runIDP.R [options]

Options:
  -h --help                     Show this screen.
  -d --data_path=<PATH>         Path to the cleaned and deduped counts data set (RNA-seq data).
  -o --outfolder=<PATH>         Folder where the results are saved.
  -v --verbose=<BOOLEAN>        If set to 1, prints messages verbously.
  -V --version
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

cts.path  <- normalizePath(arguments$data_path)
outfolder <- normalizePath(arguments$outfolder, mustWork = F)
verbose   <- as.integer(arguments$verbose)

if(as.integer(arguments$verbose) == 1){
  verbose <- T
} else{
  verbose <- F
}

#---- Arguments ----
# cts.path  <- normalizePath("data/cleaned_gene_expression_data.csv")
# outfolder <- "results/TMA_DECONVOLUTION"
# verbose   <- as.integer(arguments$verbose)

# Load your cleaned gene expression data
gene_expression_matrix <- readFile(cts.path)
dir.create(outfolder, F, T)

# Ensure gene names are in the first column and set as row names
if (!is.null(gene_expression_matrix[, 1])) {
  rownames(gene_expression_matrix) <- gene_expression_matrix[, 1]
  gene_expression_matrix <- gene_expression_matrix[, -1] # Remove gene name column from data
} else {
  stop("Gene names are not present in the data.")
}

# Check if all remaining columns are numeric
if (!all(sapply(gene_expression_matrix, is.numeric))) {
  stop("Expression data should contain only numeric values.")
}

# Print first few rows for verification
if(verbose){
  print(head(gene_expression_matrix))
}

# Run ESTIMATE of deconvolution
estimate.score <- estimate.score(gene_expression_matrix)

# Check if estimate.score is NULL (in case of an error)
if (is.null(estimate.score)) {
  message("ESTIMATE deconvolution did not complete successfully.\n")
} else {
  # Save the results to a CSV file
  write.csv(estimate.score, here(outfolder, "deconvolution_estimate_score.csv"), row.names = TRUE)
  message("ESTIMATE deconvolution completed successfully and results saved.\n")
}

# List of deconvolution methods to run
methods <- c("quantiseq", "epic", "mcp_counter", "abis", "xcell", "consensus_tme", "timer")

# Create an indication vector with appropriate lengths (for consensus_tme and timer)
indications <- rep("meso", ncol(gene_expression_matrix))

# Run deconvolution methods
decon.res <- list()
for (method in methods) {
  if (method %in% c("consensus_tme", "timer")) {
    decon.res[[method]] <- run_deconvolution(method, gene_expression_matrix, indications)
    decon.res[[method]]$Method <- method
  } else {
    decon.res[[method]] <- run_deconvolution(method, gene_expression_matrix)
    decon.res[[method]]$Method <- method
  }
}

# Merge datasets"
decon.res.merged <- Reduce(function(x, y) merge(x, y, all = TRUE), decon.res)

# Serialize results
message("Serializing the results into ", here(outfolder), "...")
writexl::write_xlsx(decon.res, here(outfolder, "decon_methods.xlsx"))
write.csv(decon.res.merged, here(outfolder,"decon_all.csv"), row.names = FALSE)
message("IDP Finished Successfully!")


