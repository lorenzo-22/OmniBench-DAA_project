#!/usr/bin/env Rscript

# Limma differential expression analysis for OmniBenchmark
library(limma)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--name", type="character", help="Dataset name"),
  make_option("--data.matrix", type="character", dest="data_matrix", help="Input data file")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

cat("Running limma analysis\n")
cat("Input:", args$data_matrix, "\n")
cat("Output dir:", args$output_dir, "\n")

# Create output directory
dir.create(args$output_dir, showWarnings=FALSE, recursive=TRUE)

# Read data
data <- read.csv(args$data_matrix, row.names=1)

# Remove the last column if it's the label column
if ("is_differentially_expressed" %in% colnames(data)) {
  data <- data[, !colnames(data) %in% "is_differentially_expressed"]
}

cat("Data dimensions:", dim(data), "\n")
cat("Columns:", colnames(data), "\n")

# Define groups (assuming A* samples are group 1, B* samples are group 2)
group <- factor(ifelse(grepl("^A", colnames(data)), "A", "B"))
cat("Groups:", as.character(group), "\n")

# Create design matrix
design <- model.matrix(~group)
colnames(design) <- c("Intercept", "B_vs_A")

# Fit linear model
fit <- lmFit(data, design)
fit <- eBayes(fit)

# Extract results
results <- topTable(fit, coef="B_vs_A", number=Inf, sort.by="none")

# Add feature names
results$Feature <- rownames(data)

# Reorder columns
results <- results[, c("Feature", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# Determine significance
results$Significant <- results$adj.P.Val < 0.05

# Save results
output_file <- file.path(args$output_dir, paste0(args$name, "_limma_results.csv"))
write.csv(results, output_file, row.names=FALSE)

cat("Results saved to:", output_file, "\n")
cat("Total features:", nrow(results), "\n")
cat("Significant (FDR < 0.05):", sum(results$Significant), "\n")
