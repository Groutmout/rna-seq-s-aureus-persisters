#!/usr/bin/env Rscript

# ----------------------------------------------
# DESeq2 differential expression analysis
# Compatible with R 3.4.1 / DESeq2 1.16
# ----------------------------------------------

suppressPackageStartupMessages({
    library("DESeq2")
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: run_deseq2.R <counts_files...> <output.csv>")
}

# Last argument = output file
output_file <- args[length(args)]
count_files <- args[-length(args)]

cat("Count files:\n")
print(count_files)
cat("\nOutput:", output_file, "\n\n")

# ----------------------------------------------
# 1) Read and merge featureCounts tables
# ----------------------------------------------

read_fc <- function(file) {
    df <- read.table(file, header = TRUE, sep = "\t", comment.char = "#")
    # Keep only gene ID and counts
    df <- df[, c("Geneid", grep("^\\S+\\.bam$", names(df), value = TRUE))]
    names(df)[2] <- file  # label column by filename
    return(df)
}

list_df <- lapply(count_files, read_fc)

# Merge on Geneid
counts <- Reduce(function(x, y) merge(x, y, by="Geneid"), list_df)

# Remove Geneid column (store separately)
gene_ids <- counts$Geneid
count_matrix <- counts[, -1]

# Clean column names
colnames(count_matrix) <- gsub(".*/|\\.bam$", "", colnames(count_matrix))

# Convert to integer matrix
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "integer"

# ----------------------------------------------
# 2) Create metadata (conditions)
# ----------------------------------------------

# For the hackathon: 6 samples
# SRR10452052-54 → persister
# SRR10452055-57 → control
# (Adjust if needed)
conditions <- ifelse(grepl("SRR1045205[2-4]", colnames(count_matrix)),
                     "persister", "control")

coldata <- data.frame(
    row.names = colnames(count_matrix),
    condition = factor(conditions)
)

cat("Sample table:\n")
print(coldata)
cat("\n")

# ----------------------------------------------
# 3) DESeq2 analysis
# ----------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

# Add gene IDs back
res_df <- as.data.frame(res)
res_df$gene_id <- gene_ids

# ----------------------------------------------
# 4) Output
# ----------------------------------------------

write.csv(res_df, file = output_file, row.names = FALSE)
cat("DESeq2 results written to:", output_file, "\n")
