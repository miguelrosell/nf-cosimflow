#!/usr/bin/env Rscript
# compute_cosine.R
#
# Usage:
# Rscript compute_cosine.R --input file.csv --out_prefix cosine --method cosine --min_gene_mean 0.1 --sample_cols "Con-1,Con-2,Sh-1"
#
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(readxl)
  library(dplyr)
  library(lsa)        # for cosine
  library(pheatmap)
})

option_list = list(
  make_option(c("-i", "--input"), type="character", help="Input expression file (csv or xlsx)", metavar="file"),
  make_option(c("-o", "--out_prefix"), type="character", default="cosine", help="Output prefix"),
  make_option(c("-m", "--method"), type="character", default="cosine", help="method: cosine|pearson|spearman"),
  make_option(c("-g", "--min_gene_mean"), type="double", default=0.0, help="Minimum mean expression across samples to keep gene"),
  make_option(c("-s", "--sample_cols"), type="character", default="auto", help="Comma-separated sample columns OR 'auto' to detect numeric columns")
)

opt = parse_args(OptionParser(option_list=option_list))

# ---- Load table (support csv/tsv/xlsx) ----
infile <- opt$input
if (grepl("\\.xlsx?$", infile, ignore.case=TRUE)) {
  df <- read_excel(infile)
} else {
  df <- read_csv(infile, col_types = cols(.default = col_guess()))
}

# If first column is gene id named 'id' or similar, set rownames
if ("id" %in% tolower(names(df)) || "gene" %in% tolower(names(df))) {
  # pick the first col as rownames
  first_col <- names(df)[1]
  rownames(df) <- df[[first_col]]
  df[[first_col]] <- NULL
}

# ---- Select sample columns ----
if (opt$sample_cols != "auto") {
  cols <- unlist(strsplit(opt$sample_cols, ","))
  cols <- trimws(cols)
  if (!all(cols %in% colnames(df))) stop("Some sample columns not found in data")
  expr <- df[, cols, drop=FALSE]
} else {
  # automatic: take numeric columns
  is_num <- sapply(df, is.numeric)
  if (sum(is_num) < 2) stop("No numeric sample columns detected. Specify --sample_cols explicitly.")
  expr <- df[, is_num, drop=FALSE]
}

# ---- Convert to numeric matrix ----
mat <- as.matrix(expr)
mode(mat) <- "numeric"

# ---- Filter low-expression genes ----
if (!is.null(opt$min_gene_mean) && opt$min_gene_mean > 0) {
  keep <- rowMeans(mat, na.rm=TRUE) >= opt$min_gene_mean
  mat <- mat[keep, , drop=FALSE]
}

# ---- Compute similarity ----
method <- tolower(opt$method)
if (method == "cosine") {
  # lsa::cosine expects columns = vectors, so transpose
  sim <- lsa::cosine(t(mat))
} else if (method %in% c("pearson","spearman")) {
  sim <- cor(mat, method = method, use = "pairwise.complete.obs")
} else {
  stop("Unsupported method: choose cosine, pearson or spearman")
}

# label rows/cols
colnames(sim) <- colnames(mat)
rownames(sim) <- colnames(mat)

# ---- Write outputs ----
out_prefix <- opt$out_prefix
write.csv(sim, paste0(out_prefix, "_matrix.tsv"), quote = FALSE)

# Save heatmap (png)
png(paste0(out_prefix, "_heatmap.png"), width = 900, height = 700)
pheatmap(sim,
         display_numbers = TRUE,
         main = paste0(toupper(substr(method,1,1)), substr(method,2,nchar(method)), " similarity"),
         color = colorRampPalette(c("white","navy"))(50),
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

cat("Outputs written:", paste0(out_prefix, "_matrix.tsv"), "and", paste0(out_prefix, "_heatmap.png"), "\n")
