#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(readxl)
  library(dplyr)
  library(lsa)
  library(pheatmap)
})

option_list = list(
  make_option(c("-i","--input"), type="character", help="Input expression file (csv or xlsx)"),
  make_option(c("-o","--out_prefix"), type="character", default="cosine"),
  make_option(c("-m","--method"), type="character", default="cosine"),
  make_option(c("-g","--min_gene_mean"), type="double", default=0.0),
  make_option(c("-s","--sample_cols"), type="character", default="auto")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---- Load table ----
if (grepl("\\.xlsx?$", opt$input, ignore.case=TRUE)) {
  df <- as.data.frame(read_excel(opt$input))
} else {
  df <- as.data.frame(read_csv(opt$input, col_types = cols(.default = col_guess())))
}

# ---- Set rownames if first column is gene/id ----
first_col <- names(df)[1]
if (tolower(first_col) %in% c("gene","id")) {
  if (is.character(df[[first_col]]) || is.factor(df[[first_col]])) {
    rownames(df) <- df[[first_col]]
    df[[first_col]] <- NULL
  }
}

# ---- Select sample columns ----
if (opt$sample_cols != "auto") {
  cols <- trimws(unlist(strsplit(opt$sample_cols,",")))
  if (!all(cols %in% colnames(df))) stop("Some sample columns not found in data")
  mat <- as.matrix(df[, cols, drop=FALSE])
} else {
  num_cols <- sapply(df, is.numeric)
  if (sum(num_cols) < 2) stop("Need at least 2 numeric columns")
  mat <- as.matrix(df[, num_cols, drop=FALSE])
}

# ---- Filter low-expression genes ----
if (opt$min_gene_mean > 0) {
  keep <- rowMeans(mat, na.rm=TRUE) >= opt$min_gene_mean
  mat <- mat[keep, , drop=FALSE]
}

# ---- Compute similarity ----
method <- tolower(opt$method)
if (method == "cosine") {
  # lsa::cosine computes similarity between columns
  sim <- lsa::cosine(mat)
} else if (method %in% c("pearson","spearman")) {
  sim <- cor(mat, method = method, use = "pairwise.complete.obs")
} else {
  stop("Unsupported method: cosine, pearson, spearman")
}

# Safety check: sim must be square and match number of samples
if (!is.matrix(sim) || nrow(sim) != ncol(sim)) {
  stop("Unexpected similarity matrix dimensions: check input matrix 'mat'")
}
if (ncol(sim) != ncol(mat)) {
  stop("Number of columns in similarity matrix does not match number of samples")
}

# ---- Label rows/cols ----
colnames(sim) <- colnames(mat)
rownames(sim) <- colnames(mat)


# ---- Save outputs ----
out_matrix <- paste0(opt$out_prefix, "_matrix.csv")
write.csv(sim, out_matrix, quote=FALSE, row.names=TRUE)

out_heat <- paste0(opt$out_prefix, "_heatmap.png")
png(out_heat, width=900, height=700)
pheatmap(sim, display_numbers=TRUE,
         main=paste0(toupper(substr(method,1,1)), substr(method,2,nchar(method)), " similarity"),
         color=colorRampPalette(c("white","navy"))(50),
         cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

cat("Outputs written:", out_matrix, "and", out_heat, "\n")
