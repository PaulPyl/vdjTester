#!/usr/bin/env Rscript

# compare_results.R
# Usage: compare_results.R <pipeline_csv1> [<pipeline_csv2> ...] <reference_csv> <outdir>
# The script reads multiple pipeline TSV/CSV files and one reference file, merges on a key column, 
# computes V‑gene, J‑gene, and CDR3 accuracy for each pipeline, writes a summary CSV and a bar
# chart plot comparing all pipelines.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: compare_results.R <pipeline_csv1> [<pipeline_csv2> ...] <reference_csv> <outdir>\n")
  quit(status = 1)
}

# Parse arguments: last is outdir, second-to-last is reference, rest are pipeline files
outdir <- args[length(args)]
ref_file <- args[length(args) - 1]
pipeline_files <- args[1:(length(args) - 2)]

# Load data -----------------------------------------------------------
read_data <- function(file_path) {
  ext <- tools::file_ext(file_path)
  if (ext %in% c("tsv", "txt")) {
    read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# Load reference data once
cat("Loading reference data: ", ref_file, "\n")
ref_df <- read_data(ref_file)

# Helper to clean v/j calls ------------------------------------------
clean_call <- function(v) {
  sapply(strsplit(v, ","), function(x) sub("\\*.*", "", x[1]))
}

# Process each pipeline file
all_results <- list()

for (i in seq_along(pipeline_files)) {
  pipeline_file <- pipeline_files[i]
  cat("\n=== Processing pipeline ", i, ": ", pipeline_file, " ===\n")
  
  pipe_df <- read_data(pipeline_file)
  
  # Determine key column
  key <- ifelse("sequence_id" %in% names(pipe_df), "sequence_id", 
                ifelse("query_id" %in% names(pipe_df), "query_id", NULL))
  if (is.null(key) || !(key %in% names(pipe_df))) {
    stop("Could not find a suitable key column in pipeline output: ", pipeline_file)
  }
  
  # Merge
  merged <- merge(pipe_df, ref_df, by = key, all = FALSE, suffixes = c("_pipe", "_ref"))
  cat(sprintf("Shared sequences: %d\n", nrow(merged)))
  
  # V‑gene accuracy
  if ("v_call_pipe" %in% names(merged) && "v_call_ref" %in% names(merged)) {
    merged$v_call_pipe_clean <- clean_call(merged$v_call_pipe)
    merged$v_call_ref_clean <- clean_call(merged$v_call_ref)
    merged$v_match <- merged$v_call_pipe_clean == merged$v_call_ref_clean
    v_accuracy <- mean(merged$v_match, na.rm = TRUE)
  } else {
    v_accuracy <- NA
  }
  
  # J‑gene accuracy
  if ("j_call_pipe" %in% names(merged) && "j_call_ref" %in% names(merged)) {
    merged$j_call_pipe_clean <- clean_call(merged$j_call_pipe)
    merged$j_call_ref_clean <- clean_call(merged$j_call_ref)
    merged$j_match <- merged$j_call_pipe_clean == merged$j_call_ref_clean
    j_accuracy <- mean(merged$j_match, na.rm = TRUE)
  } else {
    j_accuracy <- NA
  }
  
  # CDR3 accuracy
  pipe_cdr3_col <- if ("cdr3_aa_pipe" %in% names(merged)) "cdr3_aa_pipe" else if ("junction_aa_pipe" %in% names(merged)) "junction_aa_pipe" else NULL
  ref_cdr3_col <- if ("cdr3_aa_ref" %in% names(merged)) "cdr3_aa_ref" else if ("junction_aa_ref" %in% names(merged)) "junction_aa_ref" else if ("cdr3_aa" %in% names(merged)) "cdr3_aa" else NULL
  
  if (!is.null(pipe_cdr3_col) && !is.null(ref_cdr3_col)) {
    merged$cdr3_match <- merged[[pipe_cdr3_col]] == merged[[ref_cdr3_col]]
    cdr3_accuracy <- mean(merged$cdr3_match, na.rm = TRUE)
  } else {
    cdr3_accuracy <- NA
  }
  
  # Print accuracies
  cat(sprintf("V‑Gene Accuracy: %0.2f%%\n", v_accuracy * 100))
  cat(sprintf("J‑Gene Accuracy: %0.2f%%\n", j_accuracy * 100))
  cat(sprintf("CDR3 Accuracy: %0.2f%%\n", cdr3_accuracy * 100))
  
  # Store results
  pipeline_name <- sub("\\.[^.]*$", "", basename(pipeline_file))
  all_results[[i]] <- data.frame(
    Pipeline = pipeline_name,
    Shared_Sequences = nrow(merged),
    V_Gene_Accuracy = v_accuracy,
    J_Gene_Accuracy = j_accuracy,
    CDR3_Accuracy = cdr3_accuracy,
    stringsAsFactors = FALSE
  )
}

# Write summary
cat("\n=== Summary ===\n")
summary_df <- do.call(rbind, all_results)
print(summary_df)
write.csv(summary_df, file = file.path(outdir, "comparison_summary.csv"), row.names = FALSE)

# Plot
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# Reshape for plotting
plot_df <- reshape2::melt(summary_df, id.vars = "Pipeline", 
                          measure.vars = c("V_Gene_Accuracy", "J_Gene_Accuracy", "CDR3_Accuracy"),
                          variable.name = "Metric", value.name = "Accuracy")
plot_df$Metric <- sub("_Accuracy", "", plot_df$Metric)

if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2", repos = "https://cloud.r-project.org")
}

p <- ggplot(plot_df, aes(x = Metric, y = Accuracy, fill = Pipeline)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 1.05) +
  ggtitle("Pipeline Comparison with Reference") +
  ylab("Accuracy (Exact Match)") +
  xlab("Metric") +
  geom_text(aes(label = sprintf("%0.1f%%", Accuracy * 100)), vjust = -0.5, position = position_dodge(width = 0.9), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(outdir, "comparison_plot.png"), plot = p, width = 10, height = 6)

cat("Results written to ", outdir, "\n")
library(ggplot2)

metrics <- c("V_Gene", "J_Gene", "CDR3")
values <- c(v_accuracy, j_accuracy, cdr3_accuracy)
plot_df <- data.frame(Metric = metrics, Accuracy = values)

p <- ggplot(plot_df, aes(x = Metric, y = Accuracy)) +
  geom_bar(stat = "identity", fill = c("skyblue", "lightgreen", "salmon")) +
  ylim(0, 1.05) +
  ggtitle(sprintf("Comparison with Reference (n=%d)", nrow(merged))) +
  ylab("Accuracy (Exact Match)") +
  geom_text(aes(label = sprintf("%0.1f%%", Accuracy * 100)), vjust = -0.5)

ggsave(filename = file.path(outdir, "comparison_plot.png"), plot = p, width = 8, height = 6)

cat("Results written to ", outdir, "\n")

