#!/usr/bin/env Rscript

# compare_results_detailed.R
# Usage: compare_results_detailed.R <pipeline_csv1> [<pipeline_csv2> ...] <reference_csv> <outdir>
# The script reads multiple pipeline TSV/CSV files and one reference file, computes
# V-gene, J-gene, and CDR3 accuracy for each pipeline, generates visualizations
# including Venn diagrams, and saves data for Quarto report rendering.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: compare_results_detailed.R <pipeline_csv1> [<pipeline_csv2> ...] <reference_csv> <outdir>\n")
  quit(status = 1)
}

# Parse arguments
outdir <- args[length(args)]
ref_file <- args[length(args) - 1]
pipeline_files <- args[1:(length(args) - 2)]

# Create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load data -----------------------------------------------------------
read_data <- function(file_path) {
  ext <- tools::file_ext(file_path)
  if (ext %in% c("tsv", "txt")) {
    read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# Install/load required packages
install_if_needed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing ", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}

install_if_needed("ggplot2")
install_if_needed("reshape2")
install_if_needed("VennDiagram")

# Load reference data once
cat("Loading reference data: ", ref_file, "\n")
ref_df <- read_data(ref_file)

# Helper to clean v/j calls
clean_call <- function(v) {
  sapply(strsplit(v, ","), function(x) sub("\\*.*", "", x[1]))
}

# Process each pipeline file
all_results <- list()
venn_data <- list()

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
  
  pipeline_name <- sub("\\.[^.]*$", "", basename(pipeline_file))
  
  # Get unique sequence IDs
  pipeline_ids <- unique(pipe_df[[key]])
  reference_ids <- unique(ref_df[[key]])
  
  # Merge on common sequences
  merged <- merge(pipe_df, ref_df, by = key, all = FALSE, suffixes = c("_pipe", "_ref"))
  shared_ids <- merged[[key]]
  
  shared_count <- length(shared_ids)
  pipeline_only_count <- length(setdiff(pipeline_ids, reference_ids))
  reference_only_count <- length(setdiff(reference_ids, pipeline_ids))
  
  cat(sprintf("Total pipeline sequences: %d\n", length(pipeline_ids)))
  cat(sprintf("Total reference sequences: %d\n", length(reference_ids)))
  cat(sprintf("Shared sequences: %d\n", shared_count))
  cat(sprintf("Pipeline-only sequences: %d\n", pipeline_only_count))
  cat(sprintf("Reference-only sequences: %d\n", reference_only_count))
  
  # V-gene accuracy
  if ("v_call_pipe" %in% names(merged) && "v_call_ref" %in% names(merged)) {
    merged$v_call_pipe_clean <- clean_call(merged$v_call_pipe)
    merged$v_call_ref_clean <- clean_call(merged$v_call_ref)
    merged$v_match <- merged$v_call_pipe_clean == merged$v_call_ref_clean
    v_accuracy <- mean(merged$v_match, na.rm = TRUE)
  } else {
    v_accuracy <- NA
  }
  
  # J-gene accuracy
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
  cat(sprintf("V-Gene Accuracy: %0.2f%%\n", v_accuracy * 100))
  cat(sprintf("J-Gene Accuracy: %0.2f%%\n", j_accuracy * 100))
  cat(sprintf("CDR3 Accuracy: %0.2f%%\n", cdr3_accuracy * 100))
  
  # Store results
  all_results[[i]] <- data.frame(
    Pipeline = pipeline_name,
    Total_Sequences = length(pipeline_ids),
    Shared_Sequences = shared_count,
    Pipeline_Only = pipeline_only_count,
    Reference_Only = reference_only_count,
    V_Gene_Accuracy = v_accuracy,
    J_Gene_Accuracy = j_accuracy,
    CDR3_Accuracy = cdr3_accuracy,
    stringsAsFactors = FALSE
  )
  
  # Store Venn diagram data
  venn_data[[i]] <- list(
    pipeline_name = pipeline_name,
    pipeline_ids = pipeline_ids,
    reference_ids = reference_ids
  )
}

# Write summary
cat("\n=== Summary ===\n")
summary_df <- do.call(rbind, all_results)
print(summary_df)
write.csv(summary_df, file = file.path(outdir, "comparison_summary.csv"), row.names = FALSE)

# Generate Venn diagrams for each pipeline
for (i in seq_along(venn_data)) {
  data <- venn_data[[i]]
  pipeline_name <- data$pipeline_name
  
  # Create Venn diagram
  venn_file <- file.path(outdir, paste0("venn_", pipeline_name, ".png"))
  
  png(venn_file, width = 600, height = 600)
  draw.pairwise.venn(
    area1 = length(data$pipeline_ids),
    area2 = length(data$reference_ids),
    cross.area = length(intersect(data$pipeline_ids, data$reference_ids)),
    category = c("Pipeline", "Reference"),
    main = paste("Sequence Overlap -", pipeline_name),
    cex = 1.5,
    cat.cex = 1.2
  )
  dev.off()
  
  cat("Generated Venn diagram: ", venn_file, "\n")
}

# Generate accuracy comparison plot
plot_df <- melt(summary_df, id.vars = "Pipeline", 
                measure.vars = c("V_Gene_Accuracy", "J_Gene_Accuracy", "CDR3_Accuracy"),
                variable.name = "Metric", value.name = "Accuracy")
plot_df$Metric <- sub("_Accuracy", "", plot_df$Metric)

p <- ggplot(plot_df, aes(x = Metric, y = Accuracy, fill = Pipeline)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 1.05) +
  ggtitle("Pipeline Accuracy Comparison with Reference") +
  ylab("Accuracy (Exact Match)") +
  xlab("Metric") +
  geom_text(aes(label = sprintf("%0.1f%%", Accuracy * 100)), vjust = -0.5, 
            position = position_dodge(width = 0.9), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "gray90"))

ggsave(filename = file.path(outdir, "accuracy_comparison.png"), plot = p, width = 10, height = 6, dpi = 300)

# Generate sequence overlap bar plot
overlap_df <- summary_df[, c("Pipeline", "Shared_Sequences", "Pipeline_Only", "Reference_Only")]
overlap_melted <- melt(overlap_df, id.vars = "Pipeline", 
                       variable.name = "Category", value.name = "Count")
overlap_melted$Category <- factor(overlap_melted$Category, 
                                   levels = c("Shared_Sequences", "Pipeline_Only", "Reference_Only"))

p_overlap <- ggplot(overlap_melted, aes(x = Pipeline, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  ggtitle("Sequence Coverage by Pipeline") +
  ylab("Sequence Count") +
  xlab("Pipeline") +
  scale_fill_manual(values = c("Shared_Sequences" = "green", 
                               "Pipeline_Only" = "orange", 
                               "Reference_Only" = "red"),
                    labels = c("Shared with Reference", "Pipeline Only", "Reference Only")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "gray90"))

ggsave(filename = file.path(outdir, "sequence_overlap.png"), plot = p_overlap, width = 10, height = 6, dpi = 300)

# Save summary data as RDS for Quarto report
saveRDS(summary_df, file = file.path(outdir, "summary_data.rds"))

cat("\nResults written to ", outdir, "\n")
