#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 7){
stop("Usage: filter_candidates.R <result/high_result.csv> <result/low_result.csv> <inc_taxa> <exc_taxa> <inc_pct> <exc_pct> <filter_mode>")
}
high_file <- args[1]
low_file <-args[2]
inc_taxa <- strsplit(args[3], ",")[[1]]
exc_taxa <- strsplit(args[4], ",")[[1]]
inc_pct <- as.numeric(args[5])
exc_pct <- as.numeric(args[6])
filter_mode <- args[7]
read_results <- function(file){
df <- read_csv(file, show_col_types = FALSE)
required_cols <-c("Gene_ID", "Organism", "Percentage_with_target_genes_complete_genomes")
if (!all(required_cols %in% colnames(df))) {
stop("Required columns not found in the input file: ", paste(setdiff(required_cols, colnames(df)), collapse = ", "))
}
has_draft <- "Percentage_with_target_genes_draft_genomes" %in% colnames(df)
df <- df %>%
select(Gene_ID, Organism, Percentage_with_target_genes_complete_genomes,
all_of(if (has_draft) "Percentage_with_target_genes_draft_genomes" else character(0))) %>%
rename(gene = Gene_ID, organism = Organism, pct_complete = Percentage_with_target_genes_complete_genomes) %>%
mutate(pct_complete = as.numeric(str_remove(pct_complete, "%")),
pct_draft = if (has_draft) as.numeric(str_remove(Percentage_with_target_genes_draft_genomes, "%")) else NA_real_
) %>%
select(-contains("draft_genomes"))
return(df)
}

# select the candidates genes based on different mode (complete for light, optional draft and for heavy
filter_by_mode <- function(df_high, df_low, mode) {
# ensure all inc/exc taxa are present
check_taxa_presence <- function(df, taxa, label) {
missing_taxa <- setdiff(taxa, unique(df$organism))
if (length(missing_taxa) > 0) {
stop(paste("Missing", label, "taxa in results:", paste(missing_taxa, collapse = ", ")))
}
}
check_taxa_presence(df_high, inc_taxa, "inclusivity (high)")
check_taxa_presence(df_low, inc_taxa, "inclusivity (low)")
check_taxa_presence(df_high, exc_taxa, "exclusivity (high)")
check_taxa_presence(df_low, exc_taxa, "exclusivity (low)")

pivot_filter <- function(df, taxa, value_col, threshold, direction = "gte") {
df_wide = df %>%
filter(organism %in% taxa) %>%
select(gene, organism, !!sym(value_col)) %>%
pivot_wider(names_from = organism, values_from = !!sym(value_col), values_fill = 0)
if (direction == "gte") {
df_wide %>% filter(if_all(all_of(taxa), ~ . >= threshold))
} else {
df_wide %>% filter(if_all(all_of(taxa), ~ . <= threshold))
}
}
complete_high_inc <- pivot_filter(df_high, inc_taxa, "pct_complete", inc_pct, "gte")
complete_low_inc <- pivot_filter(df_low, inc_taxa, "pct_complete", inc_pct, "gte")
complete_high_exc <- pivot_filter(df_high, exc_taxa, "pct_complete", exc_pct, "lte")
complete_low_exc <- pivot_filter(df_low, exc_taxa, "pct_complete", exc_pct, "lte")
has_draft <- any(!is.na(df_high$pct_draft)) && any(!is.na(df_low$pct_draft))
if (has_draft) {
check_taxa_presence(df_high, inc_taxa, "inclusivity (draft high)")
check_taxa_presence(df_low, inc_taxa, "inclusivity (draft low)")
check_taxa_presence(df_high, exc_taxa, "exclusivity (draft high)")
check_taxa_presence(df_low, exc_taxa, "exclusivity (draft low)")
draft_high_inc <- pivot_filter(df_high, inc_taxa, "pct_draft", inc_pct, "gte")
draft_low_inc <- pivot_filter(df_low, inc_taxa, "pct_draft", inc_pct, "gte")
draft_high_exc <- pivot_filter(df_high, exc_taxa, "pct_draft", exc_pct, "lte")
draft_low_exc <- pivot_filter(df_low, exc_taxa, "pct_draft", exc_pct, "lte")
}
if (mode == "complete_only") {
keeper_genes <- reduce(list(complete_high_inc, complete_low_inc, complete_high_exc, complete_low_exc),
inner_join, by = "gene")
} else if (mode == "draft_only") {
if (!has_draft) stop("Draft genome data not found. Run GeTPrev in heavy mode to use 'draft_only'.")
keeper_genes <- reduce(list(draft_high_inc, draft_low_inc, draft_high_exc, draft_low_exc),
inner_join, by = "gene")
} else if (mode == "overlap") {
if (!has_draft) stop("Draft genome data not found. Run GeTPrev in heavy mode to use 'overlap'.")
keeper_genes <- reduce(list(complete_high_inc, complete_low_inc, complete_high_exc, complete_low_exc, 
draft_high_inc, draft_low_inc, draft_high_exc, draft_low_exc),
inner_join, by = "gene")
} else {
stop("Invalid filter mode. Please use complete_only, draft_only, or overlap.")
}
return(keeper_genes$gene)
}
# main execution section
df_high <- read_results(high_file)
df_low <- read_results(low_file)
candidates <- filter_by_mode(df_high, df_low, filter_mode)
all_rows <- bind_rows(df_high, df_low) %>% filter(gene %in% candidates)
out_dir <- "result"
if (!dir.exists(out_dir)) dir.create(out_dir)
out_csv <- file.path(out_dir, paste0("final_candidates_", filter_mode, ".csv"))
out_ids <- file.path(out_dir, paste0("final_candidates_", filter_mode, "_ids.txt"))
write_csv(all_rows, out_csv)
write_lines(candidates, out_ids)
cat(length(candidates), "genes written to", out_csv, "and", out_ids, "\n")
