library(tidyverse)
library(corrplot)  
library(Hmisc)

#Figure path 
figure_path <- ""

# Path to the enrichment analysis for all dataset compositions created by following the script in Figure 2
gene_level_enrichment_dataset1 <- ""
gene_level_enrichment_dataset2 <- ""
gene_level_enrichment_dataset3 <- ""
gene_level_enrichment_dataset4 <- ""
gene_level_enrichment_dataset5 <- ""

# Preload and merge dataset descriptors
data_complete <- list(
  gene_level_enrichment_dataset1,
  gene_level_enrichment_dataset2,
  gene_level_enrichment_dataset3,
  gene_level_enrichment_dataset4,
  gene_level_enrichment_dataset5
) %>% 
  map_dfr(~ read_delim(.x, delim = "\t") %>% distinct(Gene, label, label2, family)) %>% 
  distinct()

# Function to compute enrichment metrics for each dataset
process_df <- function(df, newname){
  df %>% 
    filter((pathogenic_match + pathogenic_out) >9,
           (control_out + control_match) >9,
           (pathogenic_match + control_match) >0) %>%
    mutate(sens = pathogenic_match/(pathogenic_match + pathogenic_out),
           specf = control_out /(control_out + control_match),
           PLR = sens/(1-specf),
           NLR = (1-sens)/specf) %>% 
    mutate(log_plr = log(PLR),
           se_log_plr = sqrt(1/pathogenic_match - 1/(pathogenic_match+pathogenic_match) + 1/control_match - 1/(control_match+control_out)),
           ciL = log_plr -1.96*se_log_plr,
           ciU = log_plr +1.96*se_log_plr,
           plr_L = exp(ciL),
           plr_U = exp(ciU),
           z_score = (log_plr - 0) / se_log_plr,
           p_value = 2 * (1 - pnorm(abs(z_score)))) %>% 
    dplyr::rename(., {{ newname }} := "PLR") %>% 
    select(Gene, label, label2, family, all_of(newname)) %>% 
    return()
}

# Merge all gene level enrichment datasets
data_complete <- data_complete %>%
  left_join(read_delim(gene_level_enrichment_dataset1, delim = "\t") %>%
              process_df(.,"ClinVar one star + vs Regeneron")) %>% 
  left_join(read_delim(gene_level_enrichment_dataset2, delim = "\t") %>%
              process_df(.,"ClinVar pathogenic & HGMD vs gnomAD")) %>%
  left_join(read_delim(gene_level_enrichment_dataset3, delim = "\t") %>%
              process_df(.,"ClinVar pathogenic vs gnomAD")) %>%
  left_join(read_delim(gene_level_enrichment_dataset4, delim = "\t") %>%
              process_df(.,"HGMD vs gnomAD")) %>%
  left_join(read_delim(gene_level_enrichment_dataset5, delim = "\t") %>%
              process_df(.,"ClinVar one star + vs gnomAD"))

# Compute correlation matrix and p-values for patient subgroup
cor_test_results <- rcorr(as.matrix(
  data_complete %>% 
    filter(label == "patient") %>% 
    select(5:9) %>% 
    mutate(across(everything(), as.numeric)) %>%
    mutate(across(everything(), ~ replace(.x, is.infinite(.x), NA)))
), type = "pearson")

cor_matrix <- cor_test_results$r  # Correlation values
p_matrix <- cor_test_results$P %>% replace(is.na(.), 0)  # Replace NA p-values with 0

# Bonferroni correction for multiple comparisons
alpha <- c(0.05, 0.001)
bonferroni_threshold <- alpha / (ncol(cor_matrix) * (ncol(cor_matrix) - 1) / 2)

# Create a matrix of significance indicators
star_matrix <- ifelse(p_matrix < bonferroni_threshold | is.na(p_matrix), "**", "")

# Generate and save correlation plot
jpeg(figure_path, width = 8, height = 8, units = "in", res = 600)
corrplot(cor_matrix, method = "square", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         mar = c(0, 0, 2, 0),
         addCoef.col = NULL, # Hide correlation coefficients
         p.mat = p_matrix, sig.level = bonferroni_threshold, insig = "label_sig",
         pch.cex = 2, pch.col = "white")  # Add stars for significant correlations
dev.off()
