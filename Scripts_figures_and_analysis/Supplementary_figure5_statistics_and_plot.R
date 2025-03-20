### Load Required Libraries
library(tidyverse)
library(seqinr)
library(patchwork)  

### Define File Paths
clinvar_path <- ""
gene_family_path <- ""

figure_path <- ""

### Function: Extracts First Character of a Sequence
mod_aln <- function(gene_fam_select.df, i) {
  gene_fam_select.df %>% 
    select(all_of(i)) %>% 
    setNames("seq") %>% 
    mutate(seq = str_sub(seq, 1, 1)) %>% 
    pull(seq)
}

### Function: Counts Non-Gap Positions in a Sequence
count_pos <- function(seq) {
  count <- 0
  sapply(seq, function(i) {
    if (i == "-") return(0)
    count <<- count + 1
    return(count)
  })
}


### Function: Process Gene Alignment Data
process_aln <- function(name, gene_fam_select.df, i) {
  tibble(Gene = name, AA_ref = mod_aln(gene_fam_select.df, i)) %>% 
    mutate(Aln_pos = row_number(),
           AA_pos = count_pos(AA_ref))
}

### Load and Filter Variant Data
all_var.df <- read_delim(clinvar_path, delim = "\t") %>% 
  filter(Stars != "No star", AA_ref != AA_alt)

### Initialize Data Frame for Variant Classification Counts
variants.df <- tibble(
  label = c("One paralog", "Two paralogs", "Three paralogs", "Four or more paralogs"),
  Pathogenic = rep(0, 4),
  VUS = rep(0, 4),
  Benign = rep(0, 4)
)

### Process Each Gene Family
gene_families <- list.files(gene_family_path)
gene_fam_more_2 <- 0

for (i in seq_along(gene_families)) {
  gene_fam_select.df <- read_delim(file.path(gene_family_path, gene_families[i]), delim = "\t", show_col_types = FALSE)
  gene_index <- 2:ncol(gene_fam_select.df)
  genes_covered <- colnames(gene_fam_select.df)[2:ncol(gene_fam_select.df)]
  
  ### Process Alignments
  all_aln.df <- map_dfr(gene_index, function(gene_i) {
    process_aln(colnames(gene_fam_select.df)[gene_i], gene_fam_select.df, gene_i)
  }) %>% filter(AA_pos != 0)
  
  ### Function to filter Variants by Classification
  filter_variants <- function(classification) {
    all_var.df %>%
      filter(Gene %in% genes_covered) %>%
      left_join(all_aln.df, by = c("Gene", "AA_pos", "AA_ref")) %>%
      filter(!is.na(Aln_pos)) %>%
      distinct(Aln_pos, Classification, Gene) %>%
      mutate(cov = paste0(Classification,Gene)) %>% 
      filter(Classification == classification)
  }
  
  pathogenic_select.df <- filter_variants("Lp/P")
  vus_select.df <- filter_variants("VUS")
  benign_select.df <- filter_variants("Lb/B")
  
  path_genes <- unique(pathogenic_select.df$Gene)
  
  # Only considers genes/gene families where pathogenic variants have been observed in at least two genes
  if (length(path_genes) > 2) {
    gene_fam_more_2 <- gene_fam_more_2 + 1
    
    update_variant_counts <- function(df, col_name) {
      
      # Calculate number of variants
      update_var <- df %>% 
        filter(!is.na(cov), !is.na(Classification)) %>% 
        group_by(Aln_pos) %>%
        mutate(n = n()) %>% 
        ungroup() %>% 
        distinct(n,Aln_pos) %>% 
        mutate(label = case_when(n == 1 ~ "One paralog",
                                 n == 2 ~ "Two paralogs",
                                 n == 3 ~ "Three paralogs",
                                 n > 3 ~ "Four or more paralogs")) %>% 
        group_by(label) %>% 
        summarise(n = n()) %>% 
        mutate(n = replace_na(n, 0)) 
      
      variants.df <<- variants.df %>%
        left_join(update_var, by = "label") %>%
        mutate(n = replace_na(n, 0)) %>%
        mutate(!!sym(col_name) := !!sym(col_name) + n) %>%
        select(-n)
    }
    
    for (gene in path_genes) {
      patient_variant_data <- function(variant_df) {
        all_aln.df %>%
          filter(Gene == gene) %>%
          left_join(pathogenic_select.df %>% filter(Gene != gene) %>% select(cov, Aln_pos), by = "Aln_pos") %>%
          left_join(variant_df %>% filter(Gene == gene) %>% select(Classification, Aln_pos) %>% distinct(), by = "Aln_pos")
      }
      
      update_variant_counts(patient_variant_data(pathogenic_select.df), "Pathogenic")
      update_variant_counts(patient_variant_data(vus_select.df), "VUS")
      update_variant_counts(patient_variant_data(benign_select.df), "Benign")
      
    }
  }
}

# Generate Supplementary Figure 5 ####
plot_A <- tibble(
  label = rep(c("Lp/P variants vs VUS", "Lp/P variants vs Lb/B variants"), each = 4),
  paralogs = rep(variants.df$label, 2),
  ratio = c(variants.df$Pathogenic / variants.df$VUS, 
            variants.df$Pathogenic / variants.df$Benign),
  Total = c(variants.df$Pathogenic + variants.df$VUS, 
            variants.df$Pathogenic + variants.df$Benign)
) %>% 
  mutate(paralogs = case_when(paralogs == "One paralog"~"1",
                              paralogs == "Two paralogs"~"2",
                              paralogs == "Three paralogs"~"3",
                              paralogs == "Four or more paralogs"~"4+",
  )) %>% 
  filter(label == "Lp/P variants vs Lb/B variants") %>% 
  ggplot(aes(x = paralogs, y = ratio, fill = paralogs))+
  geom_bar(stat = "identity", position = position_dodge2(0.1))+
  geom_text(aes(y = ratio + 0.9, label = Total), size = 5, position = position_dodge2(0.9))+
  scale_fill_manual(values = c("#FDBB63","#FBA304","#FC6B63","#C92A2A"), name = "Number of pathogenic\nvariants in paralogous genes") +
  theme_bw(base_size = 25)+
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  labs(y = "Ratio of Lp/P to Lb/B variants",
       x = "Number of pathogenic variants\nin paralogous genes")

plot_B <- tibble(
  label = rep(c("Lp/P variants vs VUS", "Lp/P variants vs Lb/B variants"), each = 4),
  paralogs = rep(variants.df$label, 2),
  ratio = c(variants.df$Pathogenic / variants.df$VUS, 
            variants.df$Pathogenic / variants.df$Benign),
  Total = c(variants.df$Pathogenic + variants.df$VUS, 
            variants.df$Pathogenic + variants.df$Benign)
) %>% 
  mutate(paralogs = case_when(paralogs == "One paralog"~"1",
                              paralogs == "Two paralogs"~"2",
                              paralogs == "Three paralogs"~"3",
                              paralogs == "Four or more paralogs"~"4+",
  )) %>% 
  filter(label == "Lp/P variants vs VUS") %>% 
  ggplot(aes(x = paralogs, y = ratio, fill = paralogs))+
  geom_bar(stat = "identity", position = position_dodge2(0.1))+
  geom_text(aes(y = ratio + 0.05, label = Total), size = 5, position = position_dodge2(0.9))+
  scale_fill_manual(values = c("#FDBB63","#FBA304","#FC6B63","#C92A2A"), name = "Number of pathogenic\nvariants in paralogous genes") +
  theme_bw(base_size = 25)+
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  labs(y = "Ratio of Lp/P to VUS variants",
       x = "Number of pathogenic variants\nin paralogous genes")


Figure_complete <- (plot_A + plot_B) + 
  plot_layout(ncol = 2) + 
  plot_annotation(tag_levels = 'A')

# Save figure 
ggsave(figure_path, Figure_complete,  width = 13, height = 7.0, units = "in", dpi = 600)

