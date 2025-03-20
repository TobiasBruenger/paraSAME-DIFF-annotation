# Load required libraries
library(tidyverse)

# Define pathways to pathogenic dataset and figure save
pathogenic_dataset_path <- ""
figure_path <- ""
gene_family_path <- ""

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

# Load pathogenic variant dataset 
pathogenic.df <- rbind(read_delim(pathogenic_dataset_path, delim = "\t")) %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) # Ensure only missense variants

# Load gene families 
gene_families <- list.files(gene_family_path)

# Initialize output file 
output.df <- tibble(Number_of_paralogs = c(1:50),
                    Number_of_residues_covered_paralog_or_same_gene = rep(0,50),
                    Number_of_residues_covered_same_gene = rep(0,50))

# Count overlap of pathogenic paralogous variants per gene family size 
for(i in 1:length(gene_families)){ 
  
  gene_fam_select.df <- read_delim(paste0(gene_family_path,gene_families[i]), delim = "\t", show_col_types = FALSE)
  
  gene_index <- 2:ncol(gene_fam_select.df)
  
  genes_covered <- colnames(gene_fam_select.df)[gene_index]
  
  
  # Process alignment data for each gene
  all_aln.df <- map_dfr(gene_index, ~ process_aln(colnames(gene_fam_select.df)[.x], gene_fam_select.df, .x)) %>%
    filter(AA_pos != 0)  # Remove gaps
  
  # Select patient variants in gene family 
  patient_select.df <- pathogenic.df %>% 
    filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>% 
    left_join(all_aln.df, by = c("Gene" = "Gene","AA_pos" = "AA_pos","AA_ref" = "AA_ref")) %>% 
    filter(!is.na(Aln_pos)) %>% 
    mutate(cov = Gene) 
  
  path_genes <- unique(patient_select.df$Gene)

  if (length(path_genes) > 1) {
    merged.df <- left_join(all_aln.df, patient_select.df %>% distinct(cov, Aln_pos), by = "Aln_pos")
    
    num_paralog_or_same_gene_overlap <- merged.df %>% 
      filter(!is.na(cov)) %>% 
      distinct(Aln_pos,Gene) %>% nrow()
    
    num_same_gene_overlap <- merged.df %>% 
      filter(!is.na(cov),
             cov == Gene) %>% 
      distinct(Aln_pos,Gene) %>% nrow()
    
    output.df$Number_of_residues_covered_paralog_or_same_gene[length(gene_index)] <- output.df$Number_of_residues_covered_paralog_or_same_gene[length(gene_index)]+num_paralog_or_same_gene_overlap
    output.df$Number_of_residues_covered_same_gene[length(gene_index)] <- output.df$Number_of_residues_covered_same_gene[length(gene_index)]+num_same_gene_overlap
    
  }
}

#Prepare data for plot and correlation analysis 
Figure_input <- output.df %>% 
  #reduce number to gene families of size 2 to 15 family members
  filter(Number_of_paralogs>1,
         Number_of_paralogs<16) %>% 
  mutate(fold_change = Number_of_residues_covered_paralog_or_same_gene/Number_of_residues_covered_same_gene) 

# Generate figure 
Plot <- Figure_input %>% 
  ggplot(aes(x = Number_of_paralogs, y = fold_change))+
    geom_bar(stat = "identity")+
    theme_classic(base_size = 25)+
    labs(x = "Genes per gene family",
         y = "Additional residues covered by paralogous\npathogenic variants (Fold change)")+
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"))+
    coord_cartesian(expand = 0)

# Save figure 
ggsave(figure_path, Plot,  width = 10, height = 7.5, units = "in", dpi = 600)


#Correlation analysis 
cor.test(Figure_input$Number_of_paralogs, Figure_input$fold_change, method = "kendall")
