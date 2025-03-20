# Load required libraries
library(tidyverse)
library(seqinr)
library(ggpubr)
library(rstatix)

#Define pathways to save Figures, to the pathogenic and control datasets and gene family alignments 
figure_path <- ""
pathogenic_dataset_path <- ""
control_dataset_path <- ""
gene_family_path <- ""


# Function to extract first character from aligned sequences
mod_aln <- function(gene_fam_select.df,i){
  
  out <- gene_fam_select.df %>% 
    select(all_of(i)) %>% 
    setNames("seq") %>% 
    mutate(seq = str_sub(seq,1,1)) %>% 
    .$seq
  return(out)
}

# Function to count positions in alignment while ignoring gaps ('-')
count_pos <- function(seq) {
  count <- 0
  sapply(seq, function(i) {
    if (i == "-") 0 else {
      count <<- count + 1
      count
    }
  })
}

# Function to process sequence alignment, returning reference amino acids and positions
process_aln <- function(name,gene_fam_select.df, i){
  
  out <- tibble(Gene = name, AA_ref = mod_aln(gene_fam_select.df,i)) %>% 
    mutate(Aln_pos = 1:nrow(.),
           AA_pos = count_pos(AA_ref))
  
  return(out)
  
}

# Function to calculate enrichment of pathogenic vs control variants at aligned positions
paralog_enrichment_groupings_pos <- function(control_select.df,patient_select.df,number_pat,number_con,label_sel,family){
  
  out.df <- tibble()
  
  gene_selections <- unique(patient_select.df$Gene)

  for(gen_sel in gene_selections){
    
    control_gen.df <- control_select.df %>% 
      filter(Gene == gen_sel) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    pathogenic_gen.df <- patient_select.df %>% 
      filter(Gene == gen_sel)%>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    control_other_gen.df <- control_select.df %>% 
      filter(Gene != gen_sel) %>% 
      filter(!is.na(Aln_pos)) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    pathogenic_other_gen.df <- patient_select.df %>% 
      filter(Gene != gen_sel) %>% 
      filter(!is.na(Aln_pos)) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    if(number_pat == 1){
      pathogenic_other_gen.df <- pathogenic_other_gen.df %>% 
        group_by(Aln_pos, AA_alt, AA_ref) %>% 
        mutate(AA_id = paste0(Aln_pos,AA_alt),
               AA_ref_id = paste0(Aln_pos, AA_ref))
      
      pred_variants <- n_distinct(pathogenic_other_gen.df$AA_ref_id)
      
    }
    
    if(number_con == 1){
      control_other_gen.df <- control_other_gen.df %>% 
        group_by(Aln_pos, AA_alt,AA_ref) %>% 
        mutate(AA_id = paste0(Aln_pos,AA_alt),
               AA_ref_id = paste0(Aln_pos, AA_ref))
      
      pred_variants <- n_distinct(control_other_gen.df$AA_ref_id)
      
    }
    
    if(number_con == 1 & number_pat == 1){ 
      stop("Wrong input: Both number_pat and number_con cannot be 1 simultaneously.")
    }
    else if (number_pat == 1){
      
      pathogenic_gen_para<- pathogenic_gen.df %>% 
        mutate(AA_id = paste0(Aln_pos,AA_alt),
               AA_ref_id = paste0(Aln_pos,AA_ref)) %>% 
        filter(!(AA_id %in% pathogenic_other_gen.df$AA_id),
               (AA_ref_id %in% pathogenic_other_gen.df$AA_ref_id)) %>% nrow()
      
      
      control_gen_para<-  control_gen.df %>% 
        mutate(AA_id = paste0(Aln_pos,AA_alt),
               AA_ref_id = paste0(Aln_pos,AA_ref)) %>% 
        filter(!(AA_id %in% pathogenic_other_gen.df$AA_id),
               (AA_ref_id %in% pathogenic_other_gen.df$AA_ref_id)) %>% nrow()
      
    } else if (number_con == 1){ 
      
      pathogenic_gen_para<- pathogenic_gen.df %>% 
        mutate(AA_id = paste0(Aln_pos,AA_alt),
               AA_ref_id = paste0(Aln_pos,AA_ref)) %>% 
        filter(!(AA_id %in% control_other_gen.df$AA_id),
               (AA_ref_id %in% control_other_gen.df$AA_ref_id)) %>% nrow()
      
      
      control_gen_para<-  control_gen.df %>% 
        mutate(AA_id = paste0(Aln_pos,AA_alt),
               AA_ref_id = paste0(Aln_pos,AA_ref)) %>% 
        filter(!(AA_id %in% control_other_gen.df$AA_id),
               (AA_ref_id %in% control_other_gen.df$AA_ref_id)) %>% nrow()
      
    } else{
      stop("Wrong input")
    }
    
    out.df <- rbind(out.df,tibble(Gene = gen_sel, 
                                  label = ifelse(number_pat == 1 & number_con == 0, "patient","control"),
                                  label2 = label_sel,
                                  family = family,
                                  pathogenic_match = pathogenic_gen_para, 
                                  pathogenic_out = nrow(pathogenic_gen.df)- pathogenic_gen_para,
                                  control_match = control_gen_para, 
                                  control_out = nrow(control_gen.df)- control_gen_para,
                                  pred_vars = pred_variants))
    
  }
  
  return(out.df)
  
}

# Function to calculate enrichment with exact amino acid matches in paralogs 
paralog_enrichment_groupings_aaid <- function(control_select.df,patient_select.df,number_pat,number_con,label_sel,family){
  
  out.df <- tibble()
  
  gene_selections <- unique(patient_select.df$Gene)
  
  out.df <- tibble()
  
  for(gen_sel in gene_selections){
    
    control_gen.df <- control_select.df %>% 
      filter(Gene == gen_sel) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    pathogenic_gen.df <- patient_select.df %>% 
      filter(Gene == gen_sel)%>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    control_other_gen.df <- control_select.df %>% 
      filter(Gene != gen_sel) %>% 
      filter(!is.na(Aln_pos)) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    pathogenic_other_gen.df <- patient_select.df %>% 
      filter(Gene != gen_sel) %>% 
      filter(!is.na(Aln_pos)) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    if(number_pat == 1){
      pathogenic_other_gen.df <- pathogenic_other_gen.df %>% 
        mutate(AA_id = paste0(AA_ref,Aln_pos,AA_alt))
      
      pred_variants <- n_distinct(pathogenic_other_gen.df$AA_id)
      
    }
    
    if(number_con == 1){
      control_other_gen.df <- control_other_gen.df %>% 
        mutate(AA_id = paste0(AA_ref,Aln_pos,AA_alt))
      
      pred_variants <- n_distinct(control_other_gen.df$AA_id)
      
    }
    
    if(number_con == 1 & number_pat == 1){ #both variants present 
      stop("Wrong input: Both number_pat and number_con cannot be 1 simultaneously.")
    }else if (number_pat == 1){
      
      pathogenic_gen_para<- pathogenic_gen.df %>% 
        mutate(AA_id = paste0(AA_ref,Aln_pos,AA_alt)) %>% 
        filter(AA_id %in% pathogenic_other_gen.df$AA_id) %>% nrow()
      
      
      control_gen_para<-  control_gen.df %>% 
        mutate(AA_id = paste0(AA_ref,Aln_pos,AA_alt)) %>% 
        filter(AA_id %in% pathogenic_other_gen.df$AA_id) %>% nrow()
      
    } else if (number_con == 1){ 
      
      pathogenic_gen_para<- pathogenic_gen.df %>% 
        mutate(AA_id = paste0(AA_ref,Aln_pos,AA_alt)) %>% 
        filter(AA_id %in% control_other_gen.df$AA_id) %>% nrow()
      
      
      control_gen_para<-  control_gen.df %>% 
        mutate(AA_id = paste0(AA_ref,Aln_pos,AA_alt)) %>% 
        filter(AA_id %in% control_other_gen.df$AA_id) %>% nrow()
      
    } else{
      stop("Wrong input")
    }
    
    out.df <- rbind(out.df,tibble(Gene = gen_sel, 
                                  label = ifelse(number_pat == 1 & number_con == 0, "patient","control"),
                                  label2 = label_sel,
                                  family = family,
                                  pathogenic_match = pathogenic_gen_para, 
                                  pathogenic_out = nrow(pathogenic_gen.df)- pathogenic_gen_para,
                                  control_match = control_gen_para, 
                                  control_out = nrow(control_gen.df)- control_gen_para,
                                  pred_vars = pred_variants))
  }
  
  return(out.df)
  
}

# Load and prepare variant datasets

# Load dataset of pathogenic variants - example file provided 
pathogenic.df <- rbind(read_delim(pathogenic_dataset_path, delim = "\t")) %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) # Ensure only missense variants

# Load dataset of control variants -example file provided 
control.df <- read_delim(control_dataset_path,delim = "\t") %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) # Ensure only missense variants


# List of gene families and their alignments
gene_families <- list.files(gene_family_path) 

# Initialize storage for enrichment results
all_pat_diff_ex_in.df <- tibble()
all_pat_same_ex_in.df <- tibble()

# Process each gene family alignment
for(i in 1:length(gene_families)){ 
  
  gene_fam_select.df <- read_delim(paste0(gene_family_path,gene_families[i]), delim = "\t", show_col_types = FALSE)
  
  gene_index <- 2:ncol(gene_fam_select.df)
  
  
  all_aln.df <- map_dfr(gene_index, ~process_aln(colnames(gene_fam_select.df)[.x], gene_fam_select.df, .x)) %>%
    filter(AA_pos != 0)  # Remove gaps
  

  # Select variants for patient and control groups
  patient_select.df <- pathogenic.df %>% 
    filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>% 
    left_join(all_aln.df, by = c("Gene" = "Gene","AA_pos" = "AA_pos","AA_ref" = "AA_ref")) %>% 
    filter(!is.na(Aln_pos))
  
  control_select.df <- control.df %>% 
    filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>% 
    left_join(all_aln.df, by = c("Gene" = "Gene","AA_pos" = "AA_pos","AA_ref" = "AA_ref"))%>% 
    filter(!is.na(Aln_pos))
  
  family <- paste(colnames(gene_fam_select.df)[gene_index], collapse = "_")

  # Check for at least two affected paralogous genes
  if(length(patient_select.df$Gene %>% unique) >1){ # patient variants in at least two paralogous genes 
    
    all_pat_diff_ex_in.df <- rbind(all_pat_diff_ex_in.df,
                                   paralog_enrichment_groupings_pos(control_select.df,patient_select.df,0,1,"Paralogous population\nvariants at same AP",family),
                                   paralog_enrichment_groupings_pos(control_select.df,patient_select.df,1,0,"Paralogous patient\nvariants at same AP",family))
    
    all_pat_same_ex_in.df <- rbind(all_pat_same_ex_in.df,
                                   paralog_enrichment_groupings_aaid(control_select.df,patient_select.df,0,1,"Paralogous population\nvariants at same AP",family),
                                   paralog_enrichment_groupings_aaid(control_select.df,patient_select.df,1,0,"Paralogous patient\nvariants at same AP",family))
    
  }
}

# Combine results from para-SAME and para-DIFF
paralog_evidences_combined.df <- bind_rows(all_pat_diff_ex_in.df %>% 
                                         mutate(label2 = "Different exchange"),
                                          all_pat_same_ex_in.df %>% 
                                         mutate(label2 = "Same exchange"))


# Initialize counters for pathogenic variant coverage
num_paralog_overlap <- 0  # Residues covered by pathogenic variants in a paralogous gene
num_all_overlap <- 0  # Residues covered by pathogenic variants in the same or a paralogous gene
num_same_gene_overlap <- 0  # Residues covered by pathogenic variants in the same gene
num_pathogenic_var_in_gene <- 0 # Number of pathogenic variants in genes of gene families with >1 gene with pathogenic variants 
num_pathogenic_var_in_paralogs <- 0 # Number of paralogous pathogenic variants in genes gene families with >1 gene with pathogenic variants 
gene_fam_more_2 <- 0  # Gene families with 2+ affected genes where para-Same/Diff is applicable
genes_fam <- 0

for(i in 1:length(gene_families)){ #length(gene_families)
  
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
    gene_fam_more_2 <- gene_fam_more_2 + 1
    
    genes_fam <- genes_fam + length(gene_index)
    
    # Get total numnber of pathogenic variants observed in across all genes 
    num_pathogenic_var_in_gene <- num_pathogenic_var_in_gene + nrow(patient_select.df)
    
    # Get total number of paralogous pathogenic variants covering any other gene of the gene family 
    paralog_vars.df <- left_join(all_aln.df, patient_select.df %>% select(cov, Aln_pos,AA_ref), by = c("Aln_pos","AA_ref"), relationship = "many-to-many")
    num_pathogenic_var_in_paralogs <- num_pathogenic_var_in_paralogs+
      paralog_vars.df %>% 
      filter(!is.na(cov),
             cov != Gene) %>% nrow()
    
    merged.df <- left_join(all_aln.df, patient_select.df %>% distinct(cov, Aln_pos,AA_ref), by = c("Aln_pos","AA_ref"))
    
    num_paralog_overlap <- num_paralog_overlap + 
      merged.df %>% 
      filter(!is.na(cov),
             cov != Gene) %>% 
      distinct(Aln_pos,Gene) %>% nrow()
    
    num_same_gene_overlap <- num_same_gene_overlap + 
      merged.df %>% 
      filter(!is.na(cov),
             cov == Gene) %>% 
      distinct(Aln_pos,Gene) %>% nrow()
    
    num_all_overlap <- num_all_overlap + 
      merged.df %>% 
      filter(!is.na(cov)) %>% 
      distinct(Aln_pos,Gene) %>% nrow()

  }
}

# Generate plot to demonstrate overlap of pathogenic variants in the same gene and overlap of pathogenic paralogous variants 
plot_A <- tibble(n = c(num_same_gene_overlap,num_paralog_overlap, num_all_overlap),
              label = c("Same gene","Paralogous gene", "Same or paralogous gene")) %>% 
  mutate(label = factor(label, levels = c("Same gene","Paralogous gene","Same or paralogous gene"))) %>% 
  ggplot(aes(x = label, y = n, fill = label))+
  geom_bar(stat = "identity")+
  theme_classic(base_size = 20)+
  labs(x = "",
       y = "Number of residues with\n pathogenic variants")+
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color ="black"),
        legend.position = "None")+
  scale_x_discrete(labels = c("Same\ngene","Paralogous\ngene", "Same or\nparalogous gene"))+
  scale_y_continuous(limits = c(0,90000), breaks = c(0,25000,50000,75000), expand = c(0,0), labels = c("0","25,000","50,000","75,000"))+
  scale_fill_manual(values = c("grey","#ffa400","#ff6961"))   


# Prepare dataset for plotting
plot_B_input <- paralog_evidences_combined.df %>% 
  # Filter for genes with at least 10 pathogenic and control variants and at least one overlapping variant
  filter((pathogenic_match + pathogenic_out) >9,
         (control_out + control_match) >9,
         (pathogenic_match + control_match) >0) %>%
  # Calculate likelihood ratios
  mutate(sens = pathogenic_match/(pathogenic_match + pathogenic_out),
         specf = control_out /(control_out + control_match),
         PLR = sens/(1-specf),
         NLR = (1-sens)/specf) %>% 
  # Compute confidence intervals and statistical values
  mutate(log_plr = log(PLR),
         se_log_plr = sqrt(1/pathogenic_match - 1/(pathogenic_match+pathogenic_match) + 1/control_match - 1/(control_match+control_out)),
         ciL = log_plr -1.96*se_log_plr,
         ciU = log_plr +1.96*se_log_plr,
         plr_L = exp(ciL),
         plr_U = exp(ciU),
         z_score = (log_plr - 0) / se_log_plr,
         p_value = 2 * (1 - pnorm(abs(z_score)))) %>% 
  # Scale extreme values
  mutate(PLR = pmin(pmax(PLR, 0.1), 1000))

# Generate likelihood ratio plot
plot_B <- plot_B_input %>% 
  mutate(label_final = ifelse(label == "control","Control","Pathogenic")) %>% 
  ggboxplot(., x = "label2", y = "PLR", fill = "label_final", alpha =0.5, outlier.shape = NA)+
  geom_point(aes(color = label2, group = label),alpha = 0.5,   position=position_jitterdodge(jitter.width = 0.2), size = 0.8)+
  geom_hline(yintercept = 1, linetype = "dashed",color = "black")+
  theme_classic(base_size = 20)+
  labs(y = "LR+",
       x = "Amino acid exchange",
       title = "")+
  scale_fill_manual(values = c("#00b0f0","purple","#00b0f0","purple"), name = "Paralogous variant:")+
  scale_color_manual(values = c("black","black"))+
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 19),
        plot.title = element_text(color = "black", hjust = 0.5),
        legend.position = "bottom")+
  scale_y_log10(limits = c(0.1,500), breaks = c(0.1,1,10,100))+
  guides(color = "none")


# Perform Wilcoxon test to for differences between paralogous pathogenic and control variants 
stat.test <- plot_B_input %>%
  group_by(label2) %>%
  wilcox_test(PLR ~ label) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "label2", dodge = 0.8)

# Add statistical annotations to plot
plot_B <- plot_B + stat_pvalue_manual(stat.test,label = "p", y.position = 2.3)

# Generate final plot 6.5x14
figure2 <- ggarrange(plot_A,plot_B, ncol = 2, nrow = 1,
          labels = c("A","B"),
          font.label = list(size = 20))

# Save Figure 
ggsave(figure_path,plot = figure2, width = 14, height = 6.5, units = "in", dpi = 600)


# Calculate average likelihhod odds ratio for para-SAME and para-DIFF

# para-SAME
paralog_evidences_combined.df %>%
  filter(label == "patient",
         label2 == "Same exchange") %>%
  summarise(pathogenic_match = sum(pathogenic_match),
            pathogenic_out = sum(pathogenic_out),
            control_out = sum(control_out),
            control_match = sum(control_match)) %>%
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
         p_value = 2 * (1 - pnorm(abs(z_score))))

# para-DIFF
paralog_evidences_combined.df %>%
  filter(label == "patient",
         label2 == "Same exchange") %>%
  summarise(pathogenic_match = sum(pathogenic_match),
            pathogenic_out = sum(pathogenic_out),
            control_out = sum(control_out),
            control_match = sum(control_match)) %>%
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
         p_value = 2 * (1 - pnorm(abs(z_score))))

