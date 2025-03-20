### Load Required Libraries
library(tidyverse)
library(seqinr)
library(ggpubr)

### Define File Paths
pathogenic_variants_paraz <- ""
pathogenic_variants_per <- ""
control_variants_paraz <- ""
control_variants_per <- ""
gene_family_path <- ""
path_to_of_per_accessible_genes <- ""
figure_path <- ""

# ===========================
# Functions for Alignment Processing
# ===========================

# Extract first amino acid from sequence
mod_aln <- function(gene_fam_select.df, i) {
  gene_fam_select.df %>%
    select(all_of(i)) %>%
    setNames("seq") %>%
    mutate(seq = str_sub(seq, 1, 1)) %>%
    pull(seq)  # Equivalent to .$seq for better clarity
}

# Count amino acid positions, ignoring gaps
count_pos <- function(seq) {
  count <- 0
  sapply(seq, function(i) {
    if (i == "-") return(0)
    count <<- count + 1
    return(count)
  })
}

# Process alignments and extract positional information
process_aln <- function(name, gene_fam_select.df, i) {
  tibble(
    Gene = name,
    AA_ref = mod_aln(gene_fam_select.df, i)
  ) %>%
    mutate(
      Aln_pos = row_number(),
      AA_pos = count_pos(AA_ref)
    )
}

# ===========================
# Function for Paralog Enrichment Analysis
# ===========================

#para-DIFF
paralog_enrichment_groupings_paraz_pos <- function(control_select.df, patient_select.df, 
                                             number_pat, number_con, label_sel, 
                                             paraz_max_min, paraz_max) {
  gene_selections <- unique(patient_select.df$Gene)
  out.df <- tibble()
  
  for (gen_sel in gene_selections) {
    for (paraz_c in c(0, paraz_max_min, paraz_max)) {
      
      # Define paralog score ranges
      paraz_low <- ifelse(paraz_c == 0, -1000, ifelse(paraz_c == paraz_max_min, 0.01, paraz_max))
      paraz_label <- c("neg paraz", "pos paraz", "max paraz")[match(paraz_c, c(0, paraz_max_min, paraz_max))]
      
      # Extract relevant control and patient variants
      control_gen.df <- control_select.df %>% 
        filter(Gene == gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c) %>% 
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt,Gene)
      
      pathogenic_gen.df <- patient_select.df %>% 
        filter(Gene == gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c) %>% 
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt, Gene)
      
      # Extract variants from paralogs
      control_other_gen.df <- control_select.df %>%
        filter(Gene != gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c, !is.na(Aln_pos)) %>%
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt, Gene)
      
      pathogenic_other_gen.df <- patient_select.df %>%
        filter(Gene != gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c, !is.na(Aln_pos)) %>%
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt, Gene)
      
      # Generate IDs for filtering
      if (number_pat == 1) {
        pathogenic_other_gen.df <- pathogenic_other_gen.df %>% 
          mutate(AA_ref_id = paste0(AA_ref,Aln_pos), AA_alt_id = paste0(Aln_pos,AA_alt))
      }
      if (number_con == 1) {
        control_other_gen.df <- control_other_gen.df %>% 
          mutate(AA_ref_id = paste0(AA_ref,Aln_pos), AA_alt_id = paste0(Aln_pos,AA_alt))
      }
      
      # Ensure only one type is selected
      if (number_con == 1 & number_pat == 1) {
        stop("Error: Both number_con and number_pat cannot be 1 simultaneously.")
      } else if (number_pat == 1){
        
        pathogenic_gen_para<- pathogenic_gen.df %>% 
          mutate(AA_ref_id = paste0(AA_ref,Aln_pos),
                 AA_alt_id = paste0(Aln_pos,AA_alt)) %>% 
          filter(AA_ref_id %in% pathogenic_other_gen.df$AA_ref_id,
                 !AA_alt_id %in% pathogenic_other_gen.df$AA_alt_id) %>% nrow()
        
        
        control_gen_para<-  control_gen.df %>% 
          mutate(AA_ref_id = paste0(AA_ref,Aln_pos),
                 AA_alt_id = paste0(Aln_pos,AA_alt)) %>% 
          filter(AA_ref_id %in% pathogenic_other_gen.df$AA_ref_id,
                 !AA_alt_id %in% pathogenic_other_gen.df$AA_alt_id) %>% nrow()
        
      } else if (number_con == 1){ 
        
        pathogenic_gen_para<- pathogenic_gen.df %>% 
          mutate(AA_ref_id = paste0(AA_ref,Aln_pos),
                 AA_alt_id = paste0(Aln_pos,AA_alt)) %>% 
          filter(AA_ref_id %in% control_other_gen.df$AA_ref_id,
                 !AA_alt_id %in% control_other_gen.df$AA_alt_id) %>% nrow()
        
        
        control_gen_para<-  control_gen.df %>% 
          mutate(AA_ref_id = paste0(AA_ref,Aln_pos),
                 AA_alt_id = paste0(Aln_pos,AA_alt)) %>% 
          filter(AA_ref_id %in% control_other_gen.df$AA_ref_id,
                 !AA_alt_id %in% control_other_gen.df$AA_alt_id) %>% nrow()
        
      } 
      
      # Store results
      out.df <- bind_rows(out.df, tibble(
        Gene = gen_sel,
        label = ifelse(number_pat == 1 & number_con == 0, "patient", "control"),
        label2 = label_sel,
        pathogenic_match = pathogenic_gen_para, 
        pathogenic_out = nrow(pathogenic_gen.df) - pathogenic_gen_para,
        control_match = control_gen_para, 
        control_out = nrow(control_gen.df) - control_gen_para,
        Paraz_cutoff = paraz_label
      ))
    }
  }
  return(out.df)
}

#para-SAME
paralog_enrichment_groupings_paraz_aaid <- function(control_select.df, patient_select.df, number_pat, number_con, label_sel, paraz_max_min, paraz_max) {
  
  gene_selections <- unique(patient_select.df$Gene)
  out.df <- tibble()
  
  for (gen_sel in gene_selections) {
    for (paraz_c in c(0, paraz_max_min, paraz_max)) {
      
      # Define paralog score ranges
      paraz_low <- ifelse(paraz_c == 0, -1000, ifelse(paraz_c == paraz_max_min, 0.01, paraz_max))
      paraz_label <- c("neg paraz", "pos paraz", "max paraz")[match(paraz_c, c(0, paraz_max_min, paraz_max))]
      
      # Extract relevant variants
      control_gen.df <- control_select.df %>%
        filter(Gene == gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c) %>%
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt, Gene)
      
      pathogenic_gen.df <- patient_select.df %>%
        filter(Gene == gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c) %>%
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt, Gene)
      
      # Extract variants from paralogs
      control_other_gen.df <- control_select.df %>%
        filter(Gene != gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c, !is.na(Aln_pos)) %>%
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt, Gene)
      
      pathogenic_other_gen.df <- patient_select.df %>%
        filter(Gene != gen_sel, Parazscore >= paraz_low, Parazscore <= paraz_c, !is.na(Aln_pos)) %>%
        distinct(Aln_pos, AA_pos, AA_ref, AA_alt, Gene)
      
      # Generate unique AA identifiers
      if (number_pat == 1) {
        pathogenic_other_gen.df <- pathogenic_other_gen.df %>% 
          mutate(AA_id = paste0(AA_ref, Aln_pos, AA_alt))
      }
      if (number_con == 1) {
        control_other_gen.df <- control_other_gen.df %>% 
          mutate(AA_id = paste0(AA_ref, Aln_pos, AA_alt))
      }
      
      # Ensure only one type is selected
      if (number_con == 1 & number_pat == 1) {
        stop("Error: Both number_con and number_pat cannot be 1 simultaneously.")
      }
      
      # Compute matching variants
      if (number_pat == 1) {
        pathogenic_gen_para <- pathogenic_gen.df %>%
          mutate(AA_id = paste0(AA_ref, Aln_pos, AA_alt)) %>%
          filter(AA_id %in% pathogenic_other_gen.df$AA_id) %>% nrow()
        
        control_gen_para <- control_gen.df %>%
          mutate(AA_id = paste0(AA_ref, Aln_pos, AA_alt)) %>%
          filter(AA_id %in% pathogenic_other_gen.df$AA_id) %>% nrow()
      } else if (number_con == 1) { 
        pathogenic_gen_para <- pathogenic_gen.df %>%
          mutate(AA_id = paste0(AA_ref, Aln_pos, AA_alt)) %>%
          filter(AA_id %in% control_other_gen.df$AA_id) %>% nrow()
        
        control_gen_para <- control_gen.df %>%
          mutate(AA_id = paste0(AA_ref, Aln_pos, AA_alt)) %>%
          filter(AA_id %in% control_other_gen.df$AA_id) %>% nrow()
      } else {
        stop("Error: Invalid input configuration.")
      }
      
      # Store results
      out.df <- bind_rows(out.df, tibble(
        Gene = gen_sel,
        label = ifelse(number_pat == 1 & number_con == 0, "patient", "control"),
        label2 = label_sel,
        pathogenic_match = pathogenic_gen_para, 
        pathogenic_out = nrow(pathogenic_gen.df) - pathogenic_gen_para,
        control_match = control_gen_para, 
        control_out = nrow(control_gen.df) - control_gen_para,
        Paraz_cutoff = paraz_label
      ))
    }
  }
  return(out.df)
}


# Variant loading #
pathogenic.df <- rbind(read_delim(pathogenic_variants_paraz, delim = "\t")) %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt)

control.df <- read_delim(control_variants_paraz,delim = "\t") %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt)

# Load genes covered in Lal et al. 2020 annotating Paraz scores 
Paraz_families<- read_delim(path_to_of_per_accessible_genes, delim= "\t")


gene_families <- list.files(gene_family_path) ##all gene families and their alignments 
all_pat_diff_ex.df <- tibble()
all_pat_same_ex.df <- tibble()
check_filt = 0 
fam_para_miss <- c()

for (gene_file in gene_families) {
  gene_fam_select.df <- read_delim(
    paste0(gene_family_path, gene_file),
    delim = "\t",
    col_types = cols(Adj_bin_count = "d", DM_adj_bin_count = "d", proxy = "c", `per.tag` = "c")
  )
  
  gene_index <- 2:(which(colnames(gene_fam_select.df) == "Parazscore") - 1)
  family_old <- paste(colnames(gene_fam_select.df[gene_index]), collapse = ";")
  family_filt <- Paraz_families %>% filter(Family_old == family_old, per_paraz_comp == TRUE)

  if (nrow(family_filt) >= 2) {
    all_aln.df <- bind_rows(lapply(gene_index, function(gene_i) {
      process_aln(colnames(gene_fam_select.df %>% select(all_of(gene_i))), gene_fam_select.df, gene_i)
    })) %>% filter(AA_pos != 0)
    
    all_aln.df <- all_aln.df %>%
      left_join(select(gene_fam_select.df, Index, Parazscore), by = c("Aln_pos" = "Index")) %>%
      filter(Gene %in% family_filt$Genes_old) %>%
      left_join(select(family_filt, Genes_old, Genes_new), by = c("Gene" = "Genes_old"))
    
    paraz_max_min <- sort(unique(all_aln.df$Parazscore), decreasing = TRUE)[2]
    paraz_max <- sort(unique(all_aln.df$Parazscore), decreasing = TRUE)[1]
    
    patient_select.df <- pathogenic.df %>%
      filter(Gene %in% family_filt$Genes_new) %>%
      left_join(all_aln.df, by = c("Gene" = "Genes_new", "AA_pos", "AA_ref"))
    
    control_select.df <- control.df %>%
      filter(Gene %in% family_filt$Genes_new) %>%
      left_join(all_aln.df, by = c("Gene" = "Genes_new", "AA_pos", "AA_ref"))
    
    if (length(unique(patient_select.df$Gene)) > 1 & !is.na(paraz_max)) {
      
      all_pat_diff_ex.df <- bind_rows(all_pat_diff_ex.df,
                                      paralog_enrichment_groupings_paraz_pos(control_select.df, patient_select.df, 0, 1, "Paralogous population\nvariants at same AP", paraz_max_min, paraz_max),
                                      paralog_enrichment_groupings_paraz_pos(control_select.df, patient_select.df, 1, 0, "Paralogous patient\nvariants at same AP", paraz_max_min, paraz_max))
      
      all_pat_same_ex.df <- bind_rows(all_pat_same_ex.df,
                                      paralog_enrichment_groupings_paraz_aaid(control_select.df, patient_select.df, 0, 1, "Paralogous population\nvariants at same AP", paraz_max_min, paraz_max),
                                      paralog_enrichment_groupings_paraz_aaid(control_select.df, patient_select.df, 1, 0, "Paralogous patient\nvariants at same AP", paraz_max_min, paraz_max))
    }
  }
  
}

para_diff_same_paraz_data <- bind_rows(all_pat_diff_ex.df %>%
                                     mutate(label2 = "Different exchange"),
                                   all_pat_same_ex.df %>%
                                     mutate(label2 = "Same exchange"))


# ===========================
# Visualization of Paralogous Variant Data
# ===========================

# Function to generate a plot based on exchange type
create_paralog_plot <- function(data, exchange_type, color_values, x_limit) {
  data %>%
    group_by(label2, Paraz_cutoff, label) %>%
    summarise(
      pathogenic_match = sum(pathogenic_match),
      pathogenic_out = sum(pathogenic_out),
      control_out = sum(control_out),
      control_match = sum(control_match),
      .groups = 'drop'
    ) %>%
    mutate(
      sens = pathogenic_match / (pathogenic_match + pathogenic_out),
      specf = control_out / (control_out + control_match),
      PLR = sens / (1 - specf),
      NLR = (1 - sens) / specf,
      log_plr = log(PLR),
      se_log_plr = sqrt(1 / pathogenic_match - 1 / (pathogenic_match + pathogenic_out) + 1 / control_match - 1 / (control_match + control_out)),
      ciL = log_plr - 1.96 * se_log_plr,
      ciU = log_plr + 1.96 * se_log_plr,
      plr_L = exp(ciL),
      plr_U = exp(ciU),
      z_score = (log_plr - 0) / se_log_plr,
      p_value = 2 * (1 - pnorm(abs(z_score)))
    ) %>% 
    filter(label2 == exchange_type) %>%
    mutate(
      label = ifelse(label == "control", "Control variants", "Pathogenic variants"),
      paraz_cut = case_when(
        Paraz_cutoff == "neg paraz" ~ "Parazscore < 0",
        Paraz_cutoff == "pos paraz" ~ "Parazscore > 0 &\nnot max Parazscore",
        TRUE ~ "Max Parazscore"
      ),
      paraz_cut = factor(paraz_cut, levels = c("Parazscore < 0", "Parazscore > 0 &\nnot max Parazscore", "Max Parazscore"))
    ) %>%
    ggplot(aes(x = PLR, y = paraz_cut, group = label)) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "grey", size = 1.4) +
    geom_point(aes(color = label), size = 6, position = position_dodge(width = 0.9), size = 2) +
    theme_bw(base_size = 20) +
    scale_color_manual(values = color_values, name = "Paralogous variants at same AP") +
    scale_x_continuous(limits = c(0.1, x_limit)) +
    geom_errorbar(aes(xmin = plr_L, xmax = plr_U, group = label), position = position_dodge(width = 0.9, preserve = "single"), width = 0.5) +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      axis.text = element_text(color = "black"),
      axis.ticks.x = element_line(color = "black", size = 1),
      axis.title = element_text(size = 19),
      axis.ticks.y = element_blank(),
      plot.title = element_text(color = "black", hjust = 0.5)
    ) +
    labs(x = "LR+", y = "", title = "") +
    guides(color = guide_legend(nrow = 4, byrow = TRUE))
}

# Generate sub-plots
plot_A <- create_paralog_plot(para_diff_same_paraz_data, "Same exchange", c("#00b0f0", "purple"), 30)
plot_B <- create_paralog_plot(para_diff_same_paraz_data, "Different exchange", c("#00b0f0", "#efbbff"), 15)

# Calculate average likelihhod odds ratio for para-SAME and para-DIFF with max paraz score and negative parax scores 

# para-SAME
para_diff_same_paraz_data %>%
  filter(label == "patient",
         label2 == "Same exchange",
         Paraz_cutoff == "max paraz") %>%
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


# ===========================
# Function to calculate number of protein positions covered by PER or a paralogous residue
# ===========================

pathogenic_per.df <- read_delim(pathogenic_variants_per, delim = "\t")

gene_families <- list.files(gene_family_path) 
number_per_overlap <- 0
number_paraCriteria_overlap <- 0
number_overlap_both <- 0
gene_fam_more_2 <- 0 #gene families with 2 or more affected genes 

for (i in seq_along(gene_families)) {
  
  # Load gene family data
  gene_fam_select.df <- read_delim(
    paste0(gene_family_path, gene_families[i]),
    delim = "\t",
    col_types = cols()
  )
  
  gene_index <- 2:(which(colnames(gene_fam_select.df) == "Parazscore") - 1)
  family_old <- paste(colnames(gene_fam_select.df[gene_index]), collapse = ";")
  
  family_filt <- Paraz_families %>%
    filter(Family_old == family_old, per_paraz_comp == TRUE)
  
  # Ensure at least two genes can be annotated with paralog scores
  if (nrow(family_filt) >= 2) {
    
    # Process alignments for each gene in the family
    all_aln.df <- bind_rows(lapply(gene_index, function(gene_i) {
      process_aln(colnames(gene_fam_select.df %>% select(all_of(gene_i))), gene_fam_select.df, gene_i)
    })) %>% filter(AA_pos != 0)
    
    all_aln.df <- all_aln.df %>% mutate(AA_ref = aaa(AA_ref))
    
    # Extract PER positions
    PER_pos <- gene_fam_select.df %>%
      filter(aa.per == "PER") %>%
      pull(Index)
    
    # Annotate patient variants used to calculate PER and now used to identify residues where para-SAME/DIFF is applicable 
    patient_annotated.df <- pathogenic_per.df %>%
      filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>%
      left_join(all_aln.df %>% select(Gene, AA_pos, AA_ref, Aln_pos), by = c("Gene", "AA_pos", "AA_ref")) %>%
      mutate(cov = "Yes") %>%
      distinct(cov, Aln_pos, AA_ref, Gene)
    
    # If patient annotation exists across more than one gene in a gene family 
    if (nrow(patient_annotated.df) > 1) {
      
      gene_fam_more_2 <- gene_fam_more_2 + 1
      
      paralog_cov <- 0
      paralog_and_per_cov <- 0
      
      for(gene in unique(all_aln.df$Gene)){
        
        # Identify residues that are covered by paralogous pathogenic variants in other genes
        paralog_cov <- paralog_cov + all_aln.df %>%
          left_join(patient_annotated.df %>% filter(Gene != gene) %>% select(Aln_pos, AA_ref,cov) %>% distinct(), by = c("Aln_pos", "AA_ref")) %>%
          filter(!is.na(cov)) %>%
          filter(Gene == gene) %>% 
          distinct(Aln_pos, Gene) %>% nrow()
        
        # Identify residues that are covered by paralogous pathogenic variants in other genes and a PER
        paralog_and_per_cov <- paralog_and_per_cov +
          all_aln.df %>%
          left_join(patient_annotated.df %>% filter(Gene != gene) %>% select(Aln_pos, cov) %>% distinct(), by = "Aln_pos") %>%
          filter(!is.na(cov)) %>%
          filter(Gene == gene) %>% 
          # In addition filter that the residue is in a PER
          filter(Aln_pos %in% PER_pos) %>% 
          distinct(Aln_pos, Gene) %>% nrow()
      }
      
      number_paraCriteria_overlap <- number_paraCriteria_overlap + paralog_cov
      number_per_overlap <- number_per_overlap + nrow(all_aln.df %>% filter(Aln_pos %in% PER_pos))
      number_overlap_both <- number_overlap_both + paralog_and_per_cov
    }
  }
}

# Visualize protein residue coverage
plot_C <- tibble(n = c(sum(number_per_overlap),number_paraCriteria_overlap,number_overlap_both),
                  label = c("PER","para-SAME/DIFF","PER &\npara-SAME/DIFF")) %>% 
  mutate(label = factor(label, levels = c("PER","PER &\npara-SAME/DIFF","para-SAME/DIFF"))) %>% 
  ggplot(aes(x = label, y = n, fill = label))+
  geom_bar(stat = "identity")+
  theme_classic(base_size = 20)+
  labs(x = "",
       y = "N residues covered")+
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 19),
        axis.ticks = element_line(color ="black"),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(color = "black", hjust = 0.5),
        legend.position = "none")+
  scale_y_continuous(limits = c(0,110000), breaks = c(0,25000,50000,75000,100000), expand = c(0,0), labels = c("0","25,000","50,000","75,000","100,000"))+
  scale_fill_manual(values = c("#ffa400","grey","purple"), name = "")


# ===========================
# Function for PER enrichment and para-DIFF/SAME on same variant set as used for PERs
# ===========================

# Function to calculate enrichment of pathogenic vs control variants at aligned positions
paralog_enrichment_groupings_pos <- function(control_select.df,control_select_per_data.df,patient_select.df,patient_select_per_data.df,number_pat,number_con,label_sel,family){
  
  out.df <- tibble()
  
  gene_selections <- unique(patient_select.df$Gene)
  
  for(gen_sel in gene_selections){
    
    control_gen.df <- control_select.df %>% 
      filter(Gene == gen_sel) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    pathogenic_gen.df <- patient_select.df %>% 
      filter(Gene == gen_sel)%>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    control_other_gen.df <- control_select_per_data.df %>% 
      filter(Gene != gen_sel) %>% 
      filter(!is.na(Aln_pos)) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    pathogenic_other_gen.df <- patient_select_per_data.df %>% 
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
paralog_enrichment_groupings_aaid <- function(control_select.df,control_select_per_data.df,patient_select.df,patient_select_per_data.df,number_pat,number_con,label_sel,family){
  
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
    
    control_other_gen.df <- control_select_per_data.df %>% 
      filter(Gene != gen_sel) %>% 
      filter(!is.na(Aln_pos)) %>% 
      distinct(Aln_pos,AA_pos,AA_ref,AA_alt,Gene)
    
    pathogenic_other_gen.df <- patient_select_per_data.df %>% 
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

# Load patient and control variant datasets

# Remove variants from novel variant dataset that have been used in the original PER dataset
pathogenic_new.df <- read_delim(pathogenic_variants_paraz, delim = "\t") %>%
  distinct() %>%
  filter(AA_ref != AA_alt) %>%
  mutate(AA_ref = aaa(AA_ref), AA_alt = aaa(AA_alt)) %>% 
  #Update gene names from old dataset
  left_join(Paraz_families %>% select(Genes_old, Genes_new) %>% distinct(), by = c("Gene" = "Genes_new")) %>% 
  select(-Gene) %>%
  rename(Gene = Genes_old) %>%
  mutate(ID_pos = paste0(Gene, ":", AA_ref, AA_pos)) %>% 
  filter(!(ID_pos %in% (pathogenic_per.df %>%  mutate(ID_pos = paste0(Gene, ":", AA_ref, AA_pos)) %>% .$ID_pos)))

# Load control variant datasets
control_per.df <- read_delim(control_variants_per, delim = "\t") %>% 
  mutate(ID_pos = paste0(Gene, ":", AA_ref, AA_pos))

control_new.df <- read_delim(control_variants_paraz, delim = "\t") %>%
  distinct() %>%
  filter(AA_ref != AA_alt) %>%
  mutate(AA_ref = aaa(AA_ref), AA_alt = aaa(AA_alt)) %>%
  left_join(Paraz_families %>% select(Genes_old, Genes_new) %>% distinct(), by = c("Gene" = "Genes_new")) %>%
  select(-Gene) %>%
  rename(Gene = Genes_old) %>%
  mutate(ID_pos = paste0(Gene, ":", AA_ref, AA_pos)) %>% 
  filter(!(ID_pos %in% control_per.df$ID_pos))

# Initialize output datafiles
all_per.df <- tibble()
all_pat_diff_ex_in.df <- tibble()
all_pat_same_ex_in.df <- tibble()

# Calculate pathogenic variant presence in target gene given the presence of PER or a paralogous pathogenic variants 
for (gene_file in gene_families) {
  gene_fam_select.df <- read_delim(
    paste0(gene_family_path, gene_file),
    delim = "\t", col_types = cols()
  )
  
  gene_index <- 2:(which(colnames(gene_fam_select.df) == "Parazscore") - 1)
  family_old <- paste(colnames(gene_fam_select.df[gene_index]), collapse = ";")
  
  family_filt <- Paraz_families %>%
    filter(Family_old == family_old, per_paraz_comp == TRUE)
  
  if (nrow(family_filt) >= 2) {
    all_aln.df <- bind_rows(lapply(gene_index, function(gene_i) {
      process_aln(colnames(gene_fam_select.df %>% select(all_of(gene_i))), gene_fam_select.df, gene_i)
    })) %>%
      filter(AA_pos != 0) %>%
      mutate(AA_ref = aaa(AA_ref))
    
    PER_pos <- gene_fam_select.df %>%
      filter(aa.per == "PER") %>%
      pull(Index)
    
    # Select patient and control variants
    patient_select_new.df <- pathogenic_new.df %>%
      filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>%
      left_join(all_aln.df, by = c("Gene", "AA_pos", "AA_ref"))
    
    control_select_new.df <- control_new.df %>%
      filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>%
      left_join(all_aln.df, by = c("Gene", "AA_pos", "AA_ref"))
    
    patient_select_old.df <- pathogenic_per.df %>%
      filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>%
      left_join(all_aln.df, by = c("Gene", "AA_pos", "AA_ref"))
    
    control_select_old.df <- control_per.df %>%
      filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>%
      left_join(all_aln.df, by = c("Gene", "AA_pos", "AA_ref"))
    
    family_sel <- paste(colnames(gene_fam_select.df %>% select(all_of(gene_index))), collapse = "_")
    
    for (gene_sel in unique(patient_select_new.df$Gene)) {
      if (length(PER_pos) > 0) {
        patient_select2_new.df <- patient_select_new.df %>% filter(Gene == gene_sel)
        control_select2_new.df <- control_select_new.df %>% filter(Gene == gene_sel)
        
        path_in <- sum(patient_select2_new.df$Aln_pos %in% PER_pos)
        path_out <- nrow(patient_select2_new.df) - path_in
        con_in <- sum(control_select2_new.df$Aln_pos %in% PER_pos)
        con_out <- nrow(control_select2_new.df) - con_in
        
        all_per.df <- bind_rows(all_per.df, tibble(
          Gene = gene_sel, family = family_sel,
          pathogenic_match = path_in,
          pathogenic_out = path_out,
          control_match = con_in,
          control_out = con_out
        ))
      }
    }
    
    if (length(unique(patient_select_new.df$Gene)) > 1) {
      all_pat_diff_ex_in.df <- bind_rows(all_pat_diff_ex_in.df,
                                         paralog_enrichment_groupings_pos(control_select_old.df, control_select_new.df, patient_select_old.df, patient_select_new.df, 1, 0, "Different exchange", family_sel))
      
      all_pat_same_ex_in.df <- bind_rows(all_pat_same_ex_in.df,
                                         paralog_enrichment_groupings_aaid(control_select_old.df, control_select_new.df, patient_select_old.df, patient_select_new.df, 1, 0, "Same exchange", family_sel))
    }
  }
}

# Function to calculate sensitivity, specificity, PLR, and NLR
calculate_metrics <- function(data) {
  data %>%
    summarise(
      pathogenic_match = sum(pathogenic_match),
      pathogenic_out = sum(pathogenic_out),
      control_out = sum(control_out),
      control_match = sum(control_match)
    ) %>%
    mutate(
      sens = pathogenic_match / (pathogenic_match + pathogenic_out),
      specf = control_out / (control_out + control_match),
      PLR = sens / (1 - specf),
      NLR = (1 - sens) / specf
    ) %>%
    mutate(
      log_plr = log(PLR),
      se_log_plr = sqrt(1 / pathogenic_match - 1 / (pathogenic_match + pathogenic_match) + 
                          1 / control_match - 1 / (control_match + control_out)),
      ciL = log_plr - 1.96 * se_log_plr,
      ciU = log_plr + 1.96 * se_log_plr,
      plr_L = exp(ciL),
      plr_U = exp(ciU),
      z_score = log_plr / se_log_plr,
      p_value = 2 * (1 - pnorm(abs(z_score)))
    )
}

# Calculate metrics for all_per.df dataset
plot_data_1 <- calculate_metrics(all_per.df) %>%
  mutate(g = "PER") %>%
  select(PLR, plr_L, plr_U, g)

# Calculate metrics for patient-related datasets
plot_data_2 <- bind_rows(all_pat_diff_ex_in.df, all_pat_same_ex_in.df) %>%
  filter(label == "patient") %>%
  group_by(label2) %>%
  calculate_metrics() %>%
  mutate(g = ifelse(label2 == "Same exchange", "para-SAME", "para-DIFF")) %>%
  select(PLR, plr_L, plr_U, g)

# Combine datasets and prepare for plotting
plot_data <- bind_rows(plot_data_1, plot_data_2) %>%
  arrange(PLR) %>%
  mutate(label = factor(g, levels = unique(g)))

# Create the plot
plot_D <- ggplot(plot_data, aes(x = PLR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey", size = 1.4) +
  geom_point(aes(color = label), size = 6, position = position_dodge(width = 0.9), size = 2) +
  geom_errorbar(aes(xmin = plr_L, xmax = plr_U, group = label), linewidth = 0.5, width = 0.2) +
  scale_color_manual(values = c("#ffa400","#efbbff" ,"purple"), name = "Method") +
  scale_x_continuous(limits = c(0.5, 10)) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.title = element_text(size = 19),
    axis.ticks.y = element_blank(),
    plot.title = element_text(color = "black", hjust = 0.5)
  ) +
  labs(
    x = "LR+",
    y = "",
    title = ""
  ) +
  guides(color = guide_legend(nrow = 4, byrow = TRUE))


# Combine plots 
plot <- ggarrange(
  plot_A,plot_B,
  plot_C,plot_D, ncol = 2, nrow = 2,
  labels = c("A","B","C","D"),
  font.label = list(size = 20))


# Save figure 
ggsave(figure_path, plot,  width = 15, height = 10, units = "in", dpi = 600)
