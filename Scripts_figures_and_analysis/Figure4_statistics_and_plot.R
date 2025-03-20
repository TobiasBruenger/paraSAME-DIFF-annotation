### Load Required Libraries
library(tidyverse)
library(seqinr)
library(bio3d)
library(Rfast)
library(corrplot)
library(readxl)

#Figure and file paths 
scn_variants_w_phenotype <- ""
scn_regeneron.df <- ""

scn_gene_alignment <- ""
nav1_2_structure_data <- ""

figure_path1 <- ""
figure_path2 <- ""

# Load variants 
scn_patient.df <- read_delim(scn_variants_w_phenotype, delim = "\t")

# Load alignment 
all_scn.aln <- read.alignment(scn_gene_alignment, format = "clustal") 
all_scn_aln.df <- rbind(process_aln("SCN11A",1),
                        process_aln("SCN10A",2),
                        process_aln("SCN5A",3),
                        process_aln("SCN4A",4),
                        process_aln("SCN8A",5),
                        process_aln("SCN9A",6),
                        process_aln("SCN3A",7),
                        process_aln("SCN1A",8),
                        process_aln("SCN2A",9)) %>% 
  filter(AA_pos != 0)

# Function to add unique identifier 
add_unique_i <- function(data.df){
  
  i <- 0
  data.df$index <- NA
  for(pos in unique(data.df$AA_pos)){
    
    i <- i+1
    data.df$index[which(data.df$AA_pos == pos)] <- i
    
  }
  return(data.df$index)
}

# Set seed for reproducibility
set.seed(23)

# Initialize empty tibble
out_all.df <- tibble()

# Define phenotype order
pheno_order <- scn_patient.df$pheno_id %>% unique() %>% sort()
pheno_order <- pheno_order[c(3:22,1,2)] # Rearrange order

# Prepare patient phenotype dataset
patient_3d_input.df <- scn_patient.df %>%
  distinct() %>%
  group_by(pheno_id, Aln_pos, SCN2A_pos) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(id_cols = c(Aln_pos, SCN2A_pos), names_from = pheno_id, values_from = count) %>%
  replace(is.na(.), 0) %>%
  mutate(SCN2A_pos = ifelse(SCN2A_pos == 0, NA, SCN2A_pos)) %>%
  select(Aln_pos, 3:24, SCN2A_pos)

# Load structural data
scn2a_struc.df <- read_delim(nav1_2_structure_data, delim = "\t") %>%
  filter(chain == "A") %>%
  mutate(inter_res = NA) %>%
  filter(!is.na(Position_in_structure)) %>%
  left_join(patient_3d_input.df, by = c("Position_in_structure" = "SCN2A_pos"))

# Load PDB file and extract atom information
scn2a.pdb <- read.pdb("6j8e")$atom %>%
  as_tibble() %>%
  filter(chain == "A", resno < 2006)

# Define cutoff distance for residue interaction
cutoff <- 5

# Compute distance matrix for residue interactions
dist.mat <- Dist(scn2a.pdb %>% select(x, y, z) %>% as.matrix())

# Identify interacting residues
for (res in unique(scn2a.pdb$resno)) {
  select_atoms <- which(scn2a.pdb$resno == res)
  in_cutoff_atoms <- which(matrixStats::rowMins(dist.mat[, select_atoms]) < cutoff)
  selected_res <- unique(scn2a.pdb$resno[in_cutoff_atoms])
  
  scn2a_struc.df$inter_res[scn2a_struc.df$Position_in_structure == res] <- 
    paste(selected_res, collapse = ";")
}

# Aggregate phenotype data by interacting residues
scn2a_struc_corr_input.df <- scn2a_struc.df %>%
  replace(is.na(.), 0) %>%
  group_by(Position_in_structure, Uniprot_position, inter_res) %>%
  summarise(across(starts_with("SCN"), sum), .groups = "drop") %>%
  select(4:25, Position_in_structure, Uniprot_position, inter_res)

scn2a_struc_bubble.df <- scn2a_struc_corr_input.df

# Compute mean phenotype values for interacting residues
for (pheno in pheno_order) {
  for (i in 1:nrow(scn2a_struc_corr_input.df)) {
    inter_res <- str_split(scn2a_struc_corr_input.df$inter_res[i], ";") %>%
      unlist() %>%
      as.numeric()
    
    scn2a_struc_bubble.df[i, pheno] <- mean(
      scn2a_struc_corr_input.df[scn2a_struc_corr_input.df$Position_in_structure %in% inter_res, pheno] %>%
        unlist()
    )
  }
}

# Initialize correlation matrices
pheno_cor.mat <- matrix(nrow = length(pheno_order), ncol = length(pheno_order))
rownames(pheno_cor.mat) <- pheno_order
colnames(pheno_cor.mat) <- pheno_order

pheno_cor_p.mat <- matrix(nrow = length(pheno_order), ncol = length(pheno_order))
rownames(pheno_cor_p.mat) <- pheno_order
colnames(pheno_cor_p.mat) <- pheno_order

# Compute correlation matrix using Kendall's correlation
for (pheno1 in pheno_order) {
  for (pheno2 in pheno_order) {
    cor_out <- cor.test(
      scn2a_struc_bubble.df[[pheno1]], 
      scn2a_struc_bubble.df[[pheno2]], 
      method = "kendall"
    )
    
    pheno1_i <- which(pheno_order == pheno1)
    pheno2_i <- which(pheno_order == pheno2)
    
    pheno_cor.mat[pheno1_i, pheno2_i] <- cor_out$estimate
    pheno_cor_p.mat[pheno1_i, pheno2_i] <- cor_out$p.value
  }
}

# Define color palette
col_sel <- colorRampPalette(c("#FFC000", "white", "#7000FE"))(200)

# Generate  and save correlation plot 
jpeg(figure_path1, width = 8, height = 8, units = "in", res = 600, quality = 100)
corrplot(
  pheno_cor.mat, 
  method = "color", 
  type = "lower", 
  insig = 'label_sig', 
  p.mat = pheno_cor_p.mat, 
  sig.level = c(0.001 / 230, 0.05 / 230), 
  pch.cex = 1.1, 
  tl.col = "black", 
  order = "original", 
  tl.cex = 1, 
  cl.cex = 1.1, 
  col = col_sel
)
dev.off()

# Generate second subfigure ####

# ---------------------------
# PREPARATION & INITIALIZATION
# ---------------------------

# Initialize empty output tibble
out.df <- tibble()

# Function to create stratified sets for each phenotype
create_stratified_sets <- function(data, set_sizes) {
  set.seed(23) # Ensure reproducibility
  random_set <- sample(max(data$index))
  split(random_set, cut(seq_along(random_set), length(set_sizes), labels = FALSE))
}

# Prepare datasets for different phenotypes
prepare_dataset <- function(pheno_name) {
  scn_patient.df %>%
    arrange(Aln_pos) %>%
    filter(pheno_id == pheno_name) %>%
    mutate(index = add_unique_i(.))
}

# Define phenotype datasets and create stratified sets
scn2a_dee.df <- prepare_dataset("SCN2A: Early onset DEE")
set_2A <- create_stratified_sets(scn2a_dee.df, c(11, 11, 11, 10))

scn1a_dravet.df <- prepare_dataset("SCN1A: Dravet")
set_1A <- create_stratified_sets(scn1a_dravet.df, c(62, 62, 62, 60))

scn5a_brugada.df <- prepare_dataset("SCN5A: Brugada Syndrome")
set_5A <- create_stratified_sets(scn5a_brugada.df, c(33, 33, 31, 33))

scn8a_dee.df <- prepare_dataset("SCN8A: DEE")
set_8A <- create_stratified_sets(scn8a_dee.df, c(20, 20, 20, 19))

# Load gnomAD dataset
control.df <- read_delim("Data_analysis/Genome_biology_2025/Datasets/Variants/regeneron_no_gnomAD_final.txt", delim = "\t") %>%
  left_join(all_scn_aln.df %>% distinct(Gene, AA_pos, AA_ref, Aln_pos),
            by = c("Gene", "AA_pos", "AA_ref")) %>%
  mutate(pathogenicity = "Control") %>%
  arrange(Aln_pos) %>%
  mutate(index = add_unique_i(.)) %>%
  select(-Transcript) %>%
  filter(Gene %in% c("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A")) %>%
  distinct()

set_control<- create_stratified_sets(control.df, c(484, 483, 483, 483))

# ---------------------------
# CROSS-VALIDATION LOOP (4-fold)
# ---------------------------
for(split in 1:4){
  
  # Display current split
  print(paste("Processing split", split))
  
  # ---------------------------
  # SET UP TRAINING & TEST SETS
  # ---------------------------
  # Define patient training data (exclude certain phenotypes)
  patient_train.df <- scn_patient.df %>% 
    filter(!(pheno_id %in% c("SCN1A: Dravet", 
                             "SCN2A: Early onset DEE", 
                             "SCN5A: Brugada Syndrome", 
                             "SCN8A: DEE")))
  
  # Initialize empty patient test data tibble
  patient_test.df <- tibble()
  
  # For each dataset (set_1A, set_2A, set_5A, set_8A, set_control),
  # define training and test indices for the current split.
  set_1A_train <- unlist(set_1A)[ !(unlist(set_1A) %in% set_1A[[split]]) ]
  set_1A_test  <- unlist(set_1A)[  (unlist(set_1A) %in% set_1A[[split]]) ]
  
  set_2A_train <- unlist(set_2A)[ !(unlist(set_2A) %in% set_2A[[split]]) ]
  set_2A_test  <- unlist(set_2A)[  (unlist(set_2A) %in% set_2A[[split]]) ]
  
  set_5A_train <- unlist(set_5A)[ !(unlist(set_5A) %in% set_5A[[split]]) ]
  set_5A_test  <- unlist(set_5A)[  (unlist(set_5A) %in% set_5A[[split]]) ]
  
  set_8A_train <- unlist(set_8A)[ !(unlist(set_8A) %in% set_8A[[split]]) ]
  set_8A_test  <- unlist(set_8A)[  (unlist(set_8A) %in% set_8A[[split]]) ]
  
  set_control_train <- unlist(set_control)[ !(unlist(set_control) %in% set_control[[split]]) ]
  set_control_test  <- unlist(set_control)[  (unlist(set_control) %in% set_control[[split]]) ]
  
  # Append phenotype-specific patient data to training set
  patient_train.df <- patient_train.df %>% 
    rbind( 
      scn1a_dravet.df %>% filter(index %in% set_1A_train) %>% select(-index),
      scn2a_dee.df    %>% filter(index %in% set_2A_train) %>% select(-index),
      scn5a_brugada.df %>% filter(index %in% set_5A_train) %>% select(-index),
      scn8a_dee.df    %>% filter(index %in% set_8A_train) %>% select(-index)
    )
  
  # Append phenotype-specific patient data to test set
  patient_test.df <- patient_test.df %>% 
    rbind(
      scn1a_dravet.df %>% filter(index %in% set_1A_test) %>% select(-index),
      scn2a_dee.df    %>% filter(index %in% set_2A_test) %>% select(-index),
      scn5a_brugada.df %>% filter(index %in% set_5A_test) %>% select(-index),
      scn8a_dee.df    %>% filter(index %in% set_8A_test) %>% select(-index)
    )
  
  # Split control data into training and test sets
  control_test.df  <- control.df %>% filter(index %in% set_control_test) %>% select(-index)
  control_train.df <- control.df %>% filter(index %in% set_control_train) %>% select(-index)
  
  # ---------------------------
  # PREPARE 3D PATIENT INPUT DATA
  # ---------------------------
  # Define phenotype order (custom ordering)
  pheno_order <- scn_patient.df$pheno_id %>% unique() %>% sort()
  pheno_order <- pheno_order[c(3:22, 1, 2)]
  
  # Build the 3D input data from patient training data:
  # 1. Remove duplicate rows.
  # 2. Group by pheno_id, alignment position, and SCN2A position.
  # 3. Summarize counts and pivot to a wide format.
  # 4. Replace missing values with 0 and adjust SCN2A_pos.
  patient_3d_input.df <- patient_train.df %>%
    distinct() %>% 
    group_by(pheno_id, Aln_pos, SCN2A_pos) %>% 
    summarise(count = sum(count), .groups = "drop") %>% 
    pivot_wider(id_cols = c(Aln_pos, SCN2A_pos), 
                names_from = pheno_id, 
                values_from = count) %>% 
    replace(is.na(.), 0) %>% 
    mutate(SCN2A_pos = ifelse(SCN2A_pos == 0, NA, SCN2A_pos)) %>%
    select(Aln_pos, 3:24, SCN2A_pos)
  
  # ---------------------------
  # LOAD STRUCTURAL DATA & COMPUTE RESIDUE INTERACTIONS
  # ---------------------------
  # Load structural data and merge with 3D input based on position
  scn2a_struc.df <- read_delim("Data_analysis/ClinGen_criteria/6j8e_without_variants.txt", 
                               delim = "\t") %>% 
    filter(chain == "A") %>% 
    mutate(inter_res = NA) %>% 
    filter(!is.na(Position_in_structure)) %>% 
    left_join(patient_3d_input.df, by = c("Position_in_structure" = "SCN2A_pos"))
  
  # Load PDB file and extract atom information for chain A with residue number < 2006
  scn2a.pdb <- read.pdb("6j8e") %>% .$atom %>% 
    as_tibble() %>% 
    filter(chain == "A", resno < 2006)
  
  # Set cutoff distance for interactions
  cutoff <- 5
  # Compute distance matrix using x,y,z coordinates
  dist.mat <- Dist(scn2a.pdb %>% select(x, y, z) %>% as.matrix())
  
  # For each unique residue, find all atoms within the cutoff and record interacting residues
  for(res in unique(scn2a.pdb$resno)){
    select_atoms <- which(scn2a.pdb$resno == res) 
    in_cutoff_atoms <- which(dist.mat[, select_atoms] %>% matrixStats::rowMins() < cutoff)
    selected_res <- unique(scn2a.pdb$resno[in_cutoff_atoms])
    scn2a_struc.df$inter_res[ scn2a_struc.df$Position_in_structure == res ] <- 
      paste(unique(selected_res), collapse = ";") %>% unlist()
  }
  
  # ---------------------------
  # PREPARE STRUCTURE CORRELATION INPUT
  # ---------------------------
  # Group by position and compute sum for each phenotype column;
  # then select key columns for further analysis.
  scn2a_struc_corr_input.df <- scn2a_struc.df %>% 
    replace(is.na(.), 0) %>% 
    group_by(Position_in_structure, Uniprot_position, inter_res) %>% 
    summarise(`SCN3A: DEE`                = sum(`SCN3A: DEE`),
              `SCN8A: DEE`                = sum(`SCN8A: DEE`),
              `SCN8A: Benign Epilepsy`     = sum(`SCN8A: Benign Epilepsy`),
              `SCN8A: NNDwoE`             = sum(`SCN8A: NNDwoE`),
              `SCN1A: Dravet`             = sum(`SCN1A: Dravet`),
              `SCN1A: GEFS+`              = sum(`SCN1A: GEFS+`),
              `SCN1A: Early onset DEE`    = sum(`SCN1A: Early onset DEE`),
              `SCN1A: FHM3`               = sum(`SCN1A: FHM3`),
              `SCN2A: Benign Epilepsy`     = sum(`SCN2A: Benign Epilepsy`),
              `SCN2A: Early onset DEE`    = sum(`SCN2A: Early onset DEE`),
              `SCN2A: Late onset DEE`     = sum(`SCN2A: Late onset DEE`),
              `SCN2A: Autism`             = sum(`SCN2A: Autism`),
              `SCN4A: Relaxation Impairment Disorders` = sum(`SCN4A: Relaxation Impairment Disorders`),
              `SCN4A: Myopathy`           = sum(`SCN4A: Myopathy`),
              `SCN4A: Hypokalemic PP`     = sum(`SCN4A: Hypokalemic PP`),
              `SCN4A: Hyperkalemic PP`     = sum(`SCN4A: Hyperkalemic PP`),
              `SCN5A: Brugada Syndrome`   = sum(`SCN5A: Brugada Syndrome`),
              `SCN5A: Long QT Syndrome`   = sum(`SCN5A: Long QT Syndrome`),
              `SCN9A: Small-Fiber Neuropathy` = sum(`SCN9A: Small-Fiber Neuropathy`),
              `SCN9A: Paroxysmal Pain Disorders` = sum(`SCN9A: Paroxysmal Pain Disorders`),
              `SCN11A: Neuropathies`      = sum(`SCN11A: Neuropathies`),
              `SCN10A: Neuropathies`      = sum(`SCN10A: Neuropathies`),
              .groups = "drop") %>% 
    ungroup() %>% 
    select(4:25, Position_in_structure, Uniprot_position, inter_res)
  
  # Copy the correlation input to a bubble data frame for further processing
  scn2a_struc_bubble.df <- scn2a_struc_corr_input.df
  
  # ---------------------------
  # CALCULATE MEAN VALUES ACROSS INTERACTING RESIDUES
  # ---------------------------
  # For each phenotype in pheno_order, update each row by taking the mean value
  # over the interacting residues listed in the "inter_res" field.
  for(pheno in pheno_order){
    print(pheno)
    for(i in 1:nrow(scn2a_struc_corr_input.df)){
      # Split the inter_res string into numeric residue positions
      inter_res <- scn2a_struc_corr_input.df$inter_res[i] %>% 
        str_split(";", simplify = FALSE) %>% unlist() %>% as.numeric()
      # Compute the mean for the phenotype over the interacting residues
      scn2a_struc_bubble.df[i, pheno] <- mean(
        scn2a_struc_corr_input.df[ scn2a_struc_corr_input.df$Position_in_structure %in% inter_res, pheno ] %>% unlist()
      )
    }
  }
  
  # ---------------------------
  # COMPUTE CORRELATION MATRICES
  # ---------------------------
  # Initialize empty matrices for correlation coefficients and p-values
  pheno_cor.mat <- matrix(nrow = length(pheno_order), ncol = length(pheno_order))
  rownames(pheno_cor.mat) <- pheno_order
  colnames(pheno_cor.mat) <- pheno_order
  
  pheno_cor_p.mat <- matrix(nrow = length(pheno_order), ncol = length(pheno_order))
  rownames(pheno_cor_p.mat) <- pheno_order
  colnames(pheno_cor_p.mat) <- pheno_order
  
  # Calculate Kendall correlation and associated p-values for each pair of phenotypes
  for(pheno1 in pheno_order){
    for(pheno2 in pheno_order){
      cor_out <- cor.test(
        scn2a_struc_bubble.df[[pheno1]] %>% unlist(),
        scn2a_struc_bubble.df[[pheno2]] %>% unlist(), 
        method = "kendall"
      )
      pheno1_i <- which(pheno_order == pheno1)
      pheno2_i <- which(pheno_order == pheno2)
      pheno_cor.mat[pheno1_i, pheno2_i] <- cor_out$estimate
      pheno_cor_p.mat[pheno1_i, pheno2_i] <- cor_out$p.value
    }
  }
  
  # ---------------------------
  # IDENTIFY SIMILAR & DIFFERENT DISEASE DISTRIBUTIONS
  # ---------------------------
  similar_diseases <- list()
  different_diseases <- list()
  
  # For each phenotype, determine which diseases are significantly positively or negatively correlated.
  for(pheno_sel in pheno_order){
    similar_diseases[pheno_sel] <- pheno_cor.mat[ , colnames(pheno_cor.mat) == pheno_sel] %>% 
      as_tibble() %>% 
      mutate(diseases = colnames(pheno_cor.mat),
             p = pheno_cor_p.mat[ , colnames(pheno_cor_p.mat) == pheno_sel]) %>% 
      filter(diseases != pheno_sel) %>% 
      filter(p < 0.05/230, value > 0) %>% .$diseases %>% list()
    
    different_diseases[pheno_sel] <- pheno_cor.mat[ , colnames(pheno_cor.mat) == pheno_sel] %>% 
      as_tibble() %>% 
      mutate(diseases = colnames(pheno_cor.mat),
             p = pheno_cor_p.mat[ , colnames(pheno_cor_p.mat) == pheno_sel]) %>% 
      filter(diseases != pheno_sel) %>% 
      filter(value < 0) %>% .$diseases %>% list()
  }
  
  # ---------------------------
  # ENRICHMENT ANALYSIS
  # ---------------------------
  out_split.df <- tibble()
  
  # Loop through the three similarity groups: "similar", "different", "control"
  for(sim in c("similar", "different", "control")){
    
    # For enrichment analysis, consider only the four specified phenotypes
    pheno_order <- c("SCN1A: Dravet", "SCN2A: Early onset DEE", "SCN5A: Brugada Syndrome", "SCN8A: DEE")
    
    for(pheno_sel in pheno_order){
      gene_sel <- sub(":.*", "", pheno_sel)
      
      if(sim == "similar"){
        scn_path_sel.df <- patient_test.df %>% 
          filter(Gene == gene_sel, pheno_id == pheno_sel)
        scn_path_other.df <- patient_train.df %>% 
          filter(Gene != gene_sel, pheno_id %in% similar_diseases[[pheno_sel]])
        scn_con_sel.df <- control_test.df %>% filter(Gene == gene_sel)
        
        in_path  <- scn_path_sel.df %>% filter(Aln_pos %in% scn_path_other.df$Aln_pos) %>% nrow()
        out_path <- scn_path_sel.df %>% filter(!(Aln_pos %in% scn_path_other.df$Aln_pos)) %>% nrow()
        in_con   <- scn_con_sel.df %>% filter(Aln_pos %in% scn_path_other.df$Aln_pos) %>% nrow()
        out_con  <- scn_con_sel.df %>% filter(!(Aln_pos %in% scn_path_other.df$Aln_pos)) %>% nrow()
        
        sens     <- in_path/(in_path + out_path)
        specf    <- out_con/(out_con + in_con)
        PLR      <- sens/(1 - specf)
        log_plr  <- log(PLR)
        se_log_plr <- sqrt(1/in_path - 1/(in_path + out_path) + 1/in_con - 1/(in_con + out_con))
        ciL      <- log_plr - 1.96 * se_log_plr
        ciU      <- log_plr + 1.96 * se_log_plr
        plr_L    <- exp(ciL)
        plr_U    <- exp(ciU)
        z_score  <- (log_plr - 0) / se_log_plr
        p_value  <- 2 * (1 - pnorm(abs(z_score)))
        
        out_split.df <- rbind(out_split.df,
                              tibble(label = "test",
                                     split_sel = split,
                                     gene = gene_sel,
                                     similarity = sim,
                                     pred_pheno = pheno_sel,
                                     PLR = PLR,
                                     pvalue = p_value,
                                     lCI = plr_L,
                                     uCI = plr_U,
                                     sens = sens,
                                     specf = specf,
                                     in_pat = in_path,
                                     in_con = in_con,
                                     out_pat = out_path,
                                     out_con = out_con))
        
      } else if(sim == "different"){
        scn_path_sel.df <- patient_test.df %>% 
          filter(Gene == gene_sel, pheno_id == pheno_sel)
        scn_path_other.df <- patient_train.df %>% 
          filter(Gene != gene_sel, pheno_id %in% different_diseases[[pheno_sel]])
        scn_con_sel.df <- control_test.df %>% filter(Gene == gene_sel)
        
        in_path  <- scn_path_sel.df %>% filter(Aln_pos %in% scn_path_other.df$Aln_pos) %>% nrow()
        out_path <- scn_path_sel.df %>% filter(!(Aln_pos %in% scn_path_other.df$Aln_pos)) %>% nrow()
        in_con   <- scn_con_sel.df %>% filter(Aln_pos %in% scn_path_other.df$Aln_pos) %>% nrow()
        out_con  <- scn_con_sel.df %>% filter(!(Aln_pos %in% scn_path_other.df$Aln_pos)) %>% nrow()
        
        sens     <- in_path/(in_path + out_path)
        specf    <- out_con/(out_con + in_con)
        PLR      <- sens/(1 - specf)
        log_plr  <- log(PLR)
        se_log_plr <- sqrt(1/in_path - 1/(in_path + out_path) + 1/in_con - 1/(in_con + out_con))
        ciL      <- log_plr - 1.96 * se_log_plr
        ciU      <- log_plr + 1.96 * se_log_plr
        plr_L    <- exp(ciL)
        plr_U    <- exp(ciU)
        z_score  <- (log_plr - 0) / se_log_plr
        p_value  <- 2 * (1 - pnorm(abs(z_score)))
        
        out_split.df <- rbind(out_split.df,
                              tibble(label = "test",
                                     split_sel = split,
                                     gene = gene_sel,
                                     similarity = sim,
                                     pred_pheno = pheno_sel,
                                     PLR = PLR,
                                     pvalue = p_value,
                                     lCI = plr_L,
                                     uCI = plr_U,
                                     sens = sens,
                                     specf = specf,
                                     in_pat = in_path,
                                     in_con = in_con,
                                     out_pat = out_path,
                                     out_con = out_con))
        
      } else if(sim == "control"){
        scn_path_sel.df <- patient_test.df %>%
          filter(Gene == gene_sel, pheno_id == pheno_sel)
        scn_con_other.df <- control_train.df %>% filter(Gene != gene_sel)
        scn_con_sel.df   <- control_test.df %>% filter(Gene == gene_sel)
        
        in_path  <- scn_path_sel.df %>% filter(Aln_pos %in% scn_con_other.df$Aln_pos) %>% nrow()
        out_path <- scn_path_sel.df %>% filter(!(Aln_pos %in% scn_con_other.df$Aln_pos)) %>% nrow()
        in_con   <- scn_con_sel.df %>% filter(Aln_pos %in% scn_con_other.df$Aln_pos) %>% nrow()
        out_con  <- scn_con_sel.df %>% filter(!(Aln_pos %in% scn_con_other.df$Aln_pos)) %>% nrow()
        
        sens     <- in_path/(in_path + out_path)
        specf    <- out_con/(out_con + in_con)
        PLR      <- sens/(1 - specf)
        log_plr  <- log(PLR)
        se_log_plr <- sqrt(1/in_path - 1/(in_path + out_path) + 1/in_con - 1/(in_con + out_con))
        ciL      <- log_plr - 1.96 * se_log_plr
        ciU      <- log_plr + 1.96 * se_log_plr
        plr_L    <- exp(ciL)
        plr_U    <- exp(ciU)
        z_score  <- (log_plr - 0) / se_log_plr
        p_value  <- 2 * (1 - pnorm(abs(z_score)))
        
        out_split.df <- rbind(out_split.df,
                              tibble(label = "test",
                                     gene = gene_sel,
                                     split_sel = split,
                                     similarity = sim,
                                     pred_pheno = pheno_sel,
                                     PLR = PLR,
                                     pvalue = p_value,
                                     lCI = plr_L,
                                     uCI = plr_U,
                                     sens = sens,
                                     specf = specf,
                                     in_pat = in_path,
                                     in_con = in_con,
                                     out_pat = out_path,
                                     out_con = out_con))
      }
    }
  }
  
  # Append the results from this split to the overall output
  out_all.df <- rbind(out_all.df, out_split.df)
}

# Group the overall results (out_all.df) by gene, similarity, and predicted phenotype.
# Then, calculate the total counts for patient in/out and control in/out.
out_all_pr.df <- out_all.df %>% 
  group_by(gene, similarity, pred_pheno) %>% 
  summarise(
    in_pat  = sum(in_pat),
    out_pat = sum(out_pat),
    in_con  = sum(in_con),
    out_con = sum(out_con),
    .groups = "drop"
  )


# COMPUTE PERFORMANCE METRICS PER GROUP
# Initialize an empty tibble to store the final results.
out_final.df <- tibble()

# Loop over each row in the aggregated dataframe to calculate performance metrics.
for(i in 1:nrow(out_all_pr.df)){
  
  # Extract counts for the current group
  in_path  <- out_all_pr.df$in_pat[i]
  out_path <- out_all_pr.df$out_pat[i]
  in_con   <- out_all_pr.df$in_con[i]
  out_con  <- out_all_pr.df$out_con[i]
  
  # Calculate sensitivity and specificity
  sens  <- in_path / (in_path + out_path)
  specf <- out_con / (out_con + in_con)
  
  # Compute the Positive Likelihood Ratio (PLR)
  PLR <- sens / (1 - specf)
  
  # Calculate the natural logarithm of PLR and its standard error
  log_plr    <- log(PLR)
  se_log_plr <- sqrt( 1/in_path - 1/(in_path + out_path) +
                        1/in_con   - 1/(in_con   + out_con) )
  
  # Determine the 95% confidence interval on the log-scale and back-transform
  ciL   <- log_plr - 1.96 * se_log_plr
  ciU   <- log_plr + 1.96 * se_log_plr
  plr_L <- exp(ciL)
  plr_U <- exp(ciU)
  
  # Compute z-score and corresponding two-sided p-value for log(PLR)
  z_score <- log_plr / se_log_plr
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  # Append the calculated metrics along with grouping variables to the output dataframe
  out_final.df <- rbind(
    out_final.df,
    tibble(
      label      = "test",
      gene       = out_all_pr.df$gene[i],
      similarity = out_all_pr.df$similarity[i],
      pred_pheno = out_all_pr.df$pred_pheno[i],
      PLR        = PLR,
      pvalue     = p_value,
      lCI        = plr_L,
      uCI        = plr_U,
      sens       = sens,
      specf      = specf,
      in_pat     = in_path,
      in_con     = in_con,
      out_pat    = out_path,
      out_con    = out_con
    )
  )
}

# FIGURE GENERATION: ENRICHMENT ANALYSIS PLOT
plot2 <- out_final.df %>% 
  mutate(pred_pheno = pred_pheno) %>% 
  arrange(desc(PLR)) %>% 
  mutate(similarity = case_when(
    similarity == "similar" ~ "Correlated",
    similarity == "different" ~ "Not correlated",
    similarity == "control" ~ "Control"
  )) %>% 
  mutate(
    similarity = factor(similarity, levels = c("Correlated", "Not correlated", "Control")),
    pred_pheno = factor(pred_pheno, levels = unique(pred_pheno))
  ) %>% 
  mutate(text = paste0("OR: ", round(PLR, 2))) %>% 

  ggplot(aes(x = PLR, y = pred_pheno, color = similarity)) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey", size = 1.4) +
  geom_point(size = 6, position = position_dodge(0.8)) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("#7000FE", "#e69f00", "darkgrey"), name = "") +
  scale_x_continuous(trans = "log10", breaks = c(0.1, 1, 10, 100)) +
  geom_errorbar(aes(xmin = lCI, xmax = uCI), position = position_dodge(0.8), width = 0.2) +
  scale_y_discrete(limits = rev) +
  theme(
    panel.border = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5)
  ) +
  labs(x = "LR+", y = "", title = "") +
  coord_cartesian(xlim = c(0.1, 1000)) +
  geom_hline(yintercept = seq(1.5, 18.5, by = 1), linetype = "dashed", size = 1, color = "grey", alpha = 0.5)


# Save figure 
ggsave(figure_path2, plot2,  width = 13, height = 6, units = "in", dpi = 600)
