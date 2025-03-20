#03.20.2025
#The script considered a pre generated alignment file of a gene family of the user's choice as well as a custom set of pathogenic variants observed in genes of the selected gene family. Based on the alignment file and the list of pathogenic variants the script generates an output table where listing all protein residues where paraDIFF/SAME criteria could be applied when a novel variant is observed.  

#Load required packages 
library(tidyverse)
library(seqinr)

#Functions to process the mutiple sequence alignment 
mod_aln <- function(gene_fam_select.df,i){
  
  out <- gene_fam_select.df %>% 
    select(all_of(i)) %>% 
    setNames("seq") %>% 
    mutate(seq = str_sub(seq,1,1)) %>% 
    .$seq
  return(out)
}

count_pos <- function(seq){
  out <- c()
  count <- 0
  for(i in seq){
    
    if(i == "-"){
      out <- c(out,0)
    }else{
      count <- count+1
      out <- c(out,count)
    }
    
  }
  return(out)
}

process_aln <- function(name,gene_fam_select.df, i){
  
  out <- tibble(Gene = name, AA_ref = mod_aln(gene_fam_select.df,i)) %>% 
    mutate(Aln_pos = 1:nrow(.),
           AA_pos = count_pos(AA_ref))
  
  return(out)
  
}

#Load your set of pathogenic variants in the gene family of interest 
pathogenic_variants.df <- read_delim("Annotate_paralogous_pathogenic_variants/Input_Pathogenic_variants_example.txt",delim = "\t")

gene_fam_select.df <- read_delim("Annotate_paralogous_pathogenic_variants/Input_Gene_famlily_alignment_example.stats", delim = "\t")

gene_index <- 2:ncol(gene_fam_select.df)

all_aln.df <- tibble()
for(gene_i in gene_index){
  
  all_aln.df <- rbind(all_aln.df,process_aln(colnames(gene_fam_select.df %>% select(all_of(gene_i))), gene_fam_select.df,gene_i)) %>% 
    filter(AA_pos != 0) 
  
  
}

all_aln.df <- all_aln.df %>% 
  mutate(AA_ref = aaa(AA_ref))

#select variants to create gnomad and patient df's
patient_variants_annotated.df <- pathogenic_variants.df %>% 
  filter(Gene %in% colnames(gene_fam_select.df %>% select(all_of(gene_index)))) %>% 
  left_join(all_aln.df, by = c("Gene" = "Gene","AA_pos" = "AA_pos","AA_ref" = "AA_ref")) %>% 
  filter(!is.na(Aln_pos)) %>% 
  #left_join(transcripts.df, by = c("Gene" = "Gene_ID")) %>% 
  #mutate(ID = paste0(Gene,"(", Transcript, "):",AA_ref,AA_pos,AA_alt)) %>% 
  mutate(IDgene = Gene,
         ID = paste0(Gene,":",AA_ref,AA_pos,AA_alt)) %>% 
  select(Aln_pos,AA_ref,ID,IDgene) %>% 
  group_by(Aln_pos,AA_ref,IDgene) %>% 
  summarise(ID = paste(ID,collapse = ";")) %>%  #merges mutiple variants with different amino acid exchanges observed at the same protein position
  distinct()

all_aln.df %>% 
  left_join(patient_variants_annotated.df, by = c("Aln_pos" = "Aln_pos","AA_ref" = "AA_ref"), relationship = "many-to-many") %>% 
  filter(!is.na(ID)) %>% 
  filter(IDgene != Gene) %>% 
  select(Gene,AA_pos,AA_ref,Aln_pos,ID) %>% 
  dplyr::rename(
    `Aminoacid position`= "AA_pos",
    `Reference aminoacid` = "AA_ref",
    `Gene-family alignment position` = "Aln_pos",
    `Paralogous pathogenic variants at same alignment position` = "ID"
  ) %>% View()
  write_delim("Annotate_paralogous_pathogenic_variants/Output_Annotated_alignment_example.txt",delim = "\t")

