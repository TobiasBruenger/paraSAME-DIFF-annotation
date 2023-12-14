#14.12.2023
#The script considers the output table generated with "Annotate_paralogous_variants.R" script together with a list of novel variants. The script checks for which of those variants one/mutiple paralogous pathogenic variants are present.

#Load required packages 
library(tidyverse)


annnotated_alignment.df <- read_delim("Annotate_paralogous_pathogenic_variants/Output_Annotated_alignment_exmaple.txt",delim = "\t") %>% 
  mutate(Paralogous_prot = str_split(`Paralogous pathogenic variants at same alignment position`,";")) %>% 
  unnest(Paralogous_prot) %>% 
  mutate(Paralogous_AA_alt = str_sub(Paralogous_prot,-3,-1))

novel_variants.df <- read_delim("Check_paraPS_PM5_application/Input_example_variants.txt",delim = "\t")


for (i in 1:nrow(novel_variants.df)) {
  gene_match <- novel_variants.df$Gene[i] %in% annnotated_alignment.df$Gene
  
  if (gene_match) {
    gene_index <- which(novel_variants.df$Gene[i] == annnotated_alignment.df$Gene)
    annnotated_alignment_sel_gene.df <- annnotated_alignment.df[gene_index,]
    
    pos_match <- novel_variants.df$`Aminoacid position`[i] %in% annnotated_alignment_sel_gene.df$`Aminoacid position`
    
    if (pos_match) {
      pos_index <- which(novel_variants.df$`Aminoacid position`[i] == annnotated_alignment_sel_gene.df$`Aminoacid position`)
      annnotated_alignment_sel_pos.df <- annnotated_alignment_sel_gene.df[pos_index,]
      
      aa_alt_match <- novel_variants.df$AA_alt[i] %in% annnotated_alignment_sel_pos.df$Paralogous_AA_alt
      
      if (aa_alt_match) {
        novel_variants.df$Evidence_from_paralogs[i] <- "Paralogous pathogenic variant at same position with the same amino acid substitution (para-PS1)"
      } else {
        novel_variants.df$Evidence_from_paralogs[i] <- "Paralogous pathogenic variant at same position with a different amino acid substitution (para-PM5)"
      }
    } else {
      novel_variants.df$Evidence_from_paralogs[i] <- "None"
    }
  } else {
    novel_variants.df$Evidence_from_paralogs[i] <- "None"
  }
}



write_delim(novel_variants.df,"Check_paraPS_PM5_application/Output_example_variants_with_paralogous_evidences_annotated.txt",delim = "\t")
