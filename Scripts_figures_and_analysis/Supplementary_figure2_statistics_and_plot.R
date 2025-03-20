### Load Required Libraries
library(tidyverse)
library(seqinr)
library(patchwork)  

### Define File Paths
clivnar_pathogenic_one_star_path <- ""
clinvar_pathogenic_all_path <- ""
hgmd_path <- ""
gnomad_path <- ""
regeneron_no_gnomad_path <- ""

# ClinVar dataset containg Lp/P, VUS and Lb/B variant calssifications 
clinvar_all_class <- ""

figure_path <- ""

#Load all genes with pathogenic variants 
genes <- read_delim(clinvar_pathogenic_all_path, delim = "\t") %>% 
  rbind(read_delim(hgmd_path, delim = "\t")) %>% 
  distinct(Gene)

# Load individual variant datasets
clinvar_variants <- read_delim(clinvar_pathogenic_all_path, delim = "\t") %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) %>% # Ensure only missense variants
  nrow()

clinvar_one_star_variants <- read_delim(clivnar_pathogenic_one_star_path, delim = "\t") %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) %>% # Ensure only missense variants
  nrow()

hgmd_variants <- read_delim(hgmd_path, delim = "\t") %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) %>% # Ensure only missense variants
  nrow()

gnomad_variants <- read_delim(gnomad_path, delim = "\t") %>% 
  distinct() %>% 
  filter(Gene %in% genes$Gene) %>% 
  filter(AA_ref != AA_alt) %>% # Ensure only missense variants
  nrow()

Regeneron_no_gnomAD_variants <- read_delim(regeneron_no_gnomad_path, delim = "\t") %>% 
  filter(Gene %in% genes$Gene) %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) %>% # Ensure only missense variants
  nrow()

#Supplementary Figure 4 ####

plot_A <- tibble(label = c("ClinVar Lp/P \n(all)",
                           "ClinVar Lp/P \n(>= 1 review gold star",
                           "HGMD",
                           "gnomAD",
                           "Regeneron \n(not overlapping with gnomAD)"),
                 n = c(clinvar_variants, clinvar_one_star_variants, hgmd_variants, gnomad_variants, Regeneron_no_gnomAD_variants),
                 status = c("Pathogenic","Pathogenic","Pathogenic","Control","Control")) %>% 
  mutate(label= factor(label, levels = c("ClinVar Lp/P \n(all)",
                                         "ClinVar Lp/P \n(>= 1 review gold star",
                                         "HGMD",
                                         "gnomAD",
                                         "Regeneron \n(not overlapping with gnomAD)"))) %>% 
  ggplot(aes(x = label, y = n, fill = status))+
  geom_bar(stat = "identity", position = position_dodge2(0.1))+
  geom_text(aes(y = n + 200000, label = n), size = 5, position = position_dodge2(0.9))+
  scale_fill_manual(values = c("#00b0f0","purple"), name = "Dataset type") +
  theme_bw(base_size = 25)+
  scale_y_continuous(breaks = c(0,1e06,2e06,3e06), labels = c("0","1,000,000","2,000,000","3,000,000"))+
  theme(axis.text = element_text(color = "black"),
        axis.label = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  
  labs(y = "Number of variants",
       x = NULL)

clinvar_all_class <- read_delim(clinvar_all_class, delim = "\t") %>% 
  distinct() %>% 
  filter(Gene %in% genes$Gene) %>% 
  filter(AA_ref != AA_alt)

hgmd_variants <- read_delim("Data_analysis/Genome_biology_2025/Datasets/Variants/HGMD_24_2_final.txt", delim = "\t") %>% 
  distinct() %>% 
  filter(AA_ref != AA_alt) # Ensure only missense variants

plot_B <- hgmd_variants %>% 
  left_join(clinvar_all_class) %>% 
  group_by(Classification) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(Classification = case_when(is.na(Classification) ~ "Not present\nin ClinVar",
                                    Classification == "Lb/B"~"Lb/B",
                                    Classification == "Lp/P"~"Lp/P",
                                    Classification == "VUS"~"VUS")) %>% 
  arrange(desc(n)) %>% 
  mutate(Classification = factor(Classification, levels = Classification)) %>% 
  ggplot(aes(x = Classification, y = n))+
  geom_bar(stat = "identity", position = position_dodge2(0.1))+
  geom_text(aes(y = n + 2000, label = n), size = 5, position = position_dodge2(0.9))+
  theme_bw(base_size = 25)+
  theme(axis.text = element_text(color = "black"),
        axis.label = element_text(color = "black"))+
  labs(y = "Number of variants",
       x = NULL)


#generate combined figure 
Figure4 <- (plot_A + plot_B) + 
  plot_layout(ncol = 2) + 
  plot_annotation(tag_levels = 'A')


# Save figure 
ggsave(figure_path, Figure4,  width = 16, height = 9, units = "in", dpi = 600)
