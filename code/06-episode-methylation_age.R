suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("broom"))
source("code/save_pheatmap_function.R")



#########
# 06. EDA
#########

### data import
df_cpg <- read.delim("data/cpg_methylation_beta_values.tsv",header = TRUE, stringsAsFactors = FALSE)

df_cpg_tidy = df_cpg %>% 
  pivot_longer(cols = - "CpG_id", 
               names_to = "sample", 
               values_to = "beta_coef")

### density of beta values
df_cpg_tidy %>% 
  ggplot(data = ., aes(x = beta_coef, fill = sample)) + 
    geom_density(alpha = 0.1) +
  guides(fill = FALSE) 

ggsave(filename = "img/06-density-cpg.png")

### ANOVA beta_value ~ age
sample_age <- read.delim("data/cpg_methylation_sample_age.tsv", stringsAsFactors = FALSE)
df_cpg_tidy_with_age =  dplyr::left_join(df_cpg_tidy, 
                                         sample_age, 
                                         by  = "sample") 

# To do = List Column Workflow
# One-factor ANOVA on each CpG site
fits <- df_cpg_tidy_with_age %>% 
  nest(.data = ., - CpG_id) %>% # categorical variable used for ANOVA
  mutate(fit = map(data, ~ aov(data = ., formula = .$age ~ .$beta_coef))) %>% 
  mutate(res = map(fit, glance)) %>% 
  unnest(res)

fits_with_fdr <- fits %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

cpg_differentials <- filter(fits_with_fdr, fdr < 0.05) %>% 
  dplyr::pull(CpG_id)

length(cpg_differentials)

df_cpg_tidy_with_age_filtered <- filter(df_cpg_tidy_with_age, 
                                        CpG_id %in% cpg_differentials) 


### Heatmap using only diff sites
df_for_heatmap <- pivot_wider(df_cpg_tidy_with_age_filtered, 
                              id_cols = CpG_id, 
                              names_from = "sample", 
                              values_from = "beta_coef") %>% 
  column_to_rownames("CpG_id")

## Version 1: no sample annotation
my_heatmap <- pheatmap(df_for_heatmap, 
                       show_colnames = TRUE, 
                       show_rownames = FALSE, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE, 
                       fontsize_col = 4)

save_pheatmap_png(my_heatmap, "img/06-heatmap-1.png")

## Version 2: with sample age annotation
# adding sample age (color scale)
sample_annot <- sample_age %>% column_to_rownames("sample") %>% select(age)

my_heatmap2 <- pheatmap(df_for_heatmap, 
                       annotation_col = sample_annot,
                       show_colnames = TRUE, 
                       show_rownames = FALSE, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE, 
                       fontsize_col = 4)

save_pheatmap_png(my_heatmap2, "img/06-heatmap-2.png")

####
### Save table with only differential CpG sites
##

# filter out unnecessary columns
# comply with wide format (same as original one but with only 228 CpG sites)
final_df <- 
  df_cpg_tidy_with_age_filtered %>% 
  dplyr::select(CpG_id, sample, beta_coef) %>% 
  pivot_wider(id_cols = sample, names_from = CpG_id, values_from = beta_coef) 


# make it wide to comply with the original data format
write.table(x = final_df, 
            file = "data/differential_cpgs.tsv",
            quote = FALSE, 
            sep = "\t")
