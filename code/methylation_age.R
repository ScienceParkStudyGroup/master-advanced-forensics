suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("broom"))



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

#nested_df$aov_res[nested_df$CpG_id == "cg00000292"]

### Heatmap using only diff sites
df_for_heatmap <- pivot_wider(df_cpg_tidy, 
                              id_cols = CpG_id, 
                              names_from = "sample", 
                              values_from = "beta_coef") %>% 
  column_to_rownames("CpG_id")


# multiple linear regression
