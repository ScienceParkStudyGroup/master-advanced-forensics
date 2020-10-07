---
title: "Finding candidate genes through feature engineering"
teaching: 30
exercises: 45
questions:
- "How do I transform my variables to build more meaningful variables?"
- "What type of transformation and checks should I perform?"
- "How to extract my list of selected tissue-specific genes?"
objectives:
- "Understand how distances are calculated between two tissues based on their gene expression profile."
- "Be able to name a few different clustering methods."
keypoints:
- ""
---

# Table of Contents

<!-- MarkdownTOC autolink="True" -->

- [1. Introduction](#1-introduction)
- [Feature engineering](#feature-engineering)
  - [](#)
- [References](#references)

<!-- /MarkdownTOC -->

# 1. Introduction

# Feature engineering

There is an alternative perhaps more intuitive method to find tissue-specific genes. 

## 

~~~

### Genes differential at least in one comparison
df_expr_tidy <- df_expr %>%
  select(- Description) %>% 
  pivot_longer(- gene_id, names_to = "tissue", values_to = "tpm")

# a lot of genes have very small TPM values
df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  with(., summary(median_tpm))

# 3rd quartile (genes with median TPM > 75th percentile)
# 14,049 genes instead of over 56,000
genes_selected = 
  df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  ungroup() %>% 
  filter(median_tpm > 1.7) %>% 
  dplyr::pull(gene_id)

df_expr_tidy_filtered <- filter(df_expr_tidy, gene_id %in% genes_selected)

## Step 1: extract gene TPM value in "Adipose - Subcutaneous"  
adipose_gene_expression <- df_expr_tidy_filtered %>% 
  filter(tissue == "Adipose - Subcutaneous") %>% 
  select(gene_id, tpm) %>% 
  rename(adipose_tpm = tpm)

## Step 2: calculate median TPM value in all other tissues
all_other_tissues_gene_expression <- 
  df_expr_tidy_filtered %>% 
  filter(tissue != "Adipose - Subcutaneous") %>% 
  group_by(gene_id) %>% 
  summarise(other_tissues_median_tpm = median(tpm))

## Merge the two dataframes
## Calculate a fold change
adipose_vs_other_tissues <- inner_join(x = adipose_gene_expression, 
                                       y = all_other_tissues_gene_expression, 
                                       by = "gene_id") %>% 
  mutate(fc = adipose_tpm / other_tissues_median_tpm) %>% 
  mutate(log2_fc = log2(fc)) 

# calculate Z-score 
mean_of_log2fc <- with(data = adipose_vs_other_tissues, mean(log2_fc))
sd_of_log2fc <- with(data = adipose_vs_other_tissues, sd(log2_fc))

adipose_vs_other_tissues$zscore <- map_dbl(
  adipose_vs_other_tissues$log2_fc, 
  function(x) (x - mean_of_log2fc) / sd_of_log2fc
  )

# normal distribution of log2FC?
ggplot(adipose_vs_other_tissues, aes(x = log2_fc)) +
  geom_density()

# test for normality
ks.test(x = adipose_vs_other_tissues$log2_fc, y = "pnorm", mean = mean_of_log2fc, sd = sd_of_log2fc)

# filter fc > 0 + calculate one-sided p-value
adipose_specific_genes = 
  adipose_vs_other_tissues %>%  
  filter(log2_fc > 0) %>% # FC superior to 1
  mutate(pval = 1 - pnorm(zscore)) %>% 
  filter(pval < 0.01) %>% 
  arrange(desc(log2_fc))

adipose_specific_genes
~~~
{: .language-r}



# References
- [Exploring gene expression patterns using clustering](https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html)