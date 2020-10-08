suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("dendextend"))

#############
# Import data
#############

df_expr <- read.delim(file = "data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.tsv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE,
                      check.names = FALSE)


df_expr_tidy <- df_expr %>%
  select(- Description) %>% 
  pivot_longer(- gene_id, names_to = "tissue", values_to = "tpm")

##############
# Feature engineering
##############


### 90th percentile
threshold = df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  with(., round(x = quantile(x = median_tpm, probs = 0.9)))

genes_selected = 
  df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  ungroup() %>% 
  filter(median_tpm > threshold) %>% 
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
  geom_density(fill = "lightblue")

ggsave(filename = "img/05-log2-fc-distribution.png")

# test for normality
ks.test(x = adipose_vs_other_tissues$log2_fc, y = "pnorm", mean = mean_of_log2fc, sd = sd_of_log2fc)

# filter fc > 0 + calculate one-sided p-value
adipose_specific_genes = 
  adipose_vs_other_tissues %>%  
  filter(log2_fc > 0) %>% # FC superior to 1
  mutate(pval = 1 - pnorm(zscore)) %>% 
  filter(pval < 0.01) %>% 
  arrange(desc(log2_fc))

head(adipose_specific_genes, n = 10)


###########
# Clustering
##########

genes2keep <- adipose_specific_genes$gene_id[1:20]

adipose_mat <- df_expr_tidy %>% 
  dplyr::filter(gene_id %in% genes2keep) %>% 
  pivot_wider(id_cols = gene_id, names_from = "tissue", values_from = "tpm") %>% 
  column_to_rownames("gene_id") %>% 
  as.matrix()

adipose_mat_scaled <- adipose_mat %>% 
  t() %>% 
  scale(center = TRUE, scale = TRUE) %>% 
  t() %>% 
  na.omit()

pheatmap::pheatmap(adipose_mat_scaled, 
                   clustering_distance_cols = "euclidean", 
                   clustering_distance_rows = "euclidean", 
                   clustering_method = "ward.D2", 
                   show_rownames = TRUE, 
                   fontsize_row = 6)
