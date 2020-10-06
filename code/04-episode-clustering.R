suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("broom"))


#############
# Import data
#############

df_expr <- read.delim(file = "data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.tsv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

mat_expr <- df_expr %>% 
  dplyr::select(- Description) %>% 
  column_to_rownames("gene_id") %>% 
  as.matrix()

###########
### Scaling
###########
min(mat_expr)
max(mat_expr)

mat_expr_scaled <- mat_expr %>% 
  t() %>% 
  scale(center = TRUE, scale = TRUE) %>% 
  t() %>% 
  na.omit()

mat_expr_scaled[1:5,1:5]

mean(mat_expr_scaled[1,]) # mean for ENSG00000223972.4  (close to 0)
sd(mat_expr_scaled[1,])   # standard deviation for ENSG00000223972.4 (unit variance of 1)

#################################
# Distance and method of linkage
###############################

# Distance measure: Euclidean = distance between two points 

distances_between_tissues <- dist(t(mat_expr_scaled), method = "euclidean")
as.matrix(distances_between_tissues)[1:5,1:5]

# ward.D2 is equivalent to AGNES method
hcl_ward <- cluster::agnes(x = distances_between_tissues,method = "ward")

par(mar=c(1,1,1,1))
png("img/04-dendogram-ward.png", width = 1000, height = 600)
plot(hcl_ward, main = "Tissue hierarchical clustering (AGNES and Ward's method)")
dev.off()

# The agglomerative coefficient  
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(d = distances_between_tissues, x) {
  agnes(d, method = x)$ac
}

# get agglomerative coefficient for each linkage method
purrr::map_dbl(m, function(x) ac(distances_between_tissues, x))
##   average    single  complete      ward 
## 0.9139303 0.8712890 0.9267750 0.9766577


##############
# HCL on genes
##############


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




