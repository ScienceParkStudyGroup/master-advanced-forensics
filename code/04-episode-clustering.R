suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("dendextend"))
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
hcl_tissue_ward <- cluster::agnes(x = distances_between_tissues,method = "ward")

par(mar=c(1,1,1,1))
png("img/04-dendogram-ward.png", width = 1000, height = 600)
plot(hcl_tissue_ward, main = "Tissue hierarchical clustering (AGNES and Ward's method)")
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

#################
# Gene clustering
#################

df_expr_tidy <- df_expr %>%
  select(- Description) %>% 
  pivot_longer(- gene_id, names_to = "tissue", values_to = "tpm")

### Filter genes based on a population percentile
percentile = 0.99

threshold = 
df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  with(., round(x = quantile(x = median_tpm, probs = percentile))) %>% 
  as.integer()

genes_selected = 
  df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  ungroup() %>% 
  filter(median_tpm > threshold) %>% 
  dplyr::pull(gene_id)


## clustering
genes_mat <- subset(x = mat_expr_scaled, 
                                 subset = rownames(mat_expr_scaled) %in% genes_selected)

distances_between_genes = dist(x = genes_mat, method = "euclidean")

# The AGNES clustering method coupled to Ward's cluster dissimilarity estimation method
hcl_genes_ward <- cluster::agnes(x = distances_between_genes, method = "ward")

# You can already create a dendrogram from this hierarchical cluster object (zoom to get additional details)
par(mar=c(4,5,2,2))
plot(hcl_genes_ward,
     labels = FALSE,
     which.plots = 2,
     main = "Gene hierarchical clustering (AGNES, Ward method)")

par(mar=c(1,1,1,1))
png("img/04-dendogram-genes-ward.png", width = 1000, height = 600)
plot(hcl_genes_ward,
     labels = FALSE,
     main = "Gene hierarchical clustering (AGNES and Ward's method)")
dev.off()

#########
# Heatmap
#########


p <- pheatmap(genes_mat, 
         show_rownames = FALSE, 
         cluster_rows = as.hclust(hcl_genes_ward), 
         cluster_cols = as.hclust(hcl_tissue_ward), 
         scale = "none")
p


png("img/04-heatmap.png", width = 1000, height = 600)
p
dev.off()

### Gene profiles per cluster

genes <- labels(as.dendrogram(hcl_genes_ward))


gene_to_cluster_group <- cutree(tree = hcl_genes_ward, k = 10) %>% 
  enframe() %>% 
  rename(gene = name, cluster = value) %>% 
  mutate(gene_id = genes)


mat_expr_scaled_with_clusters = 
  mat_expr_scaled %>%
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>% 
  pivot_longer(cols = - gene_id, names_to = "tissue", values_to = "scaled_value") %>%  
  inner_join(x = ., y = gene_to_cluster_group, by = "gene_id") 
  

ggplot(mat_expr_scaled_with_clusters, aes(x = tissue, y = scaled_value)) +
  geom_line(aes(group = gene_id)) +
  facet_wrap( ~ cluster) +
  stat_summary(aes(group = 1), fun = "median", colour = "brown", geom = "line", size = 1.5) +
  theme(axis.text.x = element_text(angle = 90, size = 1))

ggsave(filename = "img/04-gene-expression-profiles.png")

