suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggrepel"))
source("code/pca_function.R") # custom function to perform Principal Component Analysis


#############
# Import data
#############
df_expr <- read.delim(file = "data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.tsv", 
                 header = TRUE, 
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

df_expr_ten_tissues = dplyr::select(df_expr, - Description) %>% 
  dplyr::select(1:11)

#####
# PCA
#####

# It's usually a convention to have variables in columns and samples in rows.
# Here: Genes = variables 
# Here: Tissues = samples
# PCA: scores will be computed for samples/tissues and will give us their coordinates
# while loadings will be related to variables (genes)
# We need to transpose our dataset so that rows/columns are inverted
df_expr_transposed <- df_expr_ten_tissues %>% 
   column_to_rownames("gene_id") %>% 
  t()


### PCA object
# center: value of gene Xi,j - mean of that gene (Xi) where i = row index, j = column index
# scale: divide by gene standard deviation sd Xi
pca <- mypca(x = df_expr_transposed, center = TRUE, scale = TRUE)


### Screeplot 
# Create a dataframe with all PC components (n = number of tissues)
exp_var_df <- data.frame(PC = seq(1:nrow(pca$explained_var)), exp_var = pca$explained_var)

# make the complete screeplot
ggplot(exp_var_df, aes(x = PC, y = exp_var, label = PC)) +
  ylab('explained variance (%)') + 
  ggtitle('explained variance per component (all principal components') + 
  geom_bar(stat = "identity") +
  labs(x = "Principal Component number") + 
  scale_y_continuous(limits = c(0, 50)) +
  scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))

ggsave(filename = "img/03-pca-all-components.png")

### How many PCs?
exp_var_df = mutate(exp_var_df, cum_var = cumsum(exp_var)) 

ggplot(exp_var_df, aes(x = PC, y = cum_var)) +
  geom_point() + 
  geom_line(group = 1) + 
  labs(x = "Principal Component", y = "Cumulative Explained Variance (%)") +
  scale_x_continuous(breaks = 1:10) +
  geom_hline(yintercept = 80, color = "red")
  
ggsave(filename = "img/03-pca-cumulative-variance.png")

### Scores
scores <- pca$scores
scores[1:5,1:5]

# plot the scores of the first 2 components
explained_var = pca$explained_var$exp_var
tissue_names = row.names(scores)

ggplot(scores, aes(x = PC1, y = PC2, label = tissue_names)) + 
  geom_point(aes(colour = factor(tissue_names))) + 
  xlab(paste0('PC1 (',explained_var[1],'%)')) + 
  ylab(paste0('PC2 (',explained_var[2],'%)')) + 
  ggtitle('PCA score plot: PC1 and PC2') + 
  geom_text_repel() +
  guides(colour = FALSE)

ggsave(filename = "img/03-dataset-1-score-plot-with-names.png")


## PC2 and PC3
ggplot(scores, aes(x = PC3, y = PC2, label = tissue_names)) + 
  geom_point(aes(colour = factor(tissue_names))) + 
  xlab(paste0('PC3 (',explained_var[3],'%)')) + 
  ylab(paste0('PC2 (',explained_var[2],'%)')) + 
  ggtitle('PCA score plot: PC2 and PC3') + 
  geom_text_repel() +
  guides(colour = FALSE)

ggsave(filename = "img/03-dataset-1-score-plot-with-names-pc2-pc3.png")

## PC1 and PC3
ggplot(scores, aes(x = PC1, y = PC3, label = tissue_names)) + 
  geom_point(aes(colour = factor(tissue_names))) + 
  xlab(paste0('PC1 (',explained_var[1],'%)')) + 
  ylab(paste0('PC3 (',explained_var[3],'%)')) + 
  ggtitle('PCA score plot') + 
  geom_text_repel() +
  guides(colour = FALSE)

############
### Loadings
############

loadings <- pca$loadings %>%
  rownames_to_column("gene_id") 


### extract top genes  
top10genes_PC1_PC3 <- 
  loadings %>% 
  pivot_longer(cols = - gene_id, names_to = "PC", values_to = "loadings") %>% 
  dplyr::filter(PC == "PC1" | PC == "PC3") %>% 
  group_by(PC) %>% 
  dplyr::arrange(desc(abs(loadings))) %>% 
  dplyr::slice(1:10) %>% 
  left_join(x = ., y = df_expr[c("gene_id", "Description")], by = "gene_id")

top10genes_PC1_PC3


### loading plot

loadings4plot <- inner_join(loadings, top10genes_PC1_PC3, by = "gene_id") %>% 
  dplyr::select(gene_id, PC1, PC3)

ggplot(loadings4plot) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC3), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text_repel(data = loadings4plot, 
                  aes(x = PC1, y = PC3, label = gene_id),
                  size = 2) + 
  labs(x = "PC1", y = "PC3")

ggsave(filename = "img/03-pca-loading-plot-20-genes.png")




