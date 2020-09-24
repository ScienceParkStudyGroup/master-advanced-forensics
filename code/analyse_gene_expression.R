suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("GGally"))
source("code/pca_function.R") # custom function to perform Principal Component Analysis


#########
# 02. EDA
#########

df_expr <- read.delim(file = "data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.tsv", 
                 header = TRUE, 
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

# number of rows / columns
nrow(df_expr)
ncol(df_expr)

# checking how variables are encoded
glimpse(df_expr) # to verify that R has imported probes as character and values as doubles (number)

# making the df tidy
df_tidy <- df_expr %>% 
  pivot_longer(cols = - c(gene_id, Description), 
               names_to = "tissue", 
               values_to = "expression") 

# calculating simple metrics
df_tidy %>% 
  summarise(max = max(expression), 
            min = min(expression), 
            average = mean(expression), 
            median = median(expression)
            )

### Plotting the distribution of gene expression values for one tissue ###
df_tidy %>% 
  dplyr::filter(tissue == "Adipose - Subcutaneous") %>% # the dplyr:: notation makes sure the filter function comes from the dplyr package
  ggplot(data = ., aes(x = tissue, y = expression)) +
    geom_boxplot() 

# Using log10 scaling
df_tidy %>% 
  dplyr::filter(tissue == "Adipose - Subcutaneous") %>% # the dplyr:: notation makes sure the filter function comes from the dplyr package
  ggplot(data = ., aes(x = tissue, y = log10(expression + 1))) +
  geom_boxplot() 

# removing null expression values
df_tidy %>% 
  dplyr::filter(tissue == "Adipose - Subcutaneous") %>% # the dplyr:: notation makes sure the filter function comes from the dplyr package
  dplyr::filter(expression != 0) %>% 
  ggplot(data = ., aes(x = tissue, y = log10(expression + 1))) +
  geom_boxplot()  

### Plotting the distribution of gene expression values for all tissues ###
df_tidy %>%
  dplyr::filter(expression != 0) %>% 
  filter(grepl(pattern = "Brain*", x = tissue)) %>% 
  ggplot(data = ., aes(x = tissue, y = log10(expression + 1), fill = tissue)) + 
    geom_boxplot(alpha = 0.1) + 
    theme(axis.text.x = element_text(angle = 90)) +
    guides(fill=FALSE) 


df_tidy %>% 
  mutate(log_expression = log10(expression)) %>% 
  ggplot(data = ., aes(x = log_expression, fill = tissue)) + 
    geom_density(alpha = 0.1) + 
    theme(axis.text.x = element_text(angle = 90)) +
    guides(fill=FALSE) 

### scatterplot matrix

# Between similar tissues (artery)
df_tidy %>%
  dplyr::filter(expression != 0) %>% 
  filter(grepl(pattern = "Artery*", x = tissue)) %>% 
  pivot_wider(id_cols = c(gene_id, Description), names_from = tissue, values_from = expression) %>% 
  dplyr::select(- gene_id, - Description) %>%  # as ggpairs only accept numerical values
  ggpairs(title = "Pairwise plot matrix of artery tissues", upper = "blank") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))

ggsave(filename = "img/01-artery.png")

# between dissimilar tissues (e.g Pancreas, Prostate, Lung)
df_tidy %>%
  dplyr::filter(expression != 0) %>% 
  dplyr::filter(tissue %in% c("Prostate","Pancreas", "Lung")) %>% 
  pivot_wider(id_cols = c(gene_id, Description), names_from = tissue, values_from = expression) %>% 
  dplyr::select(- gene_id, - Description) %>%  # as ggpairs only accept numerical values
  ggpairs(title = "Pairwise plot matrix of contrasted tissues", upper = "blank")

#####
# PCA
#####

mat_expr <- df_expr %>% 
  column_to_rownames("probe")

pca <- mypca(mat_expr, center = TRUE, scale = TRUE)

pc1_pc2_df <- 
  
  
  # add a convenient column number for the bar plot to display
  dfev <- data.frame(PC = c(1,2,3,4), exp_var  = pca$explained_var)

# make the plot
scree_plot <- ggplot(dfev, aes(x = PC, y = exp_var)) +
  ylab('explained variance (%)') + 
  ggtitle('explained variance per component') + 
  geom_bar(stat = "identity")

# display it
scree_plot
  