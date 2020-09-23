suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("data.table"))
source("code/pca_function.R") # custom function to perform Principal Component Analysis


#############
# Import data
#############



#####
# EDA
#####

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
  pivot_longer(cols = - probe, 
               names_to = "tissue", 
               values_to = "expression") 

# calculating simple metrics
df_tidy %>% 
  summarise(max = max(expression), 
            min = min(expression), 
            average = mean(expression), 
            median = median(expression)
            )
  
# Plotting the distribution of gene expression values per tissue
df_tidy %>% 
  #filter(tissue == "Adipocyte"|tissue == "Cerebellum") %>% 
  ggplot(data = ., aes(x = log2(expression), fill = tissue)) + 
    geom_density(alpha = 0.1) + 
    #facet_wrap(~ tissue) +
    theme(axis.text.x = element_text(angle = 90)) +
    guides(fill=FALSE)


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
  