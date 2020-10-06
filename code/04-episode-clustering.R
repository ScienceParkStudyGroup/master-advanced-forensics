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

df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  with(., summary(median_tpm))

df_expr_tidy %>% 
  group_by(gene_id) %>% 
  summarise(median_tpm = median(tpm)) %>% 
  filter(median_tpm)

test <- df_expr_tidy %>% filter(gene_id == "ENSG00000223972.4")

fits <- df_expr_tidy %>% 
  nest(.data = ., - gene_id) %>%                                                  # categorical variable used for ANOVA
  mutate(fit = map(data, ~ aov(data = ., formula = .$tpm ~ .$tissue))) %>%        # One-factor ANOVA with age as the Y and beta value as X 
  mutate(res = map(fit, glance)) %>%                                             # Use broom::glance to return a tibble from a aov/lm object  
  unnest(res)



