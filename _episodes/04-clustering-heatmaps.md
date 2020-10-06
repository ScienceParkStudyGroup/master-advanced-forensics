---
title: "Finding candidate genes through clustering and feature engineering"
teaching: 30
exercises: 45
questions:
- "How can I group genes with similar profiles together?"
- "How do I represent the result of a clustering analysis?"
objectives:
- "Understand how distances are calculated between two tissues based on their gene expression profile."
- "Be able to name a few different clustering methods."
keypoints:
- ""
---

# Table of Contents

<!-- MarkdownTOC autolink="True" -->

- [1. Introduction](#1-introduction)
    - [1.1 Objective](#11-objective)
    - [1.1 Library and Data import](#11-library-and-data-import)
- [2. Hierarchical Clustering](#2-hierarchical-clustering)
    - [2.1 Scaling](#21-scaling)
    - [2.2 Distance matrix](#22-distance-matrix)
        - [Simple example](#simple-example)
        - [Real example](#real-example)
    - [2.3 Tissue hierarchical clustering](#23-tissue-hierarchical-clustering)
    - [2.4 Gene clustering](#24-gene-clustering)
    - [2.4 Extract cluster to gene correspondence](#24-extract-cluster-to-gene-correspondence)
- [Feature engineering](#feature-engineering)
    - [](#)
- [References](#references)

<!-- /MarkdownTOC -->



# 1. Introduction

## 1.1 Objective 
Since PCA did not exactly yielded obvious adipose-specific gene candidates, here we are going to implement a _clustering_ approach. This method will try to identify groups of genes with similar expression profiles in subcutaneous adipose tissues compared to other tissues.      

In this section, clustering will refer exclusively to __hierarchical clustering__, a process to identify groups in a dataset. Contrarily to to other methods such as K-means clustering for instance, one doesn't know the number of clusters in advance.   

> ## Objective
> Our objective in this episode here is to define a cluster of genes that behave similarly in subcutaneous adipose tissue as compared to the other tissues. 
{: .objectives}

## 1.1 Library and Data import 

We import the libraries and dataset once more to provide a clean start.  
~~~
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("cluster"))


df_expr <- read.delim(file = "data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.tsv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE,
                      check.names = FALSE)

# conversion to matrix for distance calculation later
mat_expr <- df_expr %>% 
  dplyr::select(- Description) %>% 
  column_to_rownames("gene_id") %>% 
  as.matrix()
~~~
{: .language-r}

# 2. Hierarchical Clustering

## 2.1 Scaling 

As seen previously, scaling is a procedure to ensure that all genes (variables) are comparable before doing further analysis. For instance, some genes will have a very high expression while others will be close to zero.   
Try this:  
~~~
min(mat_expr)
max(mat_expr)
~~~
{: .language-r}

We perform scaling (mean centering and variance scaling0 on the expression matrix using the `scale()` function. Since this function scales __columns__ rather than rows, we will have to transpose the matrix before scaling and once more to return it to its original format.
~~~
mat_expr_scaled <- mat_expr %>% 
  t() %>% 
  scale(center = TRUE, scale = TRUE) %>% 
  t() %>% 
  na.omit()

mean(mat_expr_scaled[1,]) # mean for ENSG00000223972.4  (close to 0)
sd(mat_expr_scaled[1,])   # standard deviation for ENSG00000223972.4 (unit variance of 1)
~~~
{: .language-r}


## 2.2 Distance matrix

To perform the clustering itself, we need to compute a so-called distance matrix between every tissue in a pairwise manner. 
Here, we will use the __Euclidean distance__ and compute it for each tissue pair. 


### Simple example
Let's explain how it is calculated using a simple example. Feel free to run the following code:  

~~~
# sample 1 and 2
sample_a = c(1,2,3)
sample_b = c(1,2,3)

# Source: https://stackoverflow.com/questions/5559384/euclidean-distance-of-two-vectors
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
euc.dist(sample_a, sample_b)
~~~
{: .language-r}

The distance is 0 because the two samples have the same coordinates in their 3-dimensional space (x = 1, y = 1, z = 1).

~~~
sample_a = c(1,2,6)
sample_b = c(1,2,3)
euc.dist(sample_a, sample_b) # equal to 3

# if you do it manually
sqrt(sum(1 - 1)^2 + sum(2 - 2)^2 + sum(6 -3)^2)
~~~
{: .language-r}

### Real example

In our case, for every tissue, we have 56,202 gene expression which would be impossible to calculate by hand. Fortunately, R will do this for us. 

We will use the `dist()` function for this. In the related `dist()` function help, it says: 
> This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.

Since we do not want the distance between genes but rather between tissues, we have to transpose the matrix first. 
~~~
distances_between_tissues <- dist(
    t(mat_expr_scaled),                       # notice the t() to calculate the distances between tissues, not between genes
    method = "euclidean")
as.matrix(distances_between_tissues)[1:5,1:3]
~~~
{: .language-r}

You can notice that the distances between the same tissue is equal to 0 (null distance).
~~~
                             Adipose - Subcutaneous Adipose - Visceral (Omentum) Adrenal Gland
Adipose - Subcutaneous                      0.00000                     96.05766      200.2761
Adipose - Visceral (Omentum)               96.05766                      0.00000      182.7446
Adrenal Gland                             200.27614                    182.74462        0.0000
~~~
{: .language-r}


## 2.3 Tissue hierarchical clustering 

Let's use the AGNES (AGlomerative NESting) method from the `cluster` package as it produces smaller clusters than other methods. This could help to identify small meaningful tissue-specific clusters. In addition, we will use the Ward's method to estimate the distance between two clusters of tissues. Ward's method minimizes the within-cluster variance. 

~~~
# The AGNES clustering method coupled to Ward's cluster dissimilarity estimation method
hcl_ward <- cluster::agnes(x = distances_between_tissues,method = "ward")

You can already create a dendrogram from this hierarchical cluster object (zoom to get additional details)
plot(hcl_ward)
~~~
{: .language-r}

<img src="../img/04-dendogram-ward.png">

Notice how the two _adipose_ tissues are clustered together on the left with the _breast_ tissue. 

  
> ## Note 
> Notice that an agglomerative coefficient (ac) equal to 0.84 is indicated below the plot. This coefficient measures the clustering structure of the dataset.
> For each observation i, denote by m(i) its dissimilarity to the first cluster it is merged with, divided by the dissimilarity of the merger in the final step
> of the algorithm. The ac is the average of all 1 - m(i). It can also be seen as the average width (or the percentage filled) of the banner plot. Because ac > grows with the number of observations, this measure should not be used to compare datasets of very different sizes. 
{: .callout}

Let's see how different methods of cluster dissimilarity performs relative to this agglomerative coefficient.
~~~
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(d = distances_between_tissues, x) {
  agnes(d, method = x)$ac
}

# get agglomerative coefficient for each linkage method
purrr::map_dbl(m, function(x) ac(distances_between_tissues, x))
~~~
{: .language-r}

According to the results, it seems that Ward's method gives the highest agglomerative coefficient.
~~~
##   average    single  complete      ward 
## 0.9139303 0.8712890 0.9267750 0.9766577
~~~
{: .output}

This analysis suggests that some gene expression profile can indeed drive 

## 2.4 Gene clustering
<center>In construction</center> :construction_worker: <center>In construction</center>

## 2.4 Extract cluster to gene correspondence
<center>In construction</center> :construction_worker: <center>In construction</center>

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
