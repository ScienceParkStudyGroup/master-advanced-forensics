suppressPackageStartupMessages(library("tidyverse"))



#############################
# 07. Linear regressions
#############################

### data import
df_cpg <- read.delim("data/cpg_methylation_beta_values.tsv",header = TRUE, stringsAsFactors = FALSE)
sample_age <- read.delim("data/cpg_methylation_sample_age.tsv", stringsAsFactors = FALSE)


df_cpg_tidy_with_age
