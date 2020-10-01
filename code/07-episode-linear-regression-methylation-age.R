suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("glmnet"))



#############################
# 07. Linear regressions
#############################

df <- read.delim("data/differential_cpgs_with_sample_age.tsv",
                 header = TRUE, 
                 na.strings = "NA", 
                 stringsAsFactors = FALSE)

# Format desired

# | sample | cg00000292 | cg000002426 | ...| age |
# |--------||--------|  |--------|  |--------|    
# GSM712303|
# GSM712435|

df_wide <- 
  dplyr::select(df, - age) %>% 
  pivot_wider(id_cols = CpG_id, names_from = sample, values_from = beta_coef) 


sample_age <- read.delim("data/cpg_methylation_sample_age.tsv",header = TRUE, stringsAsFactors = FALSE) %>% 
  tibble() %>% 
  dplyr::select(sample, age)

df_for_regression = inner_join(df_wide, sample_age, by  = "sample") %>% 
  column_to_rownames("sample")

#### Multiple Linear regression

#66 samples so 66 vairables
test <- df_for_regression[,c(1:50,229)]

model1 <- lm(formula = age ~ ., data = test)
model2 <- lm(formula = age ~ cg00059225 + cg00083937, data = df_for_regression)


broom::tidy(model1)  
broom::tidy(model2)  

#############################
# 08. Regularised regressions
#############################

cpgs <- df_for_regression %>% 
  dplyr::select(- age)           # will be placed in another R object 
  as.matrix() %>%                # for glmnet function to work properly, needs a matrix data structure 
  na.omit()                      # remove NA values (glmnet has no tolerance for missing values) 

age <- df_for_regression %>% 
  na.omit() %>% 
  dplyr::pull(var = age)

glmnet::glmnet(x = cpgs, 
               y = age, 
               alpha = 1, 
               family = "gaussian")





