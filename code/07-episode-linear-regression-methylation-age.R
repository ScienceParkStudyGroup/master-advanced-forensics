suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("glmnet"))



#############################
# 07. Linear regressions
#############################

df <- read.delim("data/differential_cpgs.tsv",
                 header = TRUE, 
                 na.strings = "NA", 
                 stringsAsFactors = FALSE)

sample_age <- read.delim("data/cpg_methylation_sample_age.tsv",header = TRUE, stringsAsFactors = FALSE) %>% 
  tibble() %>% 
  dplyr::select(sample, age)

df_for_regression = inner_join(df, sample_age, by  = "sample") %>% 
  column_to_rownames("sample")

#### Multiple Linear regression

model1 <- lm(formula = age ~ cg00059225, data = df_for_regression)              # does work
model2 <- lm(formula = age ~ cg00059225 + cg00083937, data = df_for_regression)

#66 samples so < 229
df_for_regression_subset <- df_for_regression[,c(1:50,229)] # 50 variables and the age Y

model50 <- lm(formula = age ~ ., data = df_for_regression_subset)
broom::tidy(model1)  
broom::tidy(model2)  
broom::tidy(model50)  

# More variables than samples: infinite number of best OLS fits
model_impossible <- lm(formula = age ~ ., data = df_for_regression)
broom::tidy(model_impossible)  

#############################
# 08. Regularised regressions
#############################

cpgs <- df_for_regression %>% 
  dplyr::select(- age) %>%       # will be placed in another R object 
  as.matrix() %>%                # for glmnet function to work properly, needs a matrix data structure 
  na.omit()                      # remove NA values (glmnet has no tolerance for missing values) 

age <- df_for_regression %>% 
  na.omit() %>% 
  dplyr::pull(var = age)

my_lasso <- glmnet::glmnet(x = cpgs, 
                           y = age, 
                           standardize = FALSE,  # since all our values are between 0 and 1. No scale issue.
                           alpha = 1,            # alpha = 1 performs a LASSO analysis. 0 performs a Ridge analysis.  
                           family = "gaussian")

my_lambda_lasso_cv <- glmnet::cv.glmnet(
  x = cpgs, 
  y = age, 
  alpha = 1,
  nfolds = 10, 
  lambda = NULL, 
  type.measure = "mse")
plot(lambda_lasso, main = "Lasso penalty choice\n\n")


