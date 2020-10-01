suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tree"))
suppressPackageStartupMessages(library("rpart"))
suppressPackageStartupMessages(library("randomForest"))



set.seed(101)
train = sample(1:nrow(boston), 300)


rf.boston = randomForest(medv~., data = boston, subset = train)

oob.err = double(13)
test.err = double(13)
for(mtry in 1:13){
  fit = randomForest(medv~., data = boston, subset=train, mtry=mtry, ntree = 350)
  oob.err[mtry] = fit$mse[350]
  pred = predict(fit, boston[-train,])
  test.err[mtry] = with(boston[-train,], mean( (medv-pred)^2 ))
}


df_cpg <- read.delim("data/cpg_methylation_beta_values.tsv",header = TRUE, stringsAsFactors = FALSE)
df_cpg_tidy = df_cpg %>% 
  pivot_longer(cols = - "CpG_id", 
               names_to = "sample", 
               values_to = "beta_coef")




# remove NAs and Inf values
df <- df_cpg_tidy_with_age[which(is.finite(df_cpg_tidy_with_age$beta_coef)),] %>% 
  na.omit(.) %>% 
  pivot_wider(id_cols = sample, names_from = CpG_id, values_from = beta_coef)

# subset 1000 columns to create a reasonable dataset
set.seed(123)
cols2sample = sample(x = 1:ncol(df), size = 1000)
df_subset <- column_to_rownames(.data = df, "sample") %>% select(all_of(cols2sample))

# add age
sample_age <- read.delim("data/cpg_methylation_sample_age.tsv", stringsAsFactors = FALSE) %>% 
  dplyr::select(sample, age)
df_for_rf <- rownames_to_column(df_subset, var = "sample") %>% 
  left_join(x = ., y = sample_age, by = "sample") %>% 
  column_to_rownames("sample") %>% 
  na.omit(.)

randomForest(age ~ ., data = df_for_rf)
