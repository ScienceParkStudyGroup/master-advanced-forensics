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



sample_age <- read.delim("data/cpg_methylation_sample_age.tsv", stringsAsFactors = FALSE)
df_cpg_tidy_with_age =  dplyr::left_join(df_cpg_tidy, 
                                         sample_age, 
                                         by  = "sample") 

# remove NAs and Inf values
df_for_rf <- na.omit(df_cpg_tidy_with_age)
df_for_rf <- df_for_rf[which(is.finite(df_for_rf$beta_coef)),]

rf_model <- randomForest(x = na.omit(df_for_rf), 
                         formula = age ~ .)
