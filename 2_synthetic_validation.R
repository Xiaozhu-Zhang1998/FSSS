rm(list = ls())
source("0_functions.R")
library(tidyverse)


# data generation ========
n_train = 200 * 3
n_test = 500
n = n_train + n_test
nclust = 3
nblock = 3
sindep = 5
p = 200
eta = 0.1
Xclust = gen_clust(n, nclust, d = 3, eta = 0.5)
Xblock1 = gen_block(n, 1, b = 2, eta = 0.01)
Xblock2 = gen_block(n, 1, b = 3, eta = 0.1)
Xblock3 = gen_block(n, 1, b = 4, eta = 0.1)
Xblock = cbind(Xblock1, Xblock2, Xblock3)
Xindep = matrix(rnorm(n * (p - ncol(Xclust) - ncol(Xblock) )), nrow = n)
X = cbind(Xclust, Xblock, Xindep)
hclustobj = hclust(as.dist(1-abs(cor(X))), method = "single")
cluster = as.vector(cutree(hclustobj, h = 0.3))

beta = c(
  rep(c(1,0,0), nclust),
  c(1,-1, 0),
  c(1, -1, 1, 0),
  c(1, -1, 1, -1, 0),
  rep(0.2, sindep),
  rep(0, p - ncol(Xclust) - ncol(Xblock) - sindep)
)
y = X %*% beta + rnorm(n, sd = 1.5)

synData = list(
  X = X, y = y, beta = beta
)

# save this dataset
saveRDS(synData, "synData.RDS")


# data processing ========
X = synData$X
y = synData$y

X = X - rep(1, nrow(X)) %*% t(colMeans(X))
X = X / ( rep(1, nrow(X)) %*% t(apply(X, 2, sd)) )
y = (y - mean(y)) / sd(y)

X = X[1:600,]
y = y[1:600]



# functions for validation ========
one_run_multi_alpha = function(X, y, s0, Alpha) {
  # divide the data
  n = length(y)
  n0 = round(n / 10 * 9)
  feature_id = sample(1:n, n0)
  test_id = setdiff(1:n, feature_id)
  cat("finished deviding data", "\n")
  
  # select
  feature_X = X[feature_id,]
  feature_y = y[feature_id]
  bags = l0_subsampling(feature_X, feature_y, s0, num_bags = 100)
  cat("finished subsamping", "\n")
  
  
  Mse_alpha = matrix(0, nrow = length(Alpha), ncol = 2)
  
  
  for(alpha in Alpha) {
    fit_FSSS = greedy_FSSS(feature_X, bags$Sinfo, bags$base_lst, alpha = alpha)
    model_FSSS = lin_reg(feature_X[, fit_FSSS$NODES, drop = FALSE], feature_y)
    
    test_X = X[test_id,]
    test_y = y[test_id]
    mse_FSSS = MSE(test_X[, fit_FSSS$NODES, drop = FALSE], model_FSSS, test_y)
    
    Mse_alpha[which(alpha == Alpha),] = c(alpha, mse_FSSS)
    cat("finished", "alpha", alpha, "\n\n")
  }
  
  return(Mse_alpha)
}



# execute ========
# repeat this section for s0 in S0 = seq(from = 2, to = 41, by = 1) 
# repeat each s0 100 times

# example of s0 = 30 and l = 1 (l is the index for repetition)
s0 = 30
l = 1

one_run_obj = one_run_multi_alpha(X, y, s0, Alpha = seq(from = 0.05, to = 0.2, by = 0.05))

saveRDS(one_run_obj, paste0("./validation/", s0, "_", l, ".RDS"))
cat("finished ...", "\n")



# collect results in the folder "validation" ========
validation = function(X, y, s0, Alpha, Lcv = 50) {
  CV_MSE = matrix(0, nrow = 0, ncol = 2)
  # collect results
  count = 1
  for(lcv in 1:Lcv) {
    one_run_obj = readRDS(paste0("./validation/", s0, "_", lcv, ".RDS"))
    CV_MSE = rbind(CV_MSE, one_run_obj)
  }
  # find the alpha
  CV_MSE = CV_MSE %>% data.frame() %>% `colnames<-`(c("alpha", "mse_FSSS")) %>% group_by(alpha) %>% summarise(mse = mean(mse_FSSS))
  return(unlist(CV_MSE[which.min(CV_MSE$mse) , ]))
}


CV_MSE = matrix(0, nrow = 0, ncol = 3)
for(s0 in S0) {
  CV_MSE = rbind(CV_MSE, c(s0, validation(X, y, s0, Alpha)))
  cat("finished", s0, "\n")
}

saveRDS(CV_MSE, "validation.RDS")

