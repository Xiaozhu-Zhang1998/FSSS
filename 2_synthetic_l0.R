rm(list = ls())
source("0_functions.R")
library(tidyverse)

# repeat the following for s0 in S0 = seq(from = 2, to = 41, by = 1) 
# repeat each s0 500 times

# example of s0 = 30 and l__ = 1 (l__ is the index for repetition)
s0 = 30
l__ = 1



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


X = synData$X
y = synData$y

X = X - rep(1, nrow(X)) %*% t(colMeans(X))
X = X / ( rep(1, nrow(X)) %*% t(apply(X, 2, sd)) )
y = (y - mean(y)) / sd(y)


id_test = sample(seq_along(y), 500)
id_train = setdiff(seq_along(y), id_test)

data = list(
  X_train = X[id_train,], 
  y_train = y[id_train], 
  X_test = X[id_test,], 
  y_test = y[id_test]
)



# functions for simulation ========
one_run_multi_alpha = function(data, s0, cutoff, Alpha) {
  X_train = as.matrix(data$X_train)
  y_train = data$y_train
  X_test = as.matrix(data$X_test)
  y_test = data$y_test
  
  
  # divide the data
  n = length(y_train)
  n0 = round(n / 3)
  feature_id = sample(1:n, n0)
  model_id = sample(setdiff(1:n, feature_id), n0)
  test_id = setdiff(1:n, union(feature_id, model_id))
  cat("finished deviding data", "\n")
  
  
  # cluster
  nontest_X = X_train[-test_id,]
  hclustobj = hclust(as.dist(1-abs(cor(nontest_X))), method = "single")
  cluster = as.vector(cutree(hclustobj, h = cutoff)) 
  cat("computed cluster", "\n")
  
  
  # select
  feature_X = X_train[feature_id,]
  feature_y = y_train[feature_id]
  bags = l0_subsampling(feature_X, feature_y, s0, num_bags = 100)
  cat("finished subsamping", "\n")
  
  
  Mse_alpha = matrix(0, nrow = length(Alpha), ncol = 5)
  
  for(alpha in Alpha) {
    
    if(length(alpha) == 1) {
      alpha_SS = alpha_CSS = alpha_FSSS = alpha
    } else {
      alpha_SS = alpha[1]
      alpha_CSS = alpha[2]
      alpha_FSSS = alpha[3]
    }
    
    fit_SS = which(selprob(bags) > 1 - alpha_SS)
    fit_CSS_wavg = cluSS_pred_alpha(bags$Selset, cluster, alpha = alpha_CSS, "wavg")
    fit_CSS_savg = cluSS_pred_alpha(bags$Selset, cluster, alpha = alpha_CSS, "savg")
    fit_CSS_sps = cluSS_pred_alpha(bags$Selset, cluster, alpha = alpha_CSS, "sps")
    fit_FSSS = greedy_FSSS(feature_X, bags$Sinfo, bags$base_lst, alpha = alpha_FSSS)
    cat("finished selection", "\n")
    
    # model
    model_X = X_train[model_id,]
    model_y = y_train[model_id]
    model_SS = lin_reg(model_X[, fit_SS, drop = FALSE], model_y)
    X_CSS_wavg = create_X_cluster(fit_CSS_wavg, model_X)
    model_CSS_wavg = lin_reg(X_CSS_wavg, model_y)
    X_CSS_savg = create_X_cluster(fit_CSS_savg, model_X)
    model_CSS_savg = lin_reg(X_CSS_savg, model_y)
    X_CSS_sps = create_X_cluster(fit_CSS_sps, model_X)
    model_CSS_sps = lin_reg(X_CSS_sps, model_y)
    model_FSSS = lin_reg(model_X[, fit_FSSS$NODES, drop = FALSE], model_y)
    cat("finished OLS", "\n")
    
    # test
    test_X = X_train[test_id,]
    test_y = y_train[test_id]
    mse_SS = MSE(test_X[, fit_SS, drop = FALSE], model_SS, test_y)
    X_CSS_wavg = create_X_cluster(fit_CSS_wavg, test_X)
    mse_CSS_wavg = MSE(X_CSS_wavg, model_CSS_wavg, test_y)
    X_CSS_savg = create_X_cluster(fit_CSS_savg, test_X)
    mse_CSS_savg = MSE(X_CSS_savg, model_CSS_savg, test_y)
    X_CSS_sps = create_X_cluster(fit_CSS_sps, test_X)
    mse_CSS_sps = MSE(X_CSS_sps, model_CSS_sps, test_y)
    mse_FSSS = MSE(test_X[, fit_FSSS$NODES, drop = FALSE], model_FSSS, test_y)
    cat("finished computing MSE", "\n")
    
    # record
    Mse_alpha[which(alpha == Alpha),] = c(cutoff, alpha, mse_SS, max(mse_CSS_wavg, mse_CSS_savg, mse_CSS_sps), mse_FSSS)
    cat("finished", "alpha", alpha, "\n\n")
  }
  
  return(Mse_alpha)
}




one_run = function(data, s0, cutoff, alpha) {
  X_train = as.matrix(data$X_train)
  y_train = data$y_train
  X_test = as.matrix(data$X_test)
  y_test = data$y_test
  
  # divide the data
  n = length(y_train)
  n0 = round(n / 3)
  feature_id = sample(1:n, n0)
  model_id = sample(setdiff(1:n, feature_id), n0)
  test_id = setdiff(1:n, union(feature_id, model_id))
  cat("finished deviding data", "\n")
  
  
  # cluster
  nontest_X = X_train[-test_id,]
  hclustobj = hclust(as.dist(1-abs(cor(nontest_X))), method = "single")
  # plot(hclustobj)
  cluster = as.vector(cutree(hclustobj, h = cutoff)) 
  cat("computed cluster", "\n")
  
  
  # select
  feature_X = X_train[feature_id,]
  feature_y = y_train[feature_id]
  bags = l0_subsampling(feature_X, feature_y, s0, num_bags = 100)
  cat("finished subsamping", "\n")
  
  
  if(length(alpha) == 1) {
    alpha_SS = alpha_CSS = alpha_FSSS = alpha
  } else {
    alpha_SS = alpha[1]
    alpha_CSS = alpha[2]
    alpha_FSSS = alpha[3]
  }
  
  
  fit_SS = which(selprob(bags) > 1 - alpha_SS)
  fit_CSS_wavg = cluSS_pred_alpha(bags$Selset, cluster, alpha = alpha_CSS, "wavg")
  fit_CSS_savg = cluSS_pred_alpha(bags$Selset, cluster, alpha = alpha_CSS, "savg")
  fit_CSS_sps = cluSS_pred_alpha(bags$Selset, cluster, alpha = alpha_CSS, "sps")
  fit_FSSS = greedy_FSSS(feature_X, bags$Sinfo, bags$base_lst, alpha = alpha_FSSS)
  cat("finished selection", "\n")
  
  
  # model
  model_X = X_train[model_id,]
  model_y = y_train[model_id]
  model_SS = lin_reg(model_X[, fit_SS, drop = FALSE], model_y)
  X_CSS_wavg = create_X_cluster(fit_CSS_wavg, model_X)
  model_CSS_wavg = lin_reg(X_CSS_wavg, model_y)
  X_CSS_savg = create_X_cluster(fit_CSS_savg, model_X)
  model_CSS_savg = lin_reg(X_CSS_savg, model_y)
  X_CSS_sps = create_X_cluster(fit_CSS_sps, model_X)
  model_CSS_sps = lin_reg(X_CSS_sps, model_y)
  model_FSSS = lin_reg(model_X[, fit_FSSS$NODES, drop = FALSE], model_y)
  cat("finished OLS", "\n")
  
  
  # lasso 
  fit_lasso = glmnet::glmnet(x = feature_X, y = feature_y)
  lasso_pw = which(sapply(1:ncol(fit_lasso$beta), function(i) {
    sum(fit_lasso$beta[,i] != 0)
  } ) <= s0)
  lambda_id = ifelse(length(lasso_pw) == 0, 1, max(lasso_pw))
  fit_Lasso = as.numeric(which(fit_lasso$beta[, lambda_id] != 0 ))
  model_Lasso = lin_reg(model_X[, fit_Lasso, drop = FALSE], model_y)
  
  
  # l0
  fit_l0 = L0Learn::L0Learn.fit(x = feature_X, y = feature_y, maxSuppSize = s0, intercept = TRUE)
  beta_l0 = coef(fit_l0)
  beta_l0 = beta_l0[2:nrow(beta_l0), ncol(beta_l0)]
  fit_L0 = as.numeric( which(beta_l0 != 0) )
  model_L0 = lin_reg(model_X[, fit_L0, drop = FALSE], model_y)
  
  
  return_list = list(
    fit_SS = fit_SS,
    fit_CSS_wavg = fit_CSS_wavg,
    fit_CSS_savg = fit_CSS_savg,
    fit_CSS_sps = fit_CSS_sps,
    fit_FSSS = fit_FSSS$NODES,
    fit_Lasso = fit_Lasso,
    fit_L0 = fit_L0,
    model_SS = model_SS,
    model_CSS_wavg = model_CSS_wavg,
    model_CSS_savg = model_CSS_savg,
    model_CSS_sps = model_CSS_sps,
    model_FSSS = model_FSSS,
    model_Lasso = model_Lasso,
    model_L0 = model_L0,
    feature_id = feature_id
  )
  return(return_list)
}




validation = function(data, s0, Cutoff = c(0.1, 0.3, 0.5), Alpha = seq(from = 0.05, to = 0.2, by = 0.05), Lcv = 10) {
  CV_MSE = matrix(0, nrow = 0, ncol = 5)
  
  # cross validation
  count = 1
  for(cutoff in Cutoff) {
    for(lcv in 1:Lcv) {
      one_run_obj = one_run_multi_alpha(data, s0, cutoff, Alpha)
      CV_MSE = rbind(CV_MSE, one_run_obj)
      cat("finished", "cutoff", cutoff," ", "lcv", lcv, "\n\n")
    }
  }
  
  # find the cutoff and corresponding estimate
  CV_MSE = CV_MSE %>% data.frame() %>% `colnames<-`(c("cutoff", "alpha", "mse_SS", "mse_CSS", "mse_FSSS")) 
  
  CV_MSE_SS = CV_MSE %>% select(cutoff, alpha, mse_SS) %>% group_by(cutoff, alpha) %>% summarise(mse = mean(mse_SS))
  selid_SS = CV_MSE_SS %>% pull(mse) %>% which.min()
  CV_MSE_CSS = CV_MSE %>% select(cutoff, alpha, mse_CSS) %>% group_by(cutoff, alpha) %>% summarise(mse = mean(mse_CSS))
  selid_CSS = CV_MSE_CSS %>% pull(mse) %>% which.min()
  CV_MSE_FSSS = CV_MSE %>% select(cutoff, alpha, mse_FSSS) %>% group_by(cutoff, alpha) %>% summarise(mse = mean(mse_FSSS))
  selid_FSSS = CV_MSE_FSSS %>% pull(mse) %>% which.min()
  
  alpha = c(CV_MSE_SS[[selid_SS, 2]], CV_MSE_SS[[selid_CSS, 2]], CV_MSE_FSSS[[selid_FSSS, 2]])
  cutoff = CV_MSE_CSS[[selid_CSS, 1]]
  
  return(list(
    alpha = alpha,
    cutoff = cutoff
  ))
}




one_experiment = function(data, s0, cutoff, alpha) {
  # run one experiment
  one_run_obj = one_run(data, s0, cutoff, alpha)
  
  fit_SS = one_run_obj$fit_SS
  fit_CSS_wavg = one_run_obj$fit_CSS_wavg
  fit_CSS_savg = one_run_obj$fit_CSS_savg
  fit_CSS_sps = one_run_obj$fit_CSS_sps
  fit_FSSS = one_run_obj$fit_FSSS
  fit_Lasso = one_run_obj$fit_Lasso
  fit_L0 = one_run_obj$fit_L0
  model_SS = one_run_obj$model_SS
  model_CSS_wavg = one_run_obj$model_CSS_wavg
  model_CSS_savg = one_run_obj$model_CSS_savg
  model_CSS_sps = one_run_obj$model_CSS_sps
  model_FSSS = one_run_obj$model_FSSS
  model_Lasso = one_run_obj$model_Lasso
  model_L0 = one_run_obj$model_L0
  
  # calculate MSE on [test/report]
  X_test = as.matrix(data$X_test)
  y_test = data$y_test
  # find newX
  X_CSS_wavg = create_X_cluster(fit_CSS_wavg, X_test)
  X_CSS_savg = create_X_cluster(fit_CSS_savg, X_test)
  X_CSS_sps = create_X_cluster(fit_CSS_sps, X_test)
  # calculate mse
  mse_SS = MSE(X_test[, fit_SS, drop = FALSE], model_SS, y_test)
  mse_CSS_wavg = MSE(X_CSS_wavg, model_CSS_wavg, y_test)
  mse_CSS_savg = MSE(X_CSS_savg, model_CSS_savg, y_test)
  mse_CSS_sps = MSE(X_CSS_sps, model_CSS_sps, y_test)
  mse_FSSS = MSE(X_test[, fit_FSSS, drop = FALSE], model_FSSS, y_test)
  mse_Lasso = MSE(X_test[, fit_Lasso, drop = FALSE], model_Lasso, y_test)
  mse_L0 = MSE(X_test[, fit_L0, drop = FALSE], model_L0, y_test)
  # nubmer of selected features
  q_SS = length(fit_SS)
  q_CSS = length(fit_CSS_wavg$Cluster)
  q_FSSS = length(fit_FSSS)
  q_Lasso = length(fit_Lasso)
  q_L0 = length(fit_L0)
  
  # interpretability
  inter_SS = ifelse(length(fit_SS) == 0, 0, 1)
  inter_FSSS = ifelse(length(fit_FSSS) == 0, 0, 1)
  inter_CSS_wavg = interpretability(fit_CSS_wavg)
  inter_CSS_savg = interpretability(fit_CSS_savg)
  inter_CSS_sps = ifelse(length(fit_CSS_sps$Cluster) == 0, 0, 1)
  inter_Lasso = ifelse(length(fit_Lasso) == 0, 0, 1)
  inter_L0 = ifelse(length(fit_L0) == 0, 0, 1)
  
  # find PS
  PStar = proj(X_test[, synData$beta != 0])
  PSperp = diag(rep(1, nrow(X_test))) - PStar
  # calculate FD
  FD_SS = ifelse(length(fit_SS)==0, 0, tr(PSperp %*% proj(X_test[, fit_SS, drop = FALSE])) )
  FD_CSS_wavg = ifelse(length(fit_CSS_wavg$Cluster) == 0, 0, tr(PSperp %*% proj(X_CSS_wavg)) ) 
  FD_CSS_savg = ifelse(length(fit_CSS_savg$Cluster) == 0, 0, tr(PSperp %*% proj(X_CSS_savg)) ) 
  FD_CSS_sps = ifelse(length(fit_CSS_sps$Cluster) == 0, 0, tr(PSperp %*% proj(X_CSS_sps)) ) 
  FD_FSSS = ifelse(length(fit_FSSS) == 0, 0, tr(PSperp %*% proj(X_test[, fit_FSSS, drop = FALSE])) ) 
  FD_Lasso = ifelse(length(fit_Lasso)==0, 0, tr(PSperp %*% proj(X_test[, fit_Lasso, drop = FALSE])) )
  FD_L0 = ifelse(length(fit_L0)==0, 0, tr(PSperp %*% proj(X_test[, fit_L0, drop = FALSE])) )
  # calculate PW
  PW_SS = ifelse(length(fit_SS)==0, 0, tr(PStar %*% proj(X_test[, fit_SS, drop = FALSE])) )
  PW_CSS_wavg = ifelse(length(fit_CSS_wavg$Cluster) == 0, 0, tr(PStar %*% proj(X_CSS_wavg)) )
  PW_CSS_savg = ifelse(length(fit_CSS_savg$Cluster) == 0, 0, tr(PStar %*% proj(X_CSS_savg)) )
  PW_CSS_sps = ifelse(length(fit_CSS_sps$Cluster) == 0, 0, tr(PStar %*% proj(X_CSS_sps)) )
  PW_FSSS = ifelse(length(fit_FSSS) == 0, 0, tr(PStar %*% proj(X_test[, fit_FSSS, drop = FALSE])) )
  PW_Lasso = ifelse(length(fit_Lasso) == 0, 0, tr(PStar %*% proj(X_test[, fit_Lasso, drop = FALSE])) )
  PW_L0 = ifelse(length(fit_L0) == 0, 0, tr(PStar %*% proj(X_test[, fit_L0, drop = FALSE])) )
  
  return(list(
    summary = c(
      alpha, cutoff, 
      mse_SS, mse_CSS_wavg, mse_CSS_savg, mse_CSS_sps, mse_FSSS, mse_Lasso, mse_L0,
      q_SS, q_CSS, q_FSSS, q_Lasso, q_L0,
      FD_SS, FD_CSS_wavg, FD_CSS_savg, FD_CSS_sps, FD_FSSS, FD_Lasso, FD_L0,
      PW_SS, PW_CSS_wavg, PW_CSS_savg, PW_CSS_sps, PW_FSSS, PW_Lasso, PW_L0,
      inter_SS, inter_CSS_wavg, inter_CSS_savg, inter_CSS_sps, inter_FSSS, inter_Lasso, inter_L0
    ),
    fit_SS = fit_SS,
    fit_CSS_wavg = fit_CSS_wavg,
    fit_CSS_savg = fit_CSS_savg,
    fit_CSS_sps = fit_CSS_sps,
    fit_FSSS = fit_FSSS,
    fit_Lasso = fit_Lasso,
    fit_L0 = fit_L0,
    feature_id = one_run_obj$feature_id
  ))
}




stability = function(X, fit1, fit2, type) {
  
  stab_metric = function(S1, S2, X) {
    if(length(S1) == 0 | length(S2) == 0) return(1)
    U = svd(X[, S1, drop = FALSE])$u
    V = svd(X[, S2, drop = FALSE])$u
    tr(U %*% t(U) %*% V %*% t(V)) / max(length(S1), length(S2))
  }
  
  if(type %in% c("SS", "Lasso", "L0", "FSSS")) {
    return( stab_metric(fit1, fit2, X) )
  }
  
  if(type == "CSS") {
    S1 = unlist( sapply(seq_along(fit1$Weights), function(i) { 
      fit1$Weights[[i]]$feature[ fit1$Weights[[i]]$weight != 0 ]
    } ))
    S2 = unlist( sapply(seq_along(fit2$Weights), function(i) { 
      fit2$Weights[[i]]$feature[ fit2$Weights[[i]]$weight != 0 ]
    } ))
    return( stab_metric(S1 = S1, S2 = S2, X = X) )
  }
  
}




interpretability = function(fit_CSS) {
  if(length(fit_CSS$Cluster) == 0) return(0)
  mean(1 / sapply(seq_along(fit_CSS$Weights), function(i) { 
    sum(fit_CSS$Weights[[i]]$weight != 0)
  } ))
}



# execute ========
validation_obj = validation(data, s0, Cutoff = c(0.1, 0.3, 0.5), Alpha = seq(from = 0.05, to = 0.3, by = 0.05) )
L = 50
RS = matrix(0, nrow = 0, ncol = 37)
STORE_OBJ = list()
for(rep in 1:L) {
  one_experiment_obj = one_experiment(data, s0, cutoff = validation_obj$cutoff, alpha = validation_obj$alpha)
  RS = rbind(RS, one_experiment_obj$summary)
  STORE_OBJ = c(STORE_OBJ, list(one_experiment_obj))
  cat('finished', rep, "out of", L, "\n\n")
}
rs_ = colMeans(RS, na.rm = TRUE)



# stability
stab_SS = 0
stab_CSS_wavg = 0
stab_CSS_savg = 0
stab_CSS_sps = 0
stab_FSSS = 0
stab_Lasso = 0
stab_L0 = 0
for(l in 1:L) {
  if (l == L)  next
  for(l_ in (l+1):L ) {
    rs1 = STORE_OBJ[[l]]
    rs2 = STORE_OBJ[[l_]]
    stab_SS = stab_SS + stability(X, rs1$fit_SS, rs2$fit_SS, "SS") 
    stab_CSS_wavg = stab_CSS_wavg + stability(X, rs1$fit_CSS_wavg, rs2$fit_CSS_wavg, "CSS") 
    stab_CSS_savg = stab_CSS_savg + stability(X, rs1$fit_CSS_savg, rs2$fit_CSS_savg, "CSS") 
    stab_CSS_sps = stab_CSS_sps + stability(X, rs1$fit_CSS_sps, rs2$fit_CSS_sps, "CSS") 
    stab_FSSS = stab_FSSS + stability(X, rs1$fit_FSSS, rs2$fit_FSSS, "FSSS")
    stab_Lasso = stab_Lasso + stability(X, rs1$fit_Lasso, rs2$fit_Lasso, "Lasso")
    stab_L0 = stab_L0 + stability(X, rs1$fit_L0, rs2$fit_L0, "L0")
  }
  cat("finished", "l", l, "out of", L, "\n")
}
stab = c(stab_SS, stab_CSS_wavg, stab_CSS_savg, stab_CSS_sps, stab_FSSS, stab_Lasso, stab_L0) / (L * (L-1)) * 2


saveRDS(c(rs_, stab), paste0("./results/", s0, "_", l__, ".RDS"))



# collect results ========
# collect all results in the folder "results" into a data frame called "Table1_l0.RDS"
# The column names are:
# [1] "s0"             "l"              "alpha_SS"       "alpha_CSS"      "alpha_FSSS"     "cutoff"        
# [7] "mse_SS"         "mse_CSS_wavg"   "mse_CSS_savg"   "mse_CSS_sps"    "mse_FSSS"       "mse_Lasso"     
# [13] "mse_L0"         "q_SS"           "q_CSS"          "q_FSSS"         "q_Lasso"        "q_L0"          
# [19] "FD_SS"          "FD_CSS_wavg"    "FD_CSS_savg"    "FD_CSS_sps"     "FD_FSSS"        "FD_Lasso"      
# [25] "FD_L0"          "PW_SS"          "PW_CSS_wavg"    "PW_CSS_savg"    "PW_CSS_sps"     "PW_FSSS"       
# [31] "PW_Lasso"       "PW_L0"          "inter_SS"       "inter_CSS_wavg" "inter_CSS_savg" "inter_CSS_sps" 
# [37] "inter_FSSS"     "inter_Lasso"    "inter_L0"       "stab_SS"        "stab_CSS_wavg"  "stab_CSS_savg" 
# [43] "stab_CSS_sps"   "stab_FSSS"      "stab_Lasso"     "stab_L0"          


