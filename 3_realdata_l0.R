rm(list = ls())
source("0_functions.R")
library(tidyverse)

# repeat the following for s0 in S0 = seq(from = 22, to = 100, by = 2) 
# repeat each s0 500 times

# example of s0 = 30 and ll = 1 (ll is the index for repetition)
s0 = 30
ll = 1

# read in data
data2990 <- readRDS("data2990.RDS")
X = data2990$X
y = data2990$y

n_test = 30
id_test = sample(1:nrow(X), n_test)
id_train = setdiff(1:nrow(X), id_test)
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
  
  
  X_train = X_train - rep(1, nrow(X_train)) %*% t(colMeans(X_train))
  y_train = y_train - mean(y_train)
  X_test = X_test - rep(1, nrow(X_test)) %*% t(colMeans(X_test))
  y_test = y_test - mean(y_test)
  
  
  # divide the data
  n = length(y_train)
  n0 = round(n * 9/10)
  feature_id = sample(1:n, n0)
  test_id = setdiff(1:n, feature_id)
  cat("finished deviding data", "\n")
  
  
  # cluster
  nontest_X = X_train[-test_id,]
  hclustobj = hclust(as.dist(1-abs(cor(nontest_X))), method = "complete")
  # plot(hclustobj)
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
    model_SS = lin_reg(feature_X[, fit_SS, drop = FALSE], feature_y)
    X_CSS_wavg = create_X_cluster(fit_CSS_wavg, feature_X)
    model_CSS_wavg = lin_reg(X_CSS_wavg, feature_y)
    X_CSS_savg = create_X_cluster(fit_CSS_savg, feature_X)
    model_CSS_savg = lin_reg(X_CSS_savg, feature_y)
    X_CSS_sps = create_X_cluster(fit_CSS_sps, feature_X)
    model_CSS_sps = lin_reg(X_CSS_sps, feature_y)
    model_FSSS = lin_reg(feature_X[, fit_FSSS$NODES, drop = FALSE], feature_y)
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
  
  
  X_train = X_train - rep(1, nrow(X_train)) %*% t(colMeans(X_train))
  y_train = y_train - mean(y_train)
  X_test = X_test - rep(1, nrow(X_test)) %*% t(colMeans(X_test))
  y_test = y_test - mean(y_test)
  
  
  # divide the data
  n = length(y_train)
  n0 = round(n * 9/10)
  feature_id = sample(1:n, n0)
  test_id = setdiff(1:n, feature_id)
  cat("finished deviding data", "\n")
  
  
  # cluster
  nontest_X = X_train[-test_id,]
  hclustobj = hclust(as.dist(1-abs(cor(nontest_X))), method = "complete")
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
  model_SS = lin_reg(feature_X[, fit_SS, drop = FALSE], feature_y)
  X_CSS_wavg = create_X_cluster(fit_CSS_wavg, feature_X)
  model_CSS_wavg = lin_reg(X_CSS_wavg, feature_y)
  X_CSS_savg = create_X_cluster(fit_CSS_savg, feature_X)
  model_CSS_savg = lin_reg(X_CSS_savg, feature_y)
  X_CSS_sps = create_X_cluster(fit_CSS_sps, feature_X)
  model_CSS_sps = lin_reg(X_CSS_sps, feature_y)
  model_FSSS = lin_reg(feature_X[, fit_FSSS$NODES, drop = FALSE], feature_y)
  cat("finished OLS", "\n")
  
  
  # lasso 
  fit_lasso = glmnet::glmnet(x = feature_X, y = feature_y)
  lasso_pw = which(sapply(1:ncol(fit_lasso$beta), function(i) {
    sum(fit_lasso$beta[,i] != 0)
  } ) <= s0)
  lambda_id = ifelse(length(lasso_pw) == 0, 1, max(lasso_pw))
  fit_Lasso = as.numeric(which(fit_lasso$beta[, lambda_id] != 0 ))
  model_Lasso = lin_reg(feature_X[, fit_Lasso, drop = FALSE], feature_y)
  
  
  # l0
  fit_l0 = L0Learn::L0Learn.fit(x = feature_X, y = feature_y, maxSuppSize = s0, intercept = TRUE)
  beta_l0 = coef(fit_l0)
  beta_l0 = beta_l0[2:nrow(beta_l0), ncol(beta_l0)]
  fit_L0 = as.numeric( which(beta_l0 != 0) )
  model_L0 = lin_reg(feature_X[, fit_L0, drop = FALSE], feature_y)
  
  
  return_list = list(
    fit_SS = fit_SS,
    fit_CSS_wavg = fit_CSS_wavg,
    fit_CSS_savg = fit_CSS_savg,
    fit_CSS_sps = fit_CSS_sps,
    fit_FSSS = fit_FSSS,
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


validation = function(data, s0, Cutoff = c(0.1, 0.3, 0.5), Alpha = seq(from = 0.05, to = 0.3, by = 0.05), Lcv = 10) {
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
  mse_FSSS = MSE(X_test[, fit_FSSS$NODES, drop = FALSE], model_FSSS, y_test)
  mse_Lasso = MSE(X_test[, fit_Lasso, drop = FALSE], model_Lasso, y_test)
  mse_L0 = MSE(X_test[, fit_L0, drop = FALSE], model_L0, y_test)
  # nubmer of selected features
  q_SS = length(fit_SS)
  q_CSS = length(fit_CSS_wavg$Cluster)
  q_FSSS = length(fit_FSSS$NODES)
  q_Lasso = length(fit_Lasso)
  q_L0 = length(fit_L0)
  
  # interpretability
  inter_SS = ifelse(length(fit_SS) == 0, 0, 1)
  inter_FSSS = ifelse(length(fit_FSSS$NODES) == 0, 0, 1)
  inter_CSS_wavg = interpretability(fit_CSS_wavg)
  inter_CSS_savg = interpretability(fit_CSS_savg)
  inter_CSS_sps = ifelse(length(fit_CSS_sps$Cluster) == 0, 0, 1)
  inter_Lasso = ifelse(length(fit_Lasso) == 0, 0, 1)
  inter_L0 = ifelse(length(fit_L0) == 0, 0, 1)
  
  
  return(list(
    summary = c(
      alpha, cutoff, 
      mse_SS, mse_CSS_wavg, mse_CSS_savg, mse_CSS_sps, mse_FSSS, mse_Lasso, mse_L0,
      q_SS, q_CSS, q_FSSS, q_Lasso, q_L0,
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
  
  if(type %in% c("SS", "L0", "Lasso")) {
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
  
  if(type == "FSSS") {
    return( stab_metric(fit1$NODES, fit2$NODES, X) )
  }
  
}


interpretability = function(fit_CSS) {
  if(length(fit_CSS$Cluster) == 0) return(0)
  mean(1 / sapply(seq_along(fit_CSS$Weights), function(i) { 
    sum(fit_CSS$Weights[[i]]$weight != 0)
  } ))
}




# execute ========
validation_obj = validation(data, s0, Cutoff = c(0.1, 0.3, 0.5), Alpha = seq(from = 0.05, to = 0.3, by = 0.05))
L = 50
RS = matrix(0, nrow = 0, ncol = 23)
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


saveRDS(c(rs_, stab), paste0("./results1/", s0, "_", ll, ".RDS"))
cat("finished ...", "\n")




# collect results ========
# collect all results in the folder "results1" into a data frame called "RS_realdata.RDS"
# The column names are:
# [1] "s0"             "l"              "alpha_SS"       "alpha_CSS"      "alpha_FSSS"     "cutoff"        
# [7] "mse_SS"         "mse_CSS_wavg"   "mse_CSS_savg"   "mse_CSS_sps"    "mse_FSSS"       "mse_Lasso"     
# [13] "mse_L0"         "q_SS"           "q_CSS"          "q_FSSS"         "q_Lasso"        "q_L0"          
# [19] "inter_SS"       "inter_CSS_wavg" "inter_CSS_savg" "inter_CSS_sps"  "inter_FSSS"     "inter_Lasso"   
# [25] "inter_L0"       "stab_SS"        "stab_CSS_wavg"  "stab_CSS_savg"  "stab_CSS_sps"   "stab_FSSS"     
# [31] "stab_Lasso"     "stab_L0"         

