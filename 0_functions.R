# generate clusters
gen_clust = function(n, s, d, eta = 0.1) {
  X = matrix(rnorm(n*s*d), nrow = n)
  
  id_rep = (1:s) * d - (d-1)
  id_prox = setdiff(1:ncol(X), id_rep)
  track_rep = 0
  for(i in 1:ncol(X)) {
    if(i %in% id_rep) {
      # X[,i] = X[,i] / norm(X[,i], "2")
      track_rep = i
    } else {
      # X[,i] = X[,i] * eta
      # X[,i] = X[,track_rep] + X[,i]
      X[,i] = X[,track_rep] + rnorm(n, sd = eta)
    }
  }
  return(X)
}


# generate complex block
gen_block = function(n, nb, b, eta = 0.1) {
  X = matrix(rnorm(n*nb*(b+1)), nrow = n)
  id_child = seq(from = b+1, to = nb*(b+1), by = b+1)
  id_parents = setdiff(1:ncol(X), id_child)
  for(i in 1:ncol(X)) {
    if(i %in% id_parents) {
      # X[,i] = X[,i] / norm(X[,i], "2")
    } else {
      # X[,i] = X[,i] / norm(X[,i], "2") * eta
      # X[,i] = apply(X[,(i-b):(i-1)], 1, sum) + X[,i]
      X[,i] = apply(X[,(i-b):(i-1)], 1, sum) + rnorm(n, sd = eta)
    }
  }
  return(X)
}


# trace of a matrix
tr = function(X) {
  sum(diag(X))
}


# projection matrix of X
proj = function(X) {
  U = svd(X)$u
  U %*% t(U)
}


# find covering pairs on the poset
add_cover = function(idx, d) {
  candidate = setdiff(1:d, idx)
  lapply(candidate, function(i) {
    c(idx, i)
  })
}


# Gram Schmidt process
gs = function(new, orth) {
  terms = lapply(seq_len(ncol(orth)), function(i) {
    sum(new * orth[,i]) / sum(orth[,i]^2) *  matrix(orth[,i])
  })
  rs = new - apply( matrix( unlist(terms), ncol = ncol(orth) ), 1, sum)
  return(rs)
}


# subsampling using l0 base procedure
l0_subsampling <- function(X, y, s0, num_bags = 100){
  n = nrow(X)
  p = ncol(X)
  avg_ss <- matrix(0, p, 1)
  avg_sss <- matrix(0, n, n)
  Selset = matrix(0, 2*num_bags, p)
  base_lst = list()
  
  count = 1
  for (bag_iter in 1:num_bags){
    
    # complementary samples
    ind <- sample(1:n, n, replace = FALSE)
    ind_1 <- ind[1:(n/2)]
    ind_2 <- ind[as.numeric(n/2+1):n]
    
    # lasso model on first dataset: updating frequency of variables and subspaces
    lasso_mod <- L0Learn::L0Learn.fit(x = X[ind_1,], y = y[ind_1], maxSuppSize = s0, intercept = FALSE)
    v <- coef(lasso_mod)
    v = v[,ncol(v)]
    sel_var <- which(abs(v) != 0)
    Selset[count,] = v
    if (length(sel_var) >0){
      U = svd(X[,sel_var])$u
      base_lst = c(base_lst, list(U))
      avg_sss <- avg_sss + U %*% t(U)
    }
    count = count + 1
    
    # lasso model on second dataset: updating frequency of variables and subspaces
    lasso_mod <- L0Learn::L0Learn.fit(x = X[ind_2,], y = y[ind_2], maxSuppSize = s0, intercept = FALSE)
    v <- coef(lasso_mod)
    v = v[,ncol(v)]
    sel_var <- which(abs(v) != 0)
    Selset[count,] = v
    if (length(sel_var) >0){
      U = svd(X[,sel_var])$u
      base_lst = c(base_lst, list(U))
      avg_sss <- avg_sss + U %*% t(U)
    }
    count = count + 1
    
  }
  avg_ss <- avg_ss/(2*num_bags)
  avg_sss <- avg_sss/(2*num_bags)
  
  return(list(Selset = Selset,
              Sinfo = avg_sss,
              base_lst = base_lst))
}


# lasso with s0 as tuning parameter
lasso.s0 = function(X, y, maxSuppSize = s0) {
  fit_lasso = glmnet::glmnet(x = X, y = y)
  lasso_pw = which(sapply(1:ncol(fit_lasso$beta), function(i) {
    sum(fit_lasso$beta[,i] != 0)
  } ) <= s0)
  lambda_id = ifelse(length(lasso_pw) == 0, 1, max(lasso_pw))
  return(fit_lasso$beta[, lambda_id])
}


# subsampling using lasso base procedure
lasso_subsampling <- function(X, y, s0, num_bags = 100){
  n = nrow(X)
  p = ncol(X)
  avg_ss <- matrix(0, p, 1)
  avg_sss <- matrix(0, n, n)
  Selset = matrix(0, 2*num_bags, p)
  base_lst = list()
  
  count = 1
  for (bag_iter in 1:num_bags){
    
    # complementary samples
    ind <- sample(1:n, n, replace = FALSE)
    ind_1 <- ind[1:(n/2)]
    ind_2 <- ind[as.numeric(n/2+1):n]
    
    # lasso model on first dataset: updating frequency of variables and subspaces
    lasso_mod <- lasso.s0(X = X[ind_1,], y = y[ind_1], maxSuppSize = s0)
    sel_var <- which(abs(lasso_mod) != 0)
    Selset[count,] = lasso_mod
    if (length(sel_var) >0){
      U = svd(X[,sel_var])$u
      base_lst = c(base_lst, list(U))
      avg_sss <- avg_sss + U %*% t(U)
    }
    count = count + 1
    
    # lasso model on second dataset: updating frequency of variables and subspaces
    lasso_mod <- lasso.s0(X = X[ind_2,], y = y[ind_2], maxSuppSize = s0)
    sel_var <- which(abs(lasso_mod) != 0)
    Selset[count,] = lasso_mod
    if (length(sel_var) >0){
      U = svd(X[,sel_var])$u
      base_lst = c(base_lst, list(U))
      avg_sss <- avg_sss + U %*% t(U)
    }
    count = count + 1
    
  }
  avg_ss <- avg_ss/(2*num_bags)
  avg_sss <- avg_sss/(2*num_bags)
  
  return(list(Selset = Selset,
              Sinfo = avg_sss,
              base_lst = base_lst))
}


# selection proportion of each feature
selprob = function(bags) {
  apply(bags$Selset, 2, function(x) { mean(x!=0) })
}


# find subspace stability
find_sigmamin = function(base_lst, Z) {
  Uc = svd(Z)$u
  Mat = matrix(0, nrow = nrow(Z), ncol = nrow(Z))
  for(i in seq_along(base_lst)) {
    Ub = base_lst[[i]]
    prod = t(Uc) %*% Ub
    Mat = Mat + Uc %*% (prod %*% t(prod)) %*% t(Uc)
  }
  Mat = Mat / length(base_lst)
  min( RSpectra::svds(Mat, ncol(Z), 0, 0)$d )
}


# the greedy FSSS
greedy_FSSS = function(X, Sinfo, base_lst, node_par = c(), orth_par = NULL, alpha = 0.1) {
  d = ncol(X)
  while(TRUE) {
    if(length(node_par) == ncol(X)) {
      break
    }
    
    node_chil = add_cover(node_par, d = d)
    if (is.null(orth_par)) {
      dir_chil = lapply(node_chil, function(x) {
        X[,x]
      })
    } else {
      dir_chil = lapply(node_chil, function(x) { # find the gs-orthogonal vector
        gs(X[, setdiff(unlist(x), node_par)],
           orth_par)
      })
    }
    psi_chil = sapply(dir_chil, function(x) {
      if(norm(unlist(x), "2") < 1e-7) {
        return(1)
      } else {
        return( 1 - t(unlist(x)) %*% Sinfo %*% unlist(x) / norm(unlist(x), "2")^2 )
      }
    })
    idx = which.min(psi_chil)
    
    sigmamin = find_sigmamin(base_lst, X[,node_chil[[idx]], drop = FALSE])
    
    if( psi_chil[idx] <= alpha & sigmamin >= 1-alpha) {
      orth_par = cbind( orth_par, dir_chil[[idx]])
      node_par = node_chil[[idx]]
    } else {
      break
    }
  }
  
  return(list(
    NODES = node_par
  ))
}


# Cluster stability selection
cluSS_pred_alpha = function(Selset, cluster, alpha, wtype = "wavg") {
  feat_sel_props = apply(Selset != 0, 2, mean)
  c = length(unique(cluster))
  clu_sel_props = sapply(1:c, function(i) {
    idx = which(cluster == i)
    mean(apply(as.matrix(Selset[, idx]) != 0, 1, sum) !=0 )
  })
  Cluster = which(clu_sel_props >= 1 - alpha)
  Weights = list()
  for(i in Cluster) {
    idx = which(cluster == i)
    if(sum(feat_sel_props[idx]) == 0) {
      weight = rep(1, length(idx)) / length(idx)
    } else {
      weight = dplyr::case_when(
        wtype == "wavg" ~ feat_sel_props[idx] / sum(feat_sel_props[idx]),
        wtype == "savg" ~ rep(1, length(idx)) / length(idx),
        wtype == "sps" ~ feat_sel_props[idx] == max(feat_sel_props[idx]) / sum(feat_sel_props[idx] == max(feat_sel_props[idx]))
      )
    }
    Weights[[length(Weights) + 1]] = list(feature = idx, weight = weight)
  }
  return(list(
    Cluster = Cluster,
    Weights = Weights,
    clu_sel_props = clu_sel_props
  ))
}


# fit a linear model
lin_reg = function(X, y) {
  if(ncol(X) <= 1) {
    dat = data.frame(y = y, X = X)
    return( coef(lm(y~., data = dat)) )
  }
  else {
    fit = glmnet::cv.glmnet(x = X, y = y, alpha = 0)
    beta0 = fit$glmnet.fit$a0[which(fit$lambda == fit$lambda.min)]
    beta = fit$glmnet.fit$beta[, which(fit$lambda == fit$lambda.min)]
    return( as.numeric(c(beta0, beta)) )
  }
}


# create predictor from CSS results
create_X_cluster = function(fit_obj, X_new) {
  X_create = matrix(0, nrow = nrow(X_new), ncol = 0)
  for(i in seq_along(fit_obj$Cluster)) {
    idx = fit_obj$Weights[[i]]$feature
    weight = fit_obj$Weights[[i]]$weight
    if(length(fit_obj$Weights[[i]]$feature) == 1) {
      X_create = cbind(X_create, X_new[,idx])
    } else {
      X_create = cbind(X_create, X_new[,idx] %*% weight)
    }
  }
  return(X_create)
}


# compute MSE
MSE = function(X, model, y) {
  beta0 = model[1]
  beta = model[-1]
  yhat = beta0 + X %*% beta
  err = y - yhat
  return( norm(err,"2")^2 / nrow(X) )
}


# all-path FSSS
allpath_fsss = function(X, Sinfo, base_lst, node_par = c(), orth_par = NULL, alpha = 0.1) {
  d = ncol(X)
  stab = 0
  while(TRUE) {
    if(length(node_par) == ncol(X)) {
      break
    }
    
    node_chil = add_cover(node_par, d = d)
    if (is.null(orth_par)) {
      dir_chil = lapply(node_chil, function(x) {
        X[,x]
      })
    } else {
      dir_chil = lapply(node_chil, function(x) { # find the gs-orthogonal vector
        gs(X[, setdiff(unlist(x), node_par)],
           orth_par)
      })
    }
    psi_chil = sapply(dir_chil, function(x) {
      if(norm(unlist(x), "2") < 1e-7) {
        return(1)
      } else {
        return( 1 - t(unlist(x)) %*% Sinfo %*% unlist(x) / norm(unlist(x), "2")^2 )
      }
    })
    candidate = which(psi_chil < alpha)
    while(length(candidate) > 0) {
      idx = ifelse(length(candidate) != 1, sample(candidate, 1), candidate)
      sigmamin = find_sigmamin(base_lst, X[,node_chil[[idx]], drop = FALSE])
      if( sigmamin >= 1-alpha & psi_chil[idx] <= alpha ) {
        break
      } else {
        candidate = setdiff(candidate, idx)
      }
    }
    
    if(length(candidate) == 0) break
    if( psi_chil[idx] <= alpha & sigmamin >= 1-alpha) {
      orth_par = cbind( orth_par, dir_chil[[idx]])
      node_par = node_chil[[idx]]
      stab = sigmamin
    } else {
      break
    }
  }
  return(list(nodes = node_par, stab = stab))
}


# subspace stability of a set
support = function(X, S, base_lst) {
  find_sigmamin(base_lst, X[,S, drop = FALSE])
}


# the substitutability metric "tau"
subs_u = function(X, y, S1, S2, S0) {
  U = svd(X[,S0])$u
  U1 = svd(X[, union(S0, S1)])$u
  U2 = svd(X[, union(S0, S2)])$u
  u1 = U1 %*% (t(U1) %*% y) - U %*% (t(U) %*% y)
  u2 = U2 %*% (t(U2) %*% y) - U %*% (t(U) %*% y)
  
  u.vec = list(u1, u2)
  
  if(norm(u1, "2") < norm(u2, "2")) {
    u1 = u.vec[[2]]
    u2 = u.vec[[1]]
  }
  
  norm(proj(u1) %*% u2, "2") / norm(u1, "2")
}


# the feature perturbation matrix "nabla tau"
nabla_tau = function(X, y, S1, S2, S0) {
  subs.mat = matrix(0, nrow = length(S1), ncol = length(S2))
  for(i in 1:nrow(subs.mat)) {
    for(j in 1:ncol(subs.mat)) {
      subs.mat[i,j] = subs_u(X, y, S1[i], S2[j], S0)
    }
  }
  
  subs.val = c()
  for(i in 1:min(length(S1), length(S2))) {
    subs.val = c(subs.val, max(subs.mat))
    idx = which(subs.mat == max(subs.mat), arr.ind = TRUE)[1,]
    subs.mat = subs.mat[-idx[1], ,drop = FALSE]
    subs.mat = subs.mat[,-idx[2], drop = FALSE]
  }
  tau_overall = subs_u(X, y, S1, S2, S0)
  return( abs(min(subs.val) - tau_overall) / max(min(subs.val), tau_overall) )
}


# the tau metric wrt selection sets
subs_u_wrt_S = function(X, y, S1, S2, Selection_set) {
  Subs_U = c()
  for(i in seq_along(Selection_set)) {
    S0 = Selection_set[[i]]
    if(all(S1 %in% S0) | all(S2 %in% S0) ) {
      S0 = setdiff(Selection_set[[i]], union(S1, S2))
      
      Subs_U = c(Subs_U, subs_u(X, y, S1, S2, S0))
    }
  }
  return(min(Subs_U))
}


# the nabla tau metric wrt selection sets
nabla_tau_wrt_S = function(X, y, S1, S2, Selection_set) {
  Nabla_Tau = c()
  Idx = c()
  for(i in seq_along(Selection_set)) {
    S0 = Selection_set[[i]]
    if(all(S1 %in% S0) | all(S2 %in% S0) ) {
      S0 = setdiff(Selection_set[[i]], union(S1, S2))
      Nabla_Tau = c(Nabla_Tau, nabla_tau(X, y, S1, S2, S0))
      Idx = c(Idx, i)
    }
  }
  return(list(
    value = min(Nabla_Tau),
    idx = Idx[which.min(Nabla_Tau)]
  ))
}


# the metric used in creating radar chart
subs_u_wrt_S_radar = function(X, y, j, S, S_, Selection_set) {
  Subs_U = c()
  for(i in seq_along(Selection_set)) {
    S0 = Selection_set[[i]]
    if(all(S %in% S0) | all(S_ %in% S0) ) {
      S0 = setdiff(Selection_set[[i]], union(S, S_))
      
      Subs_U = c(Subs_U, subs_u(X, y, j, S, S0))
    }
  }
  return(max(Subs_U))
}


