rm(list = ls())
source("0_functions.R")
library(tidyverse)
library(patchwork)


# functions ====
support_individual = function(X, Selset, type, cluster = NULL, Sinfo = NULL) {
  if(type == "SS") {
    return( apply(Selset != 0, 2, mean) )
  }
  if(type == "CSS") {
    feat_sel_props = apply(Selset != 0, 2, mean)
    clu_sel_props = sapply(1:length(cluster), function(i) {
      idx = which(cluster == cluster[i])
      mean(apply(as.matrix(Selset[, idx]) != 0, 1, sum) !=0 )
    })
    return(clu_sel_props)
  }
  if(type == "FSSS") {
    prob = sapply(1:ncol(X), function(i) {
      (t(X[,i]) %*% Sinfo %*% X[,i]) / norm(X[,i], "2")^2
    })
  }
  return(prob)
}


one_run_support = function(X, y, s0, cutoff = 0.1) {
  # subsampling
  bags = l0_subsampling(X, y, s0, num_bags = 50)
  # SS
  supp_SS = support_individual(X, bags$Selset, "SS")
  # CSS
  hclustobj = hclust(as.dist(1-abs(cor(X))), method = "single")
  cluster = as.vector(cutree(hclustobj, h = cutoff))  # cutoff is from cross-validation
  supp_CSS = support_individual(X, bags$Selset, "CSS", cluster = cluster)
  # FSSS
  supp_FSSS = support_individual(X, bags$Selset, "FSSS", Sinfo = bags$Sinfo)
  return(list(
    supp_SS = supp_SS,
    supp_CSS = supp_CSS,
    supp_FSSS = supp_FSSS
  ))
}


gen_block2 = function(n, s, eta = 0.1) {
  X = matrix(rnorm(n*s*4), nrow = n)
  beta = rep(0, 4*s)
  for(i in 1:s) {
    x1 = X[,4*i - 3]; beta[4*i - 3] = 1.5
    x2 = X[,4*i - 2]; beta[4*i - 2] = 1
    X[,4*i - 1] = x1 + x2 + rnorm(n, sd = eta)
    X[,4*i] = x1 - x2 + rnorm(n, sd = eta)
  }
  return(list(X = X, beta = beta))
}


gen_illust = function(eta = 0.2) {
  n = 100
  Xclust = gen_clust(n, s = 8, d = 3, eta = eta)
  Xblock = gen_block2(n, s = 2, eta = eta)
  Xnoise = matrix(rnorm(n*50), nrow = n)
  y = rowSums(Xclust[,3*(1:8) - 2]) + Xblock$X %*% Xblock$beta + rnorm(n, sd = eta)
  X = cbind(Xclust, Xblock$X, Xnoise); p = ncol(X)
  
  return(list(X = X, y = y))
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


# Figure 1 ====
cutoff = 0.1 
data = gen_illust()
# you may use this line to get the dataset we used:
# data = readRDS("illustration.RDS")
S0 = seq(from = 1, to = 40) 
for(s0 in S0) {
  for(l in 1:100) {
    tab = matrix(0, nrow = 82, ncol = 6)
    tab[,1] = s0
    tab[,2] = l
    tab[,3] = c(1:82)
    support_obj = one_run_support(data$X, data$y, s0 = s0, cutoff = cutoff)
    tab[,4:6] = cbind(support_obj$supp_SS,
                      support_obj$supp_CSS,
                      support_obj$supp_FSSS)
    
    write.table(tab, file = "illustrate.txt", append = TRUE, sep = ";", row.names = FALSE, col.names = FALSE)
    cat("finished","s0 =", s0, ", l =", l, "\n" )
  }
}


# you may use this line to get the `tab` we used:
# tab = read.table("illustrate.txt", sep = ";")
pointsize = 1
fontsize = 9
tab %>% 
  `colnames<-`(c("s0", "l", "feat", "Stability selection", "Cluster stability selection", "Proposed subspace framework")) %>%
  pivot_longer(-c(s0, l, feat), names_to = "method", values_to = "stab") %>%
  group_by(s0, feat, method) %>%
  summarise(stab = mean(stab)) %>%
  mutate(
    Variable = case_when(
      feat %in% c(3*(1:8) - 2,  21 + (1:2) * 4, 22 + (1:2) * 4 ) ~ "Signal",
      feat %in% c(3*(1:8)-1, 3*(1:8), 23 + (1:2) * 4, 24 + (1:2) * 4) ~ "Correlated signal",
      feat %in% c(33:82) ~ "Noise"),
    Variable = fct_inorder(Variable),
    method = factor(method, levels = c("Stability selection", "Cluster stability selection", "Proposed subspace framework"))
  ) %>%
  ggplot(aes(x = s0, y = stab, group = feat, col = Variable)) +
  geom_line(position=position_jitter(w=0.02, h=0)) +
  scale_color_manual(
    values = c("Noise" = "#999999", "Signal" = "#f0270c", "Correlated signal" = "#62c4f5", "Block correlated signal" = "#FF8000") 
  ) +
  facet_grid(cols = vars(method), switch = "y") +
  labs(x = "Tuning parameter", y = "Stability") +
  theme_bw() +
  theme(
    legend.position="bottom",
    legend.title = element_text(size=pointsize), legend.text = element_text(size=pointsize), 
    axis.title = element_text(size=pointsize), plot.title = element_text(size=pointsize),
    axis.text.x = element_text(size=pointsize), axis.text.y = element_text(size=pointsize)
  ) 



# Figure 11 ====
train_set = gen_illust()
# you may use this line to get the dataset we used:
# train_set = readRDS("illustration.RDS")
X_train = train_set$X
y_train = train_set$y
test_set = gen_illust()
X_test = test_set$X
y_test = test_set$y

beta = c(rep(c(1,0,0), 8), rep(c(1.5, 1, 0, 0),2), rep(0, 50))


# subsampling
L = 100
RS = matrix(0, nrow = 0, ncol = 13)
for(s0 in seq(from = 1, to = 40) ) {
  STORE = list()
  RS_inner = matrix(0, nrow = 0, ncol = 9)
  for(l in 1:L) {
    # subsampling
    bags = l0_subsampling(X_train, y_train, s0, num_bags = 50)
    
    # selection
    # SS
    fit_SS = which(selprob(bags) > 0.8)
    # CSS
    hclustobj = hclust(as.dist(1-abs(cor(X_train))), method = "single")
    cluster = as.vector(cutree(hclustobj, h = 0.1))  # cutoff is from cross-validation
    fit_CSS_sps = cluSS_pred_alpha(bags$Selset, cluster, alpha = 0.2, "sps")
    # FSSS
    fit_FSSS = greedy_FSSS(X_train, bags$Sinfo, bags$base_lst, alpha = 0.2)
    
    # modeling
    model_SS = lin_reg(X_train[, fit_SS, drop = FALSE], y_train)
    X_CSS_sps = create_X_cluster(fit_CSS_sps, X_train)
    model_CSS_sps = lin_reg(X_CSS_sps, y_train)
    model_FSSS = lin_reg(X_train[, fit_FSSS$NODES, drop = FALSE], y_train)
    
    # testing
    # calculate MSE
    mse_SS = MSE(X_test[, fit_SS, drop = FALSE], model_SS, y_test)
    X_CSS_sps = create_X_cluster(fit_CSS_sps, X_test)
    mse_CSS_sps = MSE(X_CSS_sps, model_CSS_sps, y_test)
    mse_FSSS = MSE(X_test[, fit_FSSS$NODES, drop = FALSE], model_FSSS, y_test)
    # find PS
    PStar = proj(X_test[, beta != 0])
    PSperp = diag(rep(1, nrow(X_test))) - PStar
    # calculate FD
    FD_SS = ifelse(length(fit_SS)==0, 0, tr(PSperp %*% proj(X_test[, fit_SS, drop = FALSE])) )
    FD_CSS_sps = ifelse(length(fit_CSS_sps$Cluster) == 0, 0, tr(PSperp %*% proj(X_CSS_sps)) ) 
    FD_FSSS = ifelse(length(fit_FSSS$NODES) == 0, 0, tr(PSperp %*% proj(X_test[, fit_FSSS$NODES, drop = FALSE])) ) 
    # calculate PW
    PW_SS = ifelse(length(fit_SS)==0, 0, tr(PStar %*% proj(X_test[, fit_SS, drop = FALSE])) )
    PW_CSS_sps = ifelse(length(fit_CSS_sps$Cluster) == 0, 0, tr(PStar %*% proj(X_CSS_sps)) )
    PW_FSSS = ifelse(length(fit_FSSS$NODES) == 0, 0, tr(PStar %*% proj(X_test[, fit_FSSS$NODES, drop = FALSE])) )
    
    rs = list(
      fit_SS = fit_SS,
      fit_CSS_sps = fit_CSS_sps,
      fit_FSSS = fit_FSSS
    )
    STORE = c(STORE, list(rs))
    RS_inner = rbind(RS_inner, c(mse_SS, mse_CSS_sps, mse_FSSS, FD_SS, FD_CSS_sps, FD_FSSS, PW_SS, PW_CSS_sps, PW_FSSS))
    cat("finished l =", l, "\n")
  }
  # summary
  summary = colMeans(RS_inner, na.rm = TRUE)
  # robustness
  stab_SS = 0
  stab_CSS_sps = 0
  stab_FSSS = 0
  for(l in 1:L) {
    if (l == L)  next
    for(l_ in (l+1):L ) {
      rs1 = STORE[[l]]
      rs2 = STORE[[l_]]
      stab_SS = stab_SS + stability(X, rs1$fit_SS, rs2$fit_SS, "SS") 
      stab_CSS_sps = stab_CSS_sps + stability(X, rs1$fit_CSS_sps, rs2$fit_CSS_sps, "CSS") 
      stab_FSSS = stab_FSSS + stability(X, rs1$fit_FSSS$NODES, rs2$fit_FSSS$NODES, "FSSS")
    }
    cat("finished", "l", l, "out of", L, "\n")
  }
  stab = c(stab_SS, stab_CSS_sps, stab_FSSS) / (L * (L-1)) * 2
  
  RS = rbind(RS, c(s0, summary, stab))
  cat("finished", s0, "\n\n")
}

odd_2_even = function(x) {
  (x - 1) %/% 2 * 2 + 1
}


# you may use this line to get the dataset we used:
# RS = readRDS("RS_illustration.RDS")
Results = RS %>% 
  data.frame() %>%
  `colnames<-`(c("s0", "mse_SS", "mse_CSS", "mse_FSSS", "FD_SS", "FD_CSS", "FD_FSSS", "PW_SS", "PW_CSS", "PW_FSSS",
                 "stab_SS", "stab_CSS", "stab_FSSS")) %>%
  mutate(s0 = odd_2_even(s0) ) 

mse.tab = Results %>%
  select(-starts_with("FD"), -starts_with("PW"), -starts_with("stab")) %>%
  pivot_longer(!c("s0"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse, na.rm = T)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("SS", "CSS", "FSSS"), labels = c("SS", "CSS (sps)", "FSSS"))
  ) 

FD.tab = Results %>%
  select(-starts_with("mse"), -starts_with("PW"), -starts_with("stab")) %>%
  pivot_longer(!c("s0"), names_to = "method", values_to = "FD") %>%
  group_by(s0, method) %>%
  summarise(FD = mean(FD, na.rm = T)) %>%
  mutate(
    Method = str_remove(method, "FD_"),
    Method = factor(Method, levels = c("SS", "CSS", "FSSS"), labels = c("SS","CSS (sps)", "FSSS"))
  ) 

PW.tab = Results %>%
  select(-starts_with("mse"), -starts_with("FD"), -starts_with("stab")) %>%
  pivot_longer(!c("s0"), names_to = "method", values_to = "PW") %>%
  group_by(s0, method) %>%
  summarise(PW = mean(PW, na.rm = T)) %>%
  mutate(
    Method = str_remove(method, "PW_"),
    Method = factor(Method, levels = c("SS", "CSS", "FSSS"), labels = c("SS","CSS (sps)", "FSSS"))
  ) 

stab.tab = Results %>%
  select(-starts_with("mse"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0"), names_to = "method", values_to = "stab") %>%
  group_by(s0, method) %>%
  summarise(stab = mean(stab, na.rm = T)) %>%
  mutate(
    Method = str_remove(method, "stab_"),
    Method = factor(Method, levels = c("SS", "CSS", "FSSS"), labels = c("SS","CSS (sps)", "FSSS"))
  ) 


p1 = mse.tab %>%
  ggplot(aes(x = s0, y = mse, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  geom_line() +
  theme_bw() +
  labs(x = "Tuning parameter", y = "Test MSE") 


p2 = FD.tab %>%
  ggplot(aes(x = s0, y = FD, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  geom_line() +
  theme_bw() +
  labs(x = "Tuning parameter", y = "FPE") 


p3 = PW.tab %>%
  ggplot(aes(x = s0, y = PW, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  geom_line() +
  theme_bw() +
  labs(x = "Tuning parameter", y = "TP") 

p4 = stab.tab %>%
  ggplot(aes(x = s0, y = stab, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  geom_line() +
  theme_bw() +
  labs(x = "Tuning parameter", y = "Output stability") 


p1 + p2 + p3  + p4 +
  plot_annotation(title = "Base procedure: L0") +
  plot_layout(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom", legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
        axis.title = element_text(size=fontsize), plot.title = element_text(size=9),
        axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize),
  )




