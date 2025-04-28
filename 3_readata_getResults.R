rm(list = ls())
source("0_functions.R")
library(tidyverse)
library(latex2exp)
library(patchwork)
library(ggdendro)

# Figure 4 ========
# read in the results (collected from 3_realdata_l0.R)
Results = readRDS("RS_realdata.RDS") %>%
  mutate(
    q_CSS_wavg = q_CSS,
    q_CSS_savg = q_CSS,
  ) %>%
  rename(q_CSS_sps = q_CSS) %>%
  select(!ends_with(c("_Lasso", "_CSS_wavg", "_CSS_savg"))) 

mse.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(s0!=-2)

q.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("inter"), -starts_with("stab")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "q") %>%
  group_by(s0, method) %>%
  summarise(q = mean(q, na.rm = T)) %>%
  mutate(
    Method = str_remove(method, "q_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(s0!=-2)


stab.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "stab") %>%
  group_by(s0, method) %>%
  summarise(stab = mean(stab, na.rm = T)) %>%
  mutate(
    Method = str_remove(method, "stab_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(s0!=-2) 

pointsize = 1
p1 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(stab.tab) %>%
  ggplot(aes(x = mse, y = stab, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  labs(x = "Test MSE", y = "Output stability")  

p2 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(q.tab) %>%
  ggplot(aes(x = q, y = mse, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  labs(x = "Number of selected features", y = "Test MSE")

p3 = stab.tab %>%
  select(s0, Method, stab) %>%
  left_join(q.tab) %>%
  ggplot(aes(x = q, y = stab, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  labs(x = "Number of selected features", y = "Output stability")

fontsize = 8
p1 + p2 + p3 +
  plot_layout(ncol = 3, guides = "collect") & 
  theme(legend.position = "bottom", legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
        axis.title = element_text(size=fontsize), plot.title = element_text(size=fontsize),
        axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize),
  )




# focus on the dataset ========
# readin data
data2990 <- readRDS("data2990.RDS")
X = as.matrix(data2990$X)
y = data2990$y
gename = data2990$gename


X = X - rep(1, nrow(X)) %*% t(colMeans(X))
X = X / ( rep(1, nrow(X)) %*% t(apply(X, 2, sd)) )
y = (y - mean(y)) / sd(y)


# selection
s0 = 50
alpha = 0.3
set.seed(1234)
bags = l0_subsampling(X, y, s0, num_bags = 100)
# you may use this line to get the bags we used:
# bags = readRDS("realdata_bags.RDS")

Selection_set = list()
Stab = c()
while(length(Selection_set) < 100) {
  Selection_set_obj = allpath_fsss(X, bags$Sinfo, bags$base_lst, alpha = alpha)
  S = Selection_set_obj$nodes
  if(! (list(sort(S)) %in% Selection_set) ) {
    Selection_set = c(Selection_set, list(sort(S)))
    Stab = c(Stab, Selection_set_obj$stab)
  }
  cat("finished", length(Selection_set), "\n")
}

# you may use this line to get the selection set we used:
# Selection_set_obj = readRDS("realdata_Selection_set.RDS")
# Selection_set = Selection_set_obj$Selection_set
# Stab = Selection_set_obj$Stab


# find the largest 45 selection sets
pw = sapply(1:length(Selection_set), function(i) {
  length(Selection_set[[i]])
})
k = 45
Selection_set_focus = Selection_set[tensr:::topK(pw, k)]
Stab_focus = Stab[tensr:::topK(pw, k)]
# stability range
Stab_focus %>% range()



# Figure 5 ========
# read in the results (collected from 3_realdata_size1subs.R)
Size1Metrics = read.table("realdata_Size1Subset.txt", sep = ";")
Size1Metrics = Size1Metrics %>%
  `colnames<-`( 
    c("set1", "set2", "stab1", "stab2", "subs_u") ) %>%
  distinct(set1, set2, .keep_all = TRUE) 
fixed.set1 = Size1Metrics %>%
  arrange(desc(subs_u)) %>%
  filter(subs_u > 0.75) %>% pull(set1)
fixed.set2 = Size1Metrics %>%
  arrange(desc(subs_u)) %>%
  filter(subs_u > 0.75) %>% pull(set2)
fixed.set = sort(unique(union(fixed.set1, fixed.set2)))

d = length(fixed.set)
subs.mat = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    if(i < j) {
      rs = 1 - Size1Metrics$subs_u[ Size1Metrics$set1 == fixed.set[i] & Size1Metrics$set2 == fixed.set[j] ]
    } else if (i > j) {
      rs = 1 - Size1Metrics$subs_u[ Size1Metrics$set1 == fixed.set[j] & Size1Metrics$set2 == fixed.set[i] ]
    } else {
      rs = 0
    }
    if(length(rs) == 0) rs = 1
    subs.mat[i,j] = rs
  }
}
colnames(subs.mat) = gename[fixed.set]
rownames(subs.mat) = gename[fixed.set]

hc <- hclust(as.dist(subs.mat), method = "complete")
ggdendrogram(hc, rotate = FALSE, theme_dendro = FALSE) +
  geom_hline(yintercept = 0.2, linetype = "dashed", col = "red") +
  labs(x = "Gene expression features", 
       y = TeX(r"($1 - \tau(\cdot ; S)$)")
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),  # Rotate x-axis labels
        legend.position = "left")



# Figure 6 ========
Size2Metrics = read.table("realdata_Size2Subset.txt", sep = ";") 
# read in the results (collected from 3_realdata_size2subs.R)
tab = Size2Metrics %>%
  `colnames<-`( 
    c("set11", "set12", "set21", "set22", "stab1", "stab2", "subs_u", "nabla_value", "inter") ) %>%
  distinct(set11, set12, set21, set22, .keep_all = TRUE) %>%
  mutate(set11_ = gename[as.numeric(set11)], set12_ = ifelse(is.na(set12), "", gename[as.numeric(set12)]),
         set21_ = gename[as.numeric(set21)], set22_ = ifelse(is.na(set22), "", gename[as.numeric(set22)])
  ) %>%
  unite("set1", set11_, set12_, sep = " | ", remove = TRUE) %>%
  unite("set2", set21_, set22_, sep = " | ", remove = TRUE) %>%
  mutate(
    label = paste0("{", set1, ", ", set2, "}"),
    inter = as.logical(inter),
    nabla_thres = ifelse(subs_u > 0.7, nabla_value, NA)
  )


all_tri_filter = tab %>% 
  filter(nabla_value > 0.5, subs_u > 0.8) %>% 
  distinct(set11, set12, set21, set22, .keep_all = TRUE) %>%
  mutate(label = paste0("{", set1, ", ", set2, "}"))

Tau = matrix(0, nrow = nrow(all_tri_filter), ncol = 5)
for(i in 1:nrow(all_tri_filter)) {
  row_focus = all_tri_filter[i,]
  S1 = c(row_focus$set11, row_focus$set12)
  S2 = c(row_focus$set21, row_focus$set22)
  
  S = c(S1, S2)
  Tau.pair = c()
  for(j in S1) {
    rs = subs_u_wrt_S_radar(X, y, j, S1, S2, Selection_set_focus)
    Tau.pair = c(Tau.pair, rs)
  }
  for(j in S2) {
    rs = subs_u_wrt_S_radar(X, y, j, S2, S1, Selection_set_focus)
    Tau.pair = c(Tau.pair, rs)
  }
  Tau[i,] = c(Tau.pair, max(Tau.pair))
}

Tau = Tau %>%
  data.frame() %>%
  `colnames<-`(c("[{a}, S1]", "[{b}, S1]",
                 "[{c}, S2]", "[{d}, S2]", "max")) 

tab1 = all_tri_filter %>%
  cbind(Tau) %>%
  arrange(max) %>%
  slice(1:7) %>%
  select(label, starts_with("["), "max") %>%
  `rownames<-`({.$label}) %>%
  select(-label) 

tab1 =rbind(
  matrix(1, 1, 5) %>% `rownames<-`("1") %>% `colnames<-`(colnames(Tau)),
  matrix(0, 1, 5) %>% `rownames<-`("2") %>% `colnames<-`(colnames(Tau)),
  tab1
) 

opar <- par() 
par(mar = rep(1,4))
par(mfrow = c(2,4))
for (i in 3:nrow(tab1)) {
  fmsb::radarchart(
    tab1[c(1:2, i), ],
    vlcex = 0.8,
    pfcol = c("#99999980",NA),
    pcol= c(NA,2), plty = 1, plwd = 2,
    title = row.names(tab1)[i]
  )
}

