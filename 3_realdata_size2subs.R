rm(list = ls())
source("0_functions.R")
library(tidyverse)

# read in data
data2990 <- readRDS("data2990.RDS")
X = data2990$X
y = data2990$y

X = X - rep(1, nrow(X)) %*% t(colMeans(X))
X = X / ( rep(1, nrow(X)) %*% t(apply(X, 2, sd)) )
y = (y - mean(y)) / sd(y)
n = nrow(X); p = ncol(X)

bags = readRDS("realdata_bags.RDS")
Selection_set = readRDS("realdata_Selection_set.RDS")$Selection_set
pw = sapply(1:length(Selection_set), function(i) {
  length(Selection_set[[i]])
})
k = 45
Selection_set_focus = Selection_set[tensr:::topK(pw, k)]

alpha = 0.3


# exe ========
# repeat the following for i in 1:990
# take i=1 as an example
i = 1

combidx_set = combinat::combn(1:length(Selection_set_focus), 2)

zoom_set = combidx_set[,i]
set1 = Selection_set_focus[[zoom_set[1]]]
set2 = Selection_set_focus[[zoom_set[2]]]
combidx_feat1 = combinat::combn(set1, 2); combidx_feat1 = cbind(combidx_feat1, rep(1, 2) %*% t(set1))
combidx_feat2 = combinat::combn(set2, 2); combidx_feat2 = cbind(combidx_feat2, rep(1, 2) %*% t(set2))

for(j in 1:ncol(combidx_feat1)) {
  for(k in 1:ncol(combidx_feat2)) {
    zoom_feat1 = combidx_feat1[,j] %>% unique()
    zoom_feat2 = combidx_feat2[,k] %>% unique()
    TS1 = sort(zoom_feat1)
    TS2 = sort(zoom_feat2)
    
    if(all(TS1 == TS2)) next
    joint_supp = support(X, union(TS1, TS2), bags$base_lst)
    if(joint_supp >= 1 - alpha) next
    if(length(intersect(TS1, TS2)) != 0) {
      inter = TRUE
    } else {
      inter = FALSE
    }
    
    if(sum(TS1) < sum(TS2)) {
      S1 = TS1
      S2 = TS2
    } else {
      S1 = TS2
      S2 = TS1
    }
    
    rs = rep(0, 9)
    rs[1] = S1[1]
    rs[2] = S1[2]
    rs[3] = S2[1]
    rs[4] = S2[2]
    rs[5] = support(X, S1, bags$base_lst)
    rs[6] = support(X, S2, bags$base_lst)
    rs[7] = subs_u_wrt_S(X, y, S1, S2, Selection_set_focus) 
    rs[8] = nabla_tau_wrt_S(X, y, S1, S2, Selection_set_focus)$value
    rs[9] = inter
    
    filename = paste0("./subsets/subset_", i, "_", ".txt")
    write.table(t(rs), file = filename, append = TRUE, sep = ";", row.names = FALSE, col.names = FALSE)
    
    cat("finished k=", which(k == sample_id2), "out of", length(sample_id2), "\n")
  }
  cat("finished j=", which(j == sample_id1), "out of", length(sample_id1), "\n\n")
}

cat("finished i=", i, "out of", 1000, "\n")
cat("finished ...", "\n")



# collect results ========
# collect all results in the folder "subsets" into a txt called "realdata_Size2Subset.txt"
# The column names are:
# [1] "set11"       "set12"       "set21"       "set22"       "stab1"       "stab2"       "subs_u"      "nabla_value"
# [9] "inter"

