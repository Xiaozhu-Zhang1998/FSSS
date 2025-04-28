rm(list = ls())
source("0_functions.R")
library(tidyverse)

# read in data
data2990 <- readRDS("data2990.RDS")
X = as.matrix(data2990$X)
y = data2990$y

X = X - rep(1, nrow(X)) %*% t(colMeans(X))
X = X / ( rep(1, nrow(X)) %*% t(apply(X, 2, sd)) )
y = (y - mean(y)) / sd(y)


bags = readRDS("realdata_bags.RDS")
Selection_set = readRDS("realdata_Selection_set.RDS")$Selection_set
pw = sapply(1:length(Selection_set), function(i) {
  length(Selection_set[[i]])
})
k = 45
Selection_set_focus = Selection_set[tensr:::topK(pw, k)]

alpha = 0.3

candidate = unique(unlist(Selection_set_focus))
combidx_set = combinat::combn(candidate, 2)

for(i in 1:ncol(combidx_set)) {
  zoom = combidx_set[,i]
  TS1 = zoom[1]
  TS2 = zoom[2]
  
  joint_supp = support(X, union(TS1, TS2), bags$base_lst)
  if(joint_supp >= 1 - alpha) next

  if(sum(TS1) < sum(TS2)) {
    S1 = TS1
    S2 = TS2
  } else {
    S1 = TS2
    S2 = TS1
  }
  
  rs = rep(0, 5)
  rs[1] = paste(S1, collapse = "|")
  rs[2] = paste(S2, collapse = "|")
  rs[3] = support(X, S1, bags$base_lst)
  rs[4] = support(X, S2, bags$base_lst)
  rs[5] = subs_u_wrt_S(X, y, S1, S2, Selection_set_focus) 
  
  filename = paste0("realdata_Size1Subset.txt")
  write.table(t(rs), file = filename, append = TRUE, sep = ";", row.names = FALSE, col.names = FALSE)
  
  cat("finished i=", i, "out of", ncol(combidx_set), "\n")
}


