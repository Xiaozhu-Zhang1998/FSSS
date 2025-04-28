rm(list = ls())
source("0_functions.R")
library(tidyverse)
library(patchwork)
library(dichromat)


# l0 base procedure =====
# read in the simulation results (collected from 2_synthetic_l0.R)
Results = readRDS("RS_Syn_l0.RDS")
Results = Results %>%
  mutate(
    q_CSS_wavg = q_CSS,
    q_CSS_savg = q_CSS,
  ) %>%
  rename(q_CSS_sps = q_CSS) %>%
  select(-ends_with("_CSS_savg"), -ends_with("_CSS_wavg")) 


pointsize = 1.2
mse.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "Lasso")


fd.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "FD") %>%
  group_by(s0, method) %>%
  summarise(FD = mean(FD)) %>%
  mutate(
    Method = str_remove(method, "FD_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "Lasso") 


pw.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "PW") %>%
  group_by(s0, method) %>%
  summarise(PW = mean(PW)) %>%
  mutate(
    Method = str_remove(method, "PW_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "Lasso")


stab.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "stab") %>%
  group_by(s0, method) %>%
  summarise(stab = mean(stab)) %>%
  mutate(
    Method = str_remove(method, "stab_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "Lasso")


q.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "q") %>%
  group_by(s0, method) %>%
  summarise(q = mean(q)) %>%
  mutate(
    Method = str_remove(method, "q_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "Lasso")


p1 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(fd.tab) %>%
  ggplot(aes(x = mse, y = FD, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0.1,1), ylim = c(0, 5)) +
  labs(x = "Test MSE", y = "FPE")  


p2 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(pw.tab) %>%
  ggplot(aes(x = mse, y = PW, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0.1, 0.5), ylim = c(5, 16)) +
  labs(x = "Test MSE", y = "TP")  


p3 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(stab.tab) %>%
  ggplot(aes(x = mse, y = stab, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0.1, 0.75), ylim = c(0.5, 1)) +
  labs(x = "Test MSE", y = "Output stability")  


p4 = fd.tab %>%
  select(s0, Method, FD) %>%
  left_join(pw.tab) %>%
  ggplot(aes(x = FD, y =PW, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 15)) +
  labs(x = "FPE", y = "TP")  


## Figure 12 ====
fontsize = 7
p1 + p2 + p3 + p4 +
  plot_annotation(title = "Base procedure: L0") +
  plot_layout(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom", legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
        axis.title = element_text(size=fontsize), plot.title = element_text(size=9),
        axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize),
  )


## Left column of Table 1 =====
mse.min = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("Lasso", "CSS (wavg)", "CSS (savg)" ))) %>%
  group_by(Method) %>%
  summarise(mse.min = min(mse)) %>%
  deframe()

# s0.cv
s0.cv = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("Lasso", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "L0" & mse == mse.min["L0"]) | (Method == "SS" & mse == mse.min["SS"]) | (Method == "CSS (sps)" & mse == mse.min["CSS (sps)"]) | (Method == "FSSS" & mse == mse.min["FSSS"]) ) %>%
  select(Method,s0) %>%
  deframe()

# mse.cv 
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(mse), mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("Lasso", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "L0" & s0 == s0.cv["L0"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 


# FD.cv
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "FD") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(FD), FD = mean(FD))  %>%
  mutate(
    Method = str_remove(method, "FD_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("Lasso", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "L0" & s0 == s0.cv["L0"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 


# TP.cv 
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "PW") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(PW), PW = mean(PW)) %>%
  mutate(
    Method = str_remove(method, "PW_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("Lasso", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "L0" & s0 == s0.cv["L0"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 

# robust.cv
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "stab") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(stab), stab = mean(stab)) %>%
  mutate(
    Method = str_remove(method, "stab_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("Lasso", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "L0" & s0 == s0.cv["L0"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 



# lasso base procedure ====
# read in the simulation results (collected from 2_synthetic_lasso.R)
Results = readRDS("RS_Syn_lasso.RDS")
Results = Results %>%
  mutate(
    q_CSS_wavg = q_CSS,
    q_CSS_savg = q_CSS,
  ) %>%
  rename(q_CSS_sps = q_CSS) %>%
  select(-ends_with("_CSS_savg"), -ends_with("_CSS_wavg")) 

pointsize = 1.2

mse.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "L0") 


fd.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "FD") %>%
  group_by(s0, method) %>%
  summarise(FD = mean(FD)) %>%
  mutate(
    Method = str_remove(method, "FD_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "L0")


pw.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "PW") %>%
  group_by(s0, method) %>%
  summarise(PW = mean(PW)) %>%
  mutate(
    Method = str_remove(method, "PW_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "L0")


stab.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "stab") %>%
  group_by(s0, method) %>%
  summarise(stab = mean(stab)) %>%
  mutate(
    Method = str_remove(method, "stab_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "L0")


q.tab = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "q") %>%
  group_by(s0, method) %>%
  summarise(q = mean(q)) %>%
  mutate(
    Method = str_remove(method, "q_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(Method != "L0")


p5 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(fd.tab) %>%
  ggplot(aes(x = mse, y = FD, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0.1,1), ylim = c(0, 5)) +
  labs(x = "Test MSE", y = "FPE")  


p6 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(pw.tab) %>%
  ggplot(aes(x = mse, y = PW, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0.1, 0.5), ylim = c(5, 16)) +
  labs(x = "Test MSE", y = "TP")  


p7 = mse.tab %>%
  select(s0, Method, mse) %>%
  left_join(stab.tab) %>%
  ggplot(aes(x = mse, y = stab, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0.1, 0.75), ylim = c(0.5, 1)) +
  labs(x = "Test MSE", y = "Output stability")  


p8 = fd.tab %>%
  select(s0, Method, FD) %>%
  left_join(pw.tab) %>%
  ggplot(aes(x = FD, y =PW, col = Method, shape = Method)) +
  geom_point(size = pointsize) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 15)) +
  labs(x = "FPE", y = "TP")  


## Figure 13 ====
p5 + p6 + p7  + p8 +
  plot_annotation(title = "Base procedure: Lasso") +
  plot_layout(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom", legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
        axis.title = element_text(size=fontsize), plot.title = element_text(size=9),
        axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize),
  )




## Right column of Table 1 ====
mse.min = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("L0", "CSS (wavg)", "CSS (savg)" ))) %>%
  group_by(Method) %>%
  summarise(mse.min = min(mse)) %>%
  deframe()

# s0.cv
s0.cv = Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("L0", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "Lasso" & mse == mse.min["Lasso"]) | (Method == "SS" & mse == mse.min["SS"]) | (Method == "CSS (sps)" & mse == mse.min["CSS (sps)"]) | (Method == "FSSS" & mse == mse.min["FSSS"]) ) %>%
  select(Method, s0) %>%
  deframe()

# mse.cv 
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "mse") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(mse), mse = mean(mse)) %>%
  mutate(
    Method = str_remove(method, "mse_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("L0", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "Lasso" & s0 == s0.cv["Lasso"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 


# FD.cv
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "FD") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(FD), FD = mean(FD))  %>%
  mutate(
    Method = str_remove(method, "FD_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("L0", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "Lasso" & s0 == s0.cv["Lasso"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 


# TP.cv 
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("stab"), -starts_with("FD")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "PW") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(PW), PW = mean(PW)) %>%
  mutate(
    Method = str_remove(method, "PW_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("L0", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "Lasso" & s0 == s0.cv["Lasso"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 

# robust.cv
Results %>%
  select(-starts_with("alph"), -cutoff, -starts_with("mse"), -starts_with("q"), -starts_with("inter"), -starts_with("FD"), -starts_with("PW")) %>%
  pivot_longer(!c("s0", "l"), names_to = "method", values_to = "stab") %>%
  group_by(s0, method) %>%
  summarise(sd = sd(stab), stab = mean(stab)) %>%
  mutate(
    Method = str_remove(method, "stab_"),
    Method = factor(Method, levels = c("L0", "Lasso", "SS", "CSS_wavg", "CSS_savg", "CSS_sps", "FSSS"), labels = c("L0", "Lasso", "SS", "CSS (wavg)", "CSS (savg)", "CSS (sps)", "FSSS"))
  ) %>%
  filter(!(Method %in% c("L0", "CSS (wavg)", "CSS (savg)" ))) %>%
  filter((Method == "Lasso" & s0 == s0.cv["Lasso"]) | (Method == "SS" & s0 == s0.cv["SS"]) | (Method == "CSS (sps)" & s0 == s0.cv["CSS (sps)"]) | (Method == "FSSS" & s0 == s0.cv["FSSS"]) ) 




# Focus on one dataset and implement FSSS =====
# readin data
synData <- readRDS("synData.RDS")
X = synData$X
y = synData$y
beta = synData$beta

X = X - rep(1, nrow(X)) %*% t(colMeans(X))
X = X / ( rep(1, nrow(X)) %*% t(apply(X, 2, sd)) )
y = (y - mean(y)) / sd(y)

X = X[1:600,]
y = y[1:600]

# subsampling and selection
s0 = 35
alpha = 0.3
bags = l0_subsampling(X, y, s0, num_bags = 100)
Selection_set = list()
Stab = c()
while(length(Selection_set) < 100) {
  Selection_set_obj = allpath_fsss(X, bags$Sinfo, bags$base_lst, alpha = alpha)
  S = Selection_set_obj$nodes
  if(! (list(sort(S)) %in% Selection_set)) {
    Selection_set = c(Selection_set, list(sort(S)))
    Stab = c(Stab, Selection_set_obj$stab)
    cat(sort(S), "\n")
  }
  cat("finished", length(Selection_set), "\n")
}


# the interest_sets need to be stable and appear in at least one selection set
interest_sets = list(
  c(1,4,7),
  c(2,5,8),
  c(3,6,9),
  c(10, 11), 
  c(10, 12), 
  c(11, 12),
  c(13, 14, 15),
  c(14, 15, 16),
  c(13, 16),
  c(17, 18, 19, 20),
  c(18, 19, 20, 21),
  c(17, 18, 21),
  c(19, 21)
)

set.levels = sapply(interest_sets, function(x) {
  paste0(unlist(x), collapse = "|")
})

combidx_set = combinat::combn(1:length(interest_sets), 2)

RS = matrix(0, nrow = ncol(combidx_set), ncol = 7)
for(i in 1:ncol(combidx_set)) {
  zoom = combidx_set[,i]
  S1 = interest_sets[[zoom[1]]]
  S2 = interest_sets[[zoom[2]]]
  joint_supp = support(X, union(S1, S2), bags$base_lst)
  if(joint_supp >= 1 - alpha){
    next
  }
  rs = rep(0, 7)
  rs[1] = set.levels[zoom[1]]
  rs[2] = set.levels[zoom[2]]
  rs[3] = support(X, S1, bags$base_lst)
  rs[4] = support(X, S2, bags$base_lst)
  rs[5] = subs_u_wrt_S(X, y, S1, S2, Selection_set)  
  nabla_obj = nabla_tau_wrt_S(X, y, S1, S2, Selection_set)
  rs[6] = nabla_obj$value
  rs[7] = nabla_obj$idx
  RS[i, ]= rs
  cat("finished", i, "\n")
}


tab = RS %>%
  data.frame() %>%
  `colnames<-`( 
    c("set1", "set2", "stab1", "stab2", "subs_u", "nabla_value", "nabla_index") ) %>%
  mutate(
    across(c(stab1, stab1, subs_u, nabla_value), as.numeric)
  )


## order by dendrogram 
d = length(set.levels)
subs.mat = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    rs = 1 - tab$subs_u[ tab$set1 == set.levels[i] & tab$set2 == set.levels[j] ]
    if(length(rs) == 0) {
      rs = 1 - tab$subs_u[ tab$set2 == set.levels[i] & tab$set1 == set.levels[j]]
      if(length(rs) == 0) {
        rs = 1
      }
    }
    subs.mat[i,j] = rs[1]
  }
  cat("finished", i, "out of", d, '\n')
}
colnames(subs.mat) = set.levels
rownames(subs.mat) = set.levels
hc <- hclust(as.dist(subs.mat), method = "single")
fixed.set = hc$labels[hc$order]


## Figure 2 ====
tab1 = tab
tab1$set1 = tab$set2
tab1$set2 = tab$set1
tab1 = rbind(tab, tab1)

tritanopia_colors <- dichromat(c("white", "#D2FFFE", "#62c4f5"), type = "tritan")
combidx = combinat::combn(fixed.set, 2)
upper_tri = rbind(t(combidx), cbind(fixed.set, fixed.set)) %>%
  data.frame() %>%
  `colnames<-`(c("set1", "set2")) %>%
  arrange(factor(set1, levels = fixed.set), 
          factor(set2, levels = fixed.set)) %>% 
  left_join(tab1) %>%
  mutate(set1 = fct_inorder( set1 ), set2 = fct_inorder( set2 ),
         text = ifelse(is.na(nabla_value), "", round(nabla_value, 3)), 
         value = nabla_value) 

lower_tri = upper_tri %>%
  mutate(set.temp = set1, set1 = set2, set2 = set.temp) %>%
  select(-set.temp) %>%
  arrange(factor(set1, levels = fixed.set), 
          factor(set2, levels = fixed.set)) %>%
  mutate(text = ifelse(is.na(subs_u), "", round(subs_u, 3) ),
         value = subs_u)

all_tri = rbind(upper_tri, lower_tri) 

p1 = all_tri %>%
  ggplot(aes(x = set1, y = set2, fill = value)) +
  geom_tile(color = "white") +  # white border
  geom_tile(data = all_tri[all_tri$set1 == all_tri$set2, ], color = "lightgrey", fill = "#525354") + # diagonal tile
  geom_tile(data = all_tri[all_tri$nabla_value >= 0.5 & all_tri$subs_u > 0.8, ] %>% drop_na(), color = "red", lwd = 0.5) + # red border
  geom_text(aes(label = text), color = "black", size = 2.5) +  # Add values
  scale_fill_gradientn(colors = tritanopia_colors, na.value = "white") +
  annotate("text", x = c("2|5|8", "18|19|20|21"), y = c("18|19|20|21", "2|5|8"), label = "Stable together", size = 3) +
  labs(fill = "Upper tri: Feature perturbation metric\nLower tri: Substitutability metric", x = "Subsets of interest", y = "\n\n\n") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),  # Rotate x-axis labels
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9, angle = 90), legend.position = "right",
        axis.title = element_text(size = 9)
  )
p1 


## Figure 3 ====
all_tri_filter = tab %>%
  filter(subs_u > 0.8, nabla_value > 0.5)

opar <- par() 
par(mar = c(1.5, 1, 1, 1))
par(mfrow = c(2,3))
for(i in 1:nrow(all_tri_filter)) {
  Tau = matrix(0, nrow = 3, ncol = 5)
  Tau[1,] = matrix(1, 1, 5) 
  Tau[2,] = matrix(0, 1, 5)
  
  row_focus = all_tri_filter[i,]
  S1 = as.numeric(strsplit(row_focus$set1, "\\|")[[1]])
  S2 = as.numeric(strsplit(row_focus$set2, "\\|")[[1]])
  
  PS1 = rje::powerSet(S1)
  PS2 = rje::powerSet(S2)
  
  Tau.pair1 = c()
  Name.pair1 = list()
  for(j in PS1[2:(length(PS1)-1)]) {
    rs = subs_u_wrt_S_radar(X, y, j, S1, S2, Selection_set)
    Tau.pair1 = c(Tau.pair1, rs)
    Name.pair1 = c(Name.pair1, list(j))
  }
  Tau.pair2 = c()
  Name.pair2 = list()
  for(j in PS2[2:(length(PS2)-1)]) {
    rs = subs_u_wrt_S_radar(X, y, j, S2, S1, Selection_set)
    Tau.pair2 = c(Tau.pair2, rs)
    Name.pair2 = c(Name.pair2, list(j))
  }
  idx1 = tensr:::topK(Tau.pair1, 2)
  
  idx2 = tensr:::topK(Tau.pair2, 2)
  Name.pair2[idx2]
  Tau.pair2[idx2]
  
  Tau[3,] = c(Tau.pair1[idx1], Tau.pair2[idx2], max(Tau.pair1, Tau.pair2))
  Tau = data.frame(Tau)
  rownames(Tau) = c("1", "2", paste0("S1={", row_focus$set1, ", S2=", row_focus$set2, "}" ))
  colnames(Tau) = c(
    sapply(1:length(idx1), function(k) {
      paste0("[{", paste0(Name.pair1[idx1][[k]], collapse = ","), "}, S1]") 
    }),
    sapply(1:length(idx2), function(k) {
      paste0("[{", paste0(Name.pair2[idx2][[k]], collapse = ","), "}, S2]") 
    }),
    "Triangle"
  )
  fmsb::radarchart(
    Tau,
    pfcol = c("#99999980",NA),
    pcol= c(NA,2), plty = 1, plwd = 2,
    axislabcol = "blue", 
    caxislabels = c(0.2, 0.4, 0.6, 0.8, 1),
    vlabels = colnames(Tau),
    vlcex = 1, calcex = 0.9,
    axistype = 1,
  )
  text = paste0("S1 = {", paste0(S1, collapse = ","), "}, S2 = {", paste0(S2, collapse = ","), "}")
  mtext(text, side=1, cex = 0.7)
  
}
