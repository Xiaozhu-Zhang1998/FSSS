library(tidyverse)
library(latex2exp)
library(patchwork)


# Figure 8 ========
alpha0 = 0.1
eta = 0.01
sstar = 10
S0 = seq(from = 10, to = 20, by = 1)
p = 2000
RSc = matrix(0, nrow = 0, ncol = 6)
for(s0 in S0) {
  for(D in c(2,3,4)) {
    for(p in c(2000, 4000, 6000)) {
      # for individual gamma
      gamma1 = (s0 - sstar) / (p - sstar)
      # for paired gamma
      gamma2 = (s0 - sstar) * (s0 - sstar - 1) / ( (p - sstar) * (p - sstar - 1) ) +
        (s0 - sstar) * (p - s0) / ( (p - sstar) * (p - sstar - 1) )
      # for individual noise
      gamma3 = 1 - prod( sapply(1:D, function(d) {
        (p - s0 - (d-1)) / (p - sstar - (d-1))
      }) )
      
      gamma = max(gamma1, gamma2, gamma3)
      fac = sstar * choose(D-1,2) + p - sstar 
      prob = 1 - gamma^2 * fac / (1 - 2 * alpha0)
      ubd = prob * sstar * eta^2 / (1 + eta^2) + (1 - prob) * s0
      
      rs = c(s0, D, p, prob, gamma, ubd)
      RSc = rbind(RSc, rs)
    }
  }
}

fontsize = 7

p1 = RSc %>%
  data.frame() %>%
  `colnames<-`(c("s0", "D", "p", "prob", "gamma", "ubd")) %>%
  mutate(D = factor(D), p = factor(p)) %>%
  ggplot(aes(x = s0, y = gamma, col = D, linetype = p)) +
  geom_line(linewidth = 0.7) +
  labs(x = TeX(r"($s_0$)"),
       y = TeX(r"($\gamma$)"),
       color = "Cluster size (D)", 
       linetype = "Data size (p)" ) +
  scale_x_continuous(breaks = S0) + 
  ylim(0, 0.03) +
  theme_bw() +
  theme(
    legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
    axis.title = element_text(size=8), plot.title = element_text(size=8),
    axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize)
  ) 


p2 = RSc %>%
  data.frame() %>%
  `colnames<-`(c("s0", "D", "p", "prob", "gamma", "ubd")) %>%
  mutate(D = factor(D), p = factor(p)) %>%
  ggplot(aes(x = s0, y = prob, col = D, linetype = p)) +
  geom_line(linewidth = 0.7) +
  labs(x = TeX(r"($s_0$)"),
       y = "Probability",
       color = "Cluster size (D)", 
       linetype = "Data size (p)" ) +
  scale_x_continuous(breaks = S0) + 
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(
    legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
    axis.title = element_text(size=8), plot.title = element_text(size=8),
    axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize)
  ) 


p3 = RSc %>%
  data.frame() %>%
  `colnames<-`(c("s0", "D", "p", "prob", "gamma", "ubd")) %>%
  mutate(D = factor(D), p = factor(p)) %>%
  ggplot(aes(x = s0, y = ubd, col = D, linetype = p)) +
  geom_line(linewidth = 0.7) +
  labs(x = TeX(r"($s_0$)"),
       y = "FPE upper bound",
       color = "Cluster size (D)", 
       linetype = "Data size (p)" ) +
  scale_x_continuous(breaks = S0) + 
  geom_abline(slope = 1, intercept = 0, col = "#999999") +
  geom_abline(slope = 1, intercept = -S0[1], col = "#999999") +
  theme_bw() +
  theme(
    legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
    axis.title = element_text(size=8), plot.title = element_text(size=8),
    axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize)
  ) 

p1 + p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')


# Figure 10  ========
alpha0 = 0.1
eta = 0.01
sindep = 5
B = 500
S0 = seq(from = 8, to = 25, by = 1)

RSb = matrix(0, nrow = 0, ncol = 6)
for(s0 in S0) {
  for(B in c(500, 1000, 1500) ) {
    
    nb_b_ratio = 0.998
    NBstar = round(B * nb_b_ratio)
    SBstar = B - NBstar
    sstar = sindep + SBstar
    for(p in c(2000, 4000, 8000)) {
      if(NBstar == 0) {
        gamma = (s0 - sstar) / (p - sstar)
      } else {
        # for Xc
        gamma1 = (s0 - sstar) / (p - sstar) + 
          sum(sapply(0:NBstar, function(k) {
            val = (SBstar + k) / (B + eta^2)
            prob = extraDistr::dmvhyper(
              x = c(k, 0, s0-sstar-k), 
              n = c(NBstar, 1, p-sstar-NBstar-1), 
              k = s0 - sstar
            )
            val * prob
          }))
        # for NBstar
        gamma2 = (s0 - sstar) / (p - sstar) + 
          sum(sapply(0:(NBstar-1), function(k) {
            val = 1 / (NBstar - k + eta^2)
            prob = extraDistr::dmvhyper(
              x = c(k, 0, 1, s0-sstar-k), 
              n = c(NBstar-1, 1, 1, p-sstar-NBstar-1), 
              k = s0 - sstar
            )
            val * prob
          }))
        # for NIstar
        gamma3 = (s0 - sstar) / (p - sstar)
        gamma = max(gamma1, gamma2, gamma3)
      }
      fac = p-sstar
      prob = 1 - gamma^2 * fac / (1 - 2 * alpha0)
      ubd = prob * eta^2 / (1 + eta^2) + (1 - prob) * s0
      
      rs = c(s0, B, p, prob, gamma, ubd)
      RSb = rbind(RSb, rs)
    }
  }
}

fontsize = 7

p1 = RSb %>%
  data.frame() %>%
  `colnames<-`(c("s0", "B", "p", "prob", "gamma", "ubd")) %>%
  mutate(B = factor(B), p = factor(p)) %>%
  ggplot(aes(x = s0, y = gamma, col = B, linetype = p)) +
  geom_line(linewidth = 0.7) +
  labs(x = TeX(r"($s_0$)"),
       y = TeX(r"($\gamma$)"),
       color = "Block size (K)", 
       linetype = "Data size (p)" ) +
  scale_x_continuous(breaks = S0) +
  ylim(0, 0.03) +
  theme_bw() +
  theme(
    legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
    axis.title = element_text(size=8), plot.title = element_text(size=8),
    axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize)
  ) 


p2 = RSb %>%
  data.frame() %>%
  `colnames<-`(c("s0", "B", "p", "prob", "gamma", "ubd")) %>%
  mutate(B = factor(B), p = factor(p)) %>%
  ggplot(aes(x = s0, y = prob, col = B, linetype = p)) +
  geom_line(linewidth = 0.7) +
  labs(x = TeX(r"($s_0$)"),
       y = "Probability",
       color = "Block size (K)", 
       linetype = "Data size (p)" ) +
  scale_x_continuous(breaks = S0) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(
    legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
    axis.title = element_text(size=8), plot.title = element_text(size=8),
    axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize)
  ) 


p3 = RSb %>%
  data.frame() %>%
  `colnames<-`(c("s0", "B", "p", "prob", "gamma", "ubd")) %>%
  mutate(B = factor(B), p = factor(p)) %>%
  ggplot(aes(x = s0, y = ubd, col = B, linetype = p)) +
  geom_line(linewidth = 0.7) +
  labs(x = TeX(r"($s_0$)"),
       y = "FPE upper bound",
       color = "Block size (K)", 
       linetype = "Data size (p)" ) +
  scale_x_continuous(breaks = S0) +
  geom_abline(slope = 1, intercept = 0, col = "#999999") +
  geom_abline(slope = 1, intercept = -S0[1], col = "#999999") +
  theme_bw() +
  theme(
    legend.title=element_text(size=fontsize), legend.text = element_text(size=fontsize), 
    axis.title = element_text(size=8), plot.title = element_text(size=8),
    axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize)
  ) 

p1 + p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')



# Figure 9 ========
alpha = 0.1
eta = 0.01
sstar = 10
S0 = seq(from = 30, to = 100, by = 4)
RSc = matrix(0, nrow = 0, ncol = 4)

for(s0 in S0) {
  for(p in c(500, 1000, 2000)) {
    for(alpha in c(0.05, 0.1, 0.15)) {
      ubd = 1/(1- 2*alpha) * s0^2 / p + 
        1/(1-2*alpha) * (2 * s0 * sqrt(eta^2 / (1 + eta^2)) + 4 * p * eta^2 / (1 + eta^2) ) +
        s0 * eta^2 / (1 + eta^2)
      rs = c(s0, p, alpha, ubd)
      RSc = rbind(RSc, rs)
    }
  }
}

RSc %>%
  data.frame() %>%
  `colnames<-`(c("s0", "p", "alpha", "ubd")) %>%
  mutate(p = factor(p), alpha = factor(alpha)) %>%
  ggplot(aes(x = s0, y = ubd, col = p, linetype = alpha)) +
  geom_line(linewidth = 0.7) +
  labs(x = TeX(r"($s_0$)"),
       y = "FPE upper bound",
       color = "Data size (p)", 
       linetype = TeX(r"($\alpha$)") ) +
  scale_x_continuous(breaks = S0) + 
  # geom_abline(slope = 1, intercept = 0, col = "#999999") +
  theme_bw()

