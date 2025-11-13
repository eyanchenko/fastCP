library(dplyr)
library(ggplot2)
library(Rcpp)
library(parallel)
library(ggpubr)
library(igraph)
library(reticulate)
library(RColorBrewer)

#py_install("cpnet")

source("~/Documents/Research/Srijan/Greedy_algorithm_R1/fastCP/sims_functions.R")
source_python('~/Documents/Research/Srijan/Greedy_algorithm_R1/fastCP/cp_be.py')

### Increasing p12

n.iters=100
n = 1000

p12.seq = seq(0.05, 0.15, length=11)
p22=0.05


df <- tibble(iter = rep(rep(1:n.iters,each=3), length(p12.seq)), 
             method=rep(c("greedyFast", "greedyNaive", "cpnet"), n.iters*length(p12.seq)), 
             p12=0, class=0, time=0, ratio=0)

cnt=1

for(p12 in p12.seq){
  p11 = 2*p12
  
  k = 0.1*n
  Cstar = c(rep(1,k), rep(0, n-k)) 
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+2), 3] <- p12
    
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    
    
    # greedyFast
    start = proc.time()[3]
    C <- greedyFast(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    df[cnt, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # greedyNaive
    start = proc.time()[3]
    C <- greedyNaive(A)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start
    df[cnt+1, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # cpnet
    start = proc.time()[3]
    C <- C_be(A, 1)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start 
    df[cnt+2, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    cnt = cnt + 3
    print(iter)
    save(df, file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_dcbm_p12_111325.RData")
  }
  
  print(p12)
}


load("~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_dcbm_p12_111325.RData")

df_plot <- df %>% group_by(p12, method) %>% summarise(accuracy = mean(class), time = mean(time))
df_plot$method <- factor(df_plot$method, levels = c("greedyFast", "greedyNaive", "cpnet"))
df_plot <- df_plot[1:33,]


p1 <- ggplot(df_plot, aes(x=p12, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  xlab(expression(p[12]))+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,5,2)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()


p2 <- ggplot(df_plot, aes(x=p12, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  xlab(expression(p[12]))+
  ylab("Time (sec)")+
  scale_y_continuous(trans='log10')+
  scale_colour_manual(name = "Method",values = myColors[c(1,5,2)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

ggsave(file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Figures/sims_dcbm_p12_111325.pdf", device="pdf", width=6, height=4, units="in")



### Increasing n

n.iters=100
n.seq = seq(500, 2000, 250)

p22 = 0.05
p12 = 0.10
p11 = 2*p12

df <- tibble(iter = rep(rep(1:n.iters,each=3), length(n.seq)), 
             method=rep(c("greedyFast", "greedyNaive", "cpnet"), n.iters*length(n.seq)), 
             n=0, class=0, time=0)

cnt=1

for(n in n.seq){
  
  k = 0.1*n
  Cstar = c(rep(1,k), rep(0, n-k)) 
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+2), 3] <- n
    
    theta <- runif(n, 0.6, 0.8)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    
    # greedyFast
    start = proc.time()[3]
    C <- greedyFast(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    df[cnt, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # greedyNaive
    start = proc.time()[3]
    C <- greedyNaive(A)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start
    df[cnt+1, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # CPnet
    start = proc.time()[3]
    C <- C_be(A, 1)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start 
    df[cnt+2, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    cnt = cnt + 3
    print(iter)
    save(df, file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_dcbm_n_111325.RData")
  }
  
  print(n)
}

load("~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_dcbm_n_111325.RData")

df_plot <- df %>% group_by(n, method) %>% summarise(accuracy = mean(class), time = mean(time))
df_plot$method <- factor(df_plot$method, levels = c("greedyFast", "greedyNaive", "cpnet"))
df_plot <- df_plot[1:21,]

p1 <- ggplot(df_plot, aes(x=n, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,5,2)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()


p2 <- ggplot(df_plot, aes(x=n, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_y_continuous(trans='log10')+
  scale_colour_manual(name = "Method",values = myColors[c(1,5,2)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

ggsave(file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Figures/sims_dcbm_n_111325.pdf", device="pdf", width=6, height=4, units="in")












