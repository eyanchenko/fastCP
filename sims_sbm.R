library(dplyr)
library(ggplot2)
library(Rcpp)
library(parallel)
library(ggpubr)
library(igraph)
library(reticulate)
library(RColorBrewer)

#py_install("cpnet")

source("~/Documents/Research/Srijan/Greedy_algorithm_R1/Code/sims_functions.R")
source_python('~/Documents/Research/Srijan/Greedy_algorithm_R1/Code/cp_be.py')

## Compare algorithm with that of CPnet and naive implementation



### Increasing p12

n.iters=100
n = 1000

p22 = 0.001
p12.seq = seq(0.002, 0.02, 0.002)


df <- tibble(iter = rep(rep(1:n.iters,each=3), length(p12.seq)), 
             method=rep(c("Algorithm 1", "Naive", "cpnet"), n.iters*length(p12.seq)), 
             p12=0, class=0, time=0, ratio=0)

cnt=1

for(p12 in p12.seq){
  p11 = 2*p12
  
  k = 0.1*n
  Cstar = c(rep(1,k), rep(0, n-k)) 
  P <- generateP(n, p11, p12, p22, prop=k/n)
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- p12
    
    A <- generateA(n, p11, p12, p22, prop=k/n)
    
    # Algorithm 1
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    df[cnt, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # Naive implementation
    start = proc.time()[3]
    C <- borgatti_naiveCpp(A)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start
    df[cnt+1, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # CPnet
    start = proc.time()[3]
    C <- C_be(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start 
    df[cnt+2, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    cnt = cnt + 3
    save(df, file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_sbm_p12_compare_103025.RData")
  }
  
  print(p12)
}


load("~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_sbm_p12_compare_103025.RData")

df_plot <- df %>% group_by(p12, method) %>% summarise(accuracy = mean(class), time = mean(time))

df_plot$method <- factor(df_plot$method, levels = c("Algorithm 1", "Naive", "cpnet"))
df_plot <- df_plot[1:20,]


p1 <- ggplot(df_plot, aes(x=p12, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  xlab(expression(p[12]))+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,3)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()


p2 <- ggplot(df_plot, aes(x=p12, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  xlab(expression(p[12]))+
  ylab("Time (sec)")+
  scale_y_continuous(trans='log10')+
  scale_colour_manual(name = "Method",values = myColors[c(1,3)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

ggsave(file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Figures/sims_sbm_p12_103025.pdf", device="pdf", width=6, height=4, units="in")



### Increasing n

n.iters=100
n.seq=seq(250, 2000, 250)

p22 = 0.001
p12 = 0.0075
p11 = 2*p12

df <- tibble(iter = rep(rep(1:n.iters,each=3), length(n.seq)), 
             method=rep(c("Algorithm 1", "Naive", "cpnet"), n.iters*length(n.seq)), 
             n=0, class=0, time=0)

cnt=1

for(n in n.seq){
  
  k = 0.1*n
  Cstar = c(rep(1,k), rep(0, n-k)) 
  P <- generateP(n, p11, p12, p22, prop=k/n)
  rho <- rho.fun(P, Cstar)
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    A <- generateA(n, p11, p12, p22, prop=k/n)
    
    # Algorithm 1
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    df[cnt, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # Naive implementation
    start = proc.time()[3]
    C <- borgatti_naiveCpp(A)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start
    df[cnt+1, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    # CPnet
    start = proc.time()[3]
    C <- C_be(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start 
    df[cnt+2, 6] <- obj.fun(A, C) / obj.fun(A, Cstar)
    
    cnt = cnt + 3
    save(df, file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_sbm_n_103025.RData")
  }
  
  print(n)
}

load("~/Documents/Research/Srijan/Greedy_algorithm_R1/Results/sims_sbm_n_103025.RData")

df_plot <- df %>% group_by(n, method) %>% summarise(accuracy = mean(class), time = mean(time))

df_plot$method <- factor(df_plot$method, levels = c("Algorithm 1", "Naive", "cpnet"))

df_plot <- df_plot[3:16,]

p1 <- ggplot(df_plot, aes(x=n, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,3)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()


p2 <- ggplot(df_plot, aes(x=n, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_y_continuous(trans='log10')+
  scale_colour_manual(name = "Method",values = myColors[c(1,3)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

ggsave(file = "~/Documents/Research/Srijan/Greedy_algorithm_R1/Figures/sims_sbm_n_103025.pdf", device="pdf", width=6, height=4, units="in")












