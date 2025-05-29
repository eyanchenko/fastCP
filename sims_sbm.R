library(dplyr)
library(ggplot2)
library(Rcpp)
library(parallel)
library(ggpubr)
library(igraph)
library(reticulate)
library(RColorBrewer)

#py_install("cpnet")

source("sims_functions.R")
source_python('cp_be.py')

### Increasing p12

n.iters=100
n = 1000

p22 = 0.001
p12.seq = seq(0.002, 0.02, 0.002)


df <- tibble(iter = rep(rep(1:n.iters,each=2), length(p12.seq)), 
             method=rep(c("BE", "BE_python"), n.iters*length(p12.seq)), 
             p12=0, class=0, time=0)

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
    
    # CPnet
    start = proc.time()[3]
    C <- C_be(A)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    cnt = cnt + 2
    save(df, file = "sims_sbm_p12.RData")
  }
  
  print(p12)
}


load("sims_sbm_p12.RData")

df_plot <- df %>% group_by(p12, method) %>% summarise(accuracy = mean(class), time = mean(time))

df_plot$method[df_plot$method=="BE"] <- "Algorithm 1"
df_plot$method[df_plot$method=="BE_python"] <- "cpnet"


df_plot$method <- factor(df_plot$method, levels = c("Algorithm 1", "cpnet"))
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

ggsave(file = "sims_sbm_p12.pdf", device="pdf", width=6, height=4, units="in")



### Increasing n

n.iters=100
n.seq=seq(250, 2000, 250)

p22 = 0.001
p12 = 0.0075
p11 = 2*p12

df <- tibble(iter = rep(rep(1:n.iters,each=2), length(n.seq)), 
             method=rep(c("BE", "BE_python"), n.iters*length(n.seq)), 
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
    
    # Cpnet
    start = proc.time()[3]
    C <- C_be(A)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    cnt = cnt + 2
    save(df, file = "sims_sbm_n.RData")
  }
  
  print(n)
}

load("sims_sbm_n.RData")

df_plot <- df %>% group_by(n, method) %>% summarise(accuracy = mean(class), time = mean(time))

df_plot$method[df_plot$method=="BE"] <- "Algorithm 1"
df_plot$method[df_plot$method=="BE_python"] <- "cpnet"

df_plot$method <- factor(df_plot$method, levels = c("Algorithm 1", "cpnet"))

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

ggsave(file = "sims_sbm_n.pdf", device="pdf", width=6, height=4, units="in")










