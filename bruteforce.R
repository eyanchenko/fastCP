library(dplyr)
library(Rcpp)
library(ggplot2)

source("sims_functions.R")

n = 20
Cbrute = expand.grid(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1),
                     c(0,1), c(0,1), c(0,1), c(0,1), c(0,1),
                     c(0,1), c(0,1), c(0,1), c(0,1), c(0,1),
                     c(0,1), c(0,1), c(0,1), c(0,1), c(0,1))
Cbrute = Cbrute[2:(2^n-1),]




niters=100
p.seq = seq(0.05, 0.95, 0.05)

df <- tibble(iter = rep(1:niters, length(p.seq)), p = 0.0,  prop=0.0)

cnt = 1
for(p in p.seq){
  for(i in 1:niters){
    
    df[cnt, 2] = p
    
    A <- generateA(n, p, p, p)
    
    apply.fun <- function(j){
      obj.fun(A, as.numeric(Cbrute[j, ]) )
    }
    
    df[cnt, 3] = obj.fun(A, borgattiCpp(A)) / max(unlist(parallel::mclapply(1:nrow(Cbrute), apply.fun, mc.cores = 8)))
    cnt = cnt +1
    print(i)
    save(df, file="~/Documents/Research/Srijan/Greedy_algorithm/Results/df_brute_force.RData")
  }
  print(p)
}

load(file="~/Documents/Research/Srijan/Greedy_algorithm/Results/df_brute_force.RData")

df_plot <- df %>% group_by(p) %>% summarize(prop = mean(prop))
df_plot

ggplot(df_plot, aes(x=p, y=prop))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  geom_hline(yintercept = 0.90, color="red")+
  ylab("Proportion of Global Maximum")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggsave(file = "~/Documents/Research/Srijan/Greedy_algorithm/Figures/brute.pdf", device="pdf", width=4, height=4, units="in")

