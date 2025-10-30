### Real-data analysis
library(dplyr)
library(ggplot2)
library(Rcpp)
library(igraph)
library(igraphdata)
library(reticulate)

source("~/Documents/Research/Srijan/Greedy_algorithm/Code/sims_functions.R")
source_python('~/Documents/Research/Srijan/Greedy_algorithm/Code/cp_be.py')

########### UK faculty
# Remove time and symmetrize
load("~/Documents/Research/Srijan/CP_clean/Data/ukfaculty_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))


system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)


########### Facebook
# http://snap.stanford.edu/data/ego-Facebook.html
load("~/Documents/Research/Srijan/CP_clean/Data/FB_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)

########### Email 4
# Remove time and symmetrize
load("~/Documents/Research/Srijan/CP_clean/Data/email4_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)


########### British MP
load("~/Documents/Research/Srijan/CP_clean/Data/BritishMP_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)


########### Twitter Congress
# https://snap.stanford.edu/data/congress-twitter.html
load("~/Documents/Research/Srijan/CP_clean/Data/congress_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)


########### Political blogs
load("~/Documents/Research/Srijan/CP_clean/Data/polblogs_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)

########### DBLP
load("~/Documents/Research/Srijan/CP_clean/Data/DBLP_adj.RData")
nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)

########### Hospital
# Remove time and symmetrize
load("~/Documents/Research/Srijan/CP_clean/Data/hosp_adj.RData")
nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)


########### School, https://networkrepository.com/primary-school-proximity.php
load("~/Documents/Research/Srijan/CP_clean/Data/school_adj.RData")
nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)

########### Copenhagen BT
# Remove time and symmetrize
load("~/Documents/Research/Srijan/CP_clean/Data/copenBT_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)

########### Biological 1, https://networkrepository.com/bio-CE-GN.php
load("~/Documents/Research/Srijan/CP_clean/Data/bio1_adj.RData")
nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)

########### Biological 2, https://networkrepository.com/bio-DR-CX.php
load("~/Documents/Research/Srijan/CP_clean/Data/bio2_adj.RData")
nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)

########### Biological 3, https://networkrepository.com/bio-DM-CX.php
load("~/Documents/Research/Srijan/CP_clean/Data/bio3_adj.RData")
nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))
system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)



########### US airports
load("~/Documents/Research/Srijan/CP_clean/Data/airports_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

system.time(C <- borgattiCpp(A))
obj.fun(A, C)
sum(C)


system.time(Cnet <- cp_be(A))
Cnet = unlist(Cnet)
obj.fun(A, Cnet)
sum(Cnet)








