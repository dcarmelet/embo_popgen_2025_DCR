sims<-read.csv("/home/carmelet/embo_popgen_2025/Me/Day 3/mosquito-task2.csv")
obs<-read.csv("/home/carmelet/embo_popgen_2025/Me/Day 3/mosquito-observed.csv")
library(abc)
ABC<-abc(obs[1,],sims[,c(1:4)],sims[,c(5:12)],tol=1e-1,method="rejection")
ABC$
