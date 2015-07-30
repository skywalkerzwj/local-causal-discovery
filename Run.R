source("function.R")

list[a,b,c]<-local_discovery(alpha=0.01,condition="all")

all.001<-local_discovery(alpha=0.01,condition="all")
all.005<-local_discovery(alpha=0.05,condition="all")

cancer.005<-local_discovery(alpha=0.05,condition="Cancer")
cancer.001<-local_discovery(alpha=0.01,condition="Cancer")

normal.005<-local_discovery(alpha=0.05,condition="Normal")
normal.001<-local_discovery(alpha=0.01,condition="Normal")