source("function.R")
source("compareResult.R")
# library(foreach)
# library(doParallel)


data_lnc<-prepareData("all",scaling=T,lnc=T)
data<-prepareData("all",scaling=T,lnc=F)


EMT <- readEMT("EMT-35-translated.csv")
EMT <- prepareEMT(EMT,T)

mmpc_local<-local_discovery(data,target="hsa-mir-200a",alpha=0.01,method = "mmpc")
mmpc_local_lnc<-local_discovery(data_lnc,target="hsa-mir-200a",alpha=0.01,method = "mmpc")
iamb_local<-local_discovery(data,target="hsa-mir-200a",alpha=0.01,method = "iamb")
iamb_local_lnc<-local_discovery(data_lnc,target="hsa-mir-200a",alpha=0.01,method = "iamb")


result <- discovery_twice(data,target="hsa-mir-200a",alpha=0.01)
result_lnc <- discovery_twice(data_lnc,target="hsa-mir-200a",alpha=0.01)


cl<-makeCluster(2)
registerDoParallel(cl)
time<- proc.time()
twice<-discovery_twice2(data_lnc,target="hsa-mir-200c",method="mmpc",alpha=0.01)
proc.time()-time
stopCluster(cl)


unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# list[a,b,c]<-local_discovery(alpha=0.01,condition="all")
# 
# all.001<-local_discovery(alpha=0.01,condition="all")
# all.005<-local_discovery(alpha=0.05,condition="all")
# 
# cancer.005<-local_discovery(alpha=0.05,condition="Cancer")
# cancer.001<-local_discovery(alpha=0.01,condition="Cancer")
# 
# normal.005<-local_discovery(alpha=0.05,condition="Normal")
# normal.001<-local_discovery(alpha=0.01,condition="Normal")