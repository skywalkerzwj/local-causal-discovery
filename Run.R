source("function.R")
source("compareResult.R")


data_lnc<-prepareData("all",scaling=T,lnc=T)
data<-prepareData("all",scaling=T,lnc=F)

local<-local_discovery(data,target="hsa-mir-200a",alpha=0.01,method = "mmpc")
local_lnc<-local_discovery(data_lnc,target="hsa-mir-200a",alpha=0.01,method = "mmpc")
result <- discovery_twice(data,target="hsa-mir-200a",alpha=0.01)
result_lnc <- discovery_twice(data_lnc,target="hsa-mir-200a",alpha=0.01)


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