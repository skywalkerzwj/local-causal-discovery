library(parallel)
makeCluster(3)

#load EMT data
data<-read.table("EMT/EMT-35-translated.csv",header=T,as.is=T,sep=",")

#load BRCA data
# condition="Normal.txt"
condition="Cancer.txt"
translation<-read.table("BRCA/coding-11687-translation.csv",as.is=T)
mRNA<-read.table(paste("BRCA/coding-11687-",condition,sep=""),header=T,as.is=T,sep="")
rownames(mRNA)<-mRNA[,1]
mRNA<-mRNA[,2:ncol(mRNA)]
mRNA<-t(mRNA)
colnames(mRNA)<-t(translation) #translate EntrezID to GeneSymbol

miRNA<-read.table(paste("BRCA/mir-1046-",condition,sep=""),header=T,as.is=T,sep="")
rownames(miRNA)<-miRNA[,1]
miRNA<-miRNA[,2:ncol(miRNA)]
miRNA<-t(miRNA)

# lncRNA<-read.table("BRCA/noncoding-1468-Normal.txt",header=T,as.is=T,sep="")
# rownames(lncRNA)<-lncRNA[,1]
# lncRNA<-lncRNA[,2:ncol(lncRNA)]
# lncRNA<-t(lncRNA)


#find the Markov Blanket of a particular miRNA: hsa-mir-200a
temp<-mRNA[,colSums(is.na(mRNA))==0] #remove columns with na and NaN
mir.200a<-miRNA[,"hsa-mir-200a"]
combine<-cbind(mir.200a,temp)
combine.df<-data.frame(combine) #combine miRNA with the expression profile of mRNAs


mmpc_cancer<-bnlearn:::nbr.backend(combine.df,"mir.200a",method='mmpc',alpha=0.01)
hiton_cancer<-bnlearn:::nbr.backend(combine.df,"mir.200a",method='si.hiton.pc',alpha=0.01)
gs_cancer<-bnlearn:::mb.backend(combine.df,"mir.200a",method="gs",alpha=0.01)
iamb<-bnlearn:::mb.backend(combine.df,"mir.200a",method="gs",alpha=0.01)

