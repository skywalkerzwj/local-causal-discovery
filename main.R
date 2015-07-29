library(parallel)
makeCluster(3)

#load EMT data
data<-read.table("EMT/EMT-35-translated.csv",header=T,as.is=T,sep=",")

#load BRCA data
mRNA<-read.table("BRCA/coding-11687-Normal.txt",header=T,as.is=T,sep="")
rownames(mRNA)<-mRNA[,1]
mRNA<-mRNA[,2:ncol(mRNA)]
mRNA<-t(mRNA)

miRNA<-read.table("BRCA/mir-1046-Normal.txt",header=T,as.is=T,sep="")
rownames(miRNA)<-miRNA[,1]
miRNA<-miRNA[,2:ncol(miRNA)]
miRNA<-t(miRNA)

lncRNA<-read.table("BRCA/noncoding-1468-Normal.txt",header=T,as.is=T,sep="")
rownames(lncRNA)<-lncRNA[,1]
lncRNA<-lncRNA[,2:ncol(lncRNA)]
lncRNA<-t(lncRNA)


#find the Markov Blanket of a particular miRNA: hsa-mir-200a
temp<-mRNA[,colSums(is.na(mRNA))==0] #remove columns with na and NaN
mir.200a<-miRNA[,"hsa-mir-200a"]
combine<-cbind(mir.200a,temp)
combine.df<-data.frame(combine)
temp<-bnlearn:::mb.backend(combine.df[,1:500],"mir.200a",method='iamb')

a<-proc.time()
temp3<-bnlearn:::mb.backend(combine.df,"mir.200a",method='fast.iamb')
proc.time()-a

