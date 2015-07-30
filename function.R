readData<-function(condition,target){
  translation<-read.table("BRCA/coding-11687-translation.csv",as.is=T)
  mRNA<-read.table(paste("BRCA/coding-11687-",condition,".txt",sep=""),header=T,as.is=T,sep="")
  rownames(mRNA)<-mRNA[,1]
  mRNA<-mRNA[,2:ncol(mRNA)]
  mRNA<-t(mRNA)
  colnames(mRNA)<-t(translation) #translate EntrezID to GeneSymbol
  
  miRNA<-read.table(paste("BRCA/mir-1046-",condition,".txt",sep=""),header=T,as.is=T,sep="")
  rownames(miRNA)<-miRNA[,1]
  miRNA<-miRNA[,2:ncol(miRNA)]
  miRNA<-t(miRNA)
  
  
  temp<-mRNA[,colSums(is.na(mRNA))==0] #remove columns with na and NaN
  mir.200a<-miRNA[,"hsa-mir-200a"]
  combine<-cbind(mir.200a,temp)
  combine.df<-data.frame(combine) #combine miRNA with the expression profile of mRNAs
  
  return (combine.df)
}

local_discovery<-function(alpha=0.05,scale=T,target="hsa-mir-200a",condition,test=NULL){
  if (condition == "all"){
    print("all")
    data1<-readData("Cancer",target)
    data2<-readData("Normal",target)
    #only retain those columns without NA in two conditions
    name1<-colnames(data1)
    name2<-colnames(data2)
    intesec<-intersect(name1,name2)
    data1<-data1[,intesec]
    data2<-data2[,intesec]
    data<-rbind(data1,data2)
  }
  else{
    print("else")
    data<-readData(condition,target)
    
  }
  if(scale){
    data.scaled<-scale(data,scale=F)
    data<-data.frame(data.scaled)
  }
  print("running mmpc")
  mmpc<-bnlearn:::nbr.backend(data,"mir.200a",method='mmpc',alpha=alpha,test=test)
  print("running hiton")
  hiton<-bnlearn:::nbr.backend(data,"mir.200a",method='si.hiton.pc',alpha=alpha,test=test)
#   print("running gs")
#   gs<-bnlearn:::mb.backend(data,"mir.200a",method='gs',alpha=alpha)
#   print("running iamb")
#   iamb<-bnlearn:::mb.backend(data,"mir.200a",method='iamb',alpha=alpha)
  
  result<-list("mmpc"=mmpc,"hiton"=hiton)
  return (result)
}