readData<-function(condition,lnc = F){
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
  
 
  
  temp.mRNA<-mRNA[,colSums(is.na(mRNA))==0] #remove columns with na and NaN
  temp.miRNA <- miRNA[,colSums(is.na(miRNA))==0]
  if (lnc){
    lncRNA<-read.table(paste("BRCA/noncoding-1468-",condition,".txt",sep=""),header=T,as.is=T,sep="")
    rownames(lncRNA)<-lncRNA[,1]
    lncRNA<-lncRNA[,2:ncol(lncRNA)]
    lncRNA<-t(lncRNA)
    temp.lncRNA <- lncRNA[,colSums(is.na(lncRNA))==0]
    combine<-cbind(temp.miRNA,temp.lncRNA,temp.mRNA)
  }
  else{
    combine<-cbind(temp.miRNA,temp.mRNA)
  }
  # print(colnames(lncRNA))

#   combine.df<-data.frame(combine,check.names=F) #combine miRNA with the expression profile of mRNAs
  
  return (combine)
}


prepareData<-function(condition,scaling=T,lnc=F){
  if (condition == "all"){
    data1<-readData("Cancer",lnc)
    data2<-readData("Normal",lnc)
    #only retain those columns without NA in two conditions
    name1<-colnames(data1)
    name2<-colnames(data2)
    intesec<-intersect(name1,name2)
    data1<-data1[,intesec]
    data2<-data2[,intesec]
    data<-rbind(data1,data2)
  }
  else{
    data<-readData(condition,lnc)
    
  }
  data.scaled<-scale(data,scale=scaling)
  data<-data.frame(data.scaled,check.names=F)

  return(data)
}


readEMT<-function(file){
  data<-read.table(paste("EMT/",file,sep=""),as.is = T, header = T,sep = ",",check.names=F)
  return(data)
}


prepareEMT<-function(data,scale=T){
  data.scaled <- scale(data,scale=scale)
  data <- data.frame(data.scaled,check.names=F)
}


local_discovery<-function(data,target = "hsa-mir-200a",alpha = 0.01,test = NULL,
                          method = c("mmpc","si.hiton.pc","gs","iamb","fast.iamb"),
                          correction = TRUE){
  #correction only useful for mmpc or hiton-pc
  # print("running local discovery")
  if (method=="mmpc"|method=="si.hiton.pc"){
    result<-bnlearn:::nbr.backend(data,target,method=method,alpha=alpha,test=test)
    result.temp<-result
    if (correction == T){
      result2<-list()
      for (i in 1:length(result)){
      # foreach (i=1:length(result)) %dopar%{
        new.target <- result[i]
        temp <- bnlearn:::nbr.backend(data,target=new.target,alpha=alpha,test=test,
                                method=method)
        if (!is.element(target,temp)){
          # remove new.target from result
          # print(paste(new.target," removed from result",sep=""))
          result.temp <- result.temp[which(result.temp!=new.target)]
          
        }
#         else{
#           # print(paste(new.target," saved into result",sep=""))
#           result2[[new.target]]<-temp
#         }
      }
      # result2[[target]] <- result.temp
      # print(result2)
      result<-result.temp
    }
  }
  else{
    result<-bnlearn:::mb.backend(data,target,method=method,alpha=alpha,test=test,debug=F)
  }
  print(paste("the pc set of ",target,"is:"))
  print(result)
  return (result)
}

local_discovery2<-function(data,target = "hsa-mir-200a",alpha = 0.01,test = NULL,
                         method = c("mmpc","si.hiton.pc","gs","iamb","fast.iamb"),
                         correction = TRUE){
  # Parameter correction is only useful for mmpc or hiton-pc
  # print("running local discovery")
  if (method=="mmpc"|method=="si.hiton.pc"){
    result<-bnlearn:::nbr.backend(data,target,method=method,alpha=alpha,test=test)
    result.temp<-result
    print("first mmpc finished")
    if (correction == T){
      result2<-list()
      # for (i in 1:length(result)){
      abc<-foreach (i=1:length(result),.combine = c) %dopar%{
        new.target <- result[i]
        temp <- bnlearn:::nbr.backend(data,target=new.target,alpha=alpha,test=test,
                                      method=method)
        print("second mmpc finished")
        if (!is.element(target,temp)){
          # remove new.target from result
          # print(paste(new.target," removed from result",sep=""))
          # result.temp <- result.temp[which(result.temp!=new.target)]
          new.target
        }
      }
      # print(result2)
      result<-setdiff(result,abc)
    }
  }
  else{
    result<-bnlearn:::mb.backend(data,target,method=method,alpha=alpha,test=test,debug=F)
  }
  return (result)
}


discovery_twice2<-function(data,target="hsa-mir-200a",alpha=0.01,test=NULL,method="mmpc",
                           correction=T){
  # print("finding the parent and children of target variable")
  result<-local_discovery2(data,target=target,alpha=alpha,test=test,
                          method=method,correction=T)
  # print(result)
  result2<-list()
  result2[[target]] <- result
  for (i in 1:length(result)){
    # print("finding the parent/children of parent/children of target variable")
    new.target <- result[i]
    print(paste("beginning finding the pc set of ",new.target))
    result2[[new.target]] <- local_discovery2(data,target=new.target,alpha=alpha,test=test,method=method,correction =T)
  }
  return (result2)
}



discovery_twice<-function(data,target="hsa-mir-200a",alpha=0.01,test=NULL,method="mmpc"){
  print("finding the parent and children of target variable")
  result<-local_discovery(data,target=target,alpha=alpha,scale=scale,test=test,method=method)
  print(result)
  result2<-list()
  for (i in 1:length(result)){
    # print("finding the parent/children of parent/children of target variable")
    new.target <- result[i]
    print(new.target)
    result2[[new.target]] <- local_discovery(data,target=new.target,alpha=alpha,scale=scale,test=test,method=method)
  }
  return (result2)
}

