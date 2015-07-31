difference<-function(list1,list2){
  # this function compare the difference of two lists of lists
  # NOTE: this is asymmetric, as setdiff.
  len1 <- length(list1)
  len2 <- length(list2)
  diff <- list()
  for (i in 1:ifelse(len1>len2,len2,len1)){
    col <- names(list1)[i]
    diff[[col]] <- setdiff(list1[[col]],list2[[col]])
  }
  return(diff)
}

toFile <- function(list,file){
  # this function output list of list to a csv, for cyctoscape draw graph
  for (i in 1:length(list)){
    for (j in 1:length(list[[i]])){
      write(paste(names(list)[i],list[[i]][j],sep=","),file=file,append=T)
    }
  }
  # unlink("data")
}