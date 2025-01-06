paraCAT<-function(testTarget,otutable,taxonomy,metric,metaData,outcomeVar,adjVar,numPerm,origR2,tree){
  taxaR2<-rep(NA,numPerm)
  for(permIter in 1:numPerm){
    set.seed(permIter)
    otutable2<-otutable
    temp<-rownames(otutable2)[apply(taxonomy, 1, function(x) sum(any(x %in% testTarget)))>0]
    otutable2[temp,]<-otutable2[temp,sample(1:ncol(otutable2),ncol(otutable2),replace = FALSE)]
    distMat2<-compDist(otutable2,metric,tree)
    if(is.null(adjVar)){
      taxaR2[permIter]<-adonis2(distMat2~metaData[,outcomeVar],permutations = 1)$R2[1]
    }else{
      newData<-metaData[,c(adjVar,outcomeVar)]
      taxaR2[permIter]<-adonis2(distMat2~.,data=newData,permutations = length(adjVar)+1)$R2[length(adjVar)+1]
    }
  }
  sum(origR2<taxaR2)/numPerm
}
