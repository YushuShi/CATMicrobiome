paraCATMiRKAT<-function(testTarget,otutable,taxonomy,metric,metaData,outcomeVar,numPerm,origR2,adjVar,tree){
  taxaR2<-rep(NA,numPerm)
  for(permIter in 1:numPerm){
    set.seed(permIter)
    otutable2<-otutable
    temp<-rownames(otutable2)[apply(taxonomy, 1, function(x) sum(any(x %in% testTarget)))>0]
    otutable2[temp,]<-otutable2[temp,sample(1:ncol(otutable2),ncol(otutable2),replace = FALSE)]
    Ks<-compDistList(otutable2,metric,tree)
    taxaR2[permIter]<-MiRKATR2(metaData,outcomeVar,adjVar,Ks)
  }
  sum(origR2<taxaR2)/numPerm
}
