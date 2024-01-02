paraCATMiRKAT<-function(testIter,testList,otutable,taxonomy,metric,metaData,outcomeVar,numBS,origR2,origBS,tree){
  otutemp<-otutable
  otuUnder<-rownames(otutable)[apply(taxonomy, 1, function(x) sum(any(x %in% testList[testIter])))>0] 
  otutemp[otuUnder,]<-0
  Ks<-compDistList(otutemp,metric,tree)
  taxaR2<-MiRKATR2(metaData,outcomeVar,Ks)
  taxaBS<-rep(NA,numBS)
  for(BSiter in 1:numBS){
    set.seed(BSiter)
    indi<-sample(1:nrow(metaData),nrow(metaData),replace = TRUE)
    outcomeBS<-metaData[indi,outcomeVar]
    KsBS<-Ks
    for(k in 1:length(Ks)){
      for(i in 1:length(outcomeBS)){
        for(j in 1:length(outcomeBS)){
          KsBS[[k]][i,j]<-Ks[[k]][indi[i],indi[j]]
        }
      }
    }
    taxaBS[BSiter]<-MiRKATR2(outcomeBS,outcomeVar,KsBS)
  }
  BSpvalue<-empiricalP(taxaBS-origBS) # supposed to be something less than 0
  c(origR2,taxaR2,BSpvalue)
}
