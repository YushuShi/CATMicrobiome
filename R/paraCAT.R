paraCAT<-function(testIter,testList,otutable,taxonomy,metric,metaData,outcomeVar,numBS,origR2,origBS,tree){
  otutemp<-otutable
  otuUnder<-rownames(otutable)[apply(taxonomy, 1, function(x) sum(any(x %in% testList[testIter])))>0] 
  otutemp[otuUnder,]<-0
  distMat<-compDist(otutemp,metric,tree)
  suppressMessages(distResult<-adonis(distMat~metaData[,outcomeVar],permutations = 1)$aov.tab)
  taxaR2<-distResult$R2[1]
  taxaBS<-rep(NA,numBS)
  for(BSiter in 1:numBS){
    set.seed(BSiter)
    indi<-sample(1:nrow(metaData),nrow(metaData),replace = TRUE)
    outcomeBS<-metaData[indi,outcomeVar]
    distMatBS<-as.matrix(distMat)
    distMatMatrix<-as.matrix(distMat)
    for(i in 1:length(outcomeBS)){
      for(j in 1:length(outcomeBS)){
        distMatBS[i,j]<-distMatMatrix[indi[i],indi[j]]
      }
    }
    suppressMessages(resultBS<-adonis(distMatBS~outcomeBS,permutations = 1)$aov.tab)
    taxaBS[BSiter]<-resultBS$R2[1]
  }
  BSpvalue<-empiricalP(taxaBS-origBS) # supposed to be something less than 0
  c(origR2,taxaR2,BSpvalue)
}
