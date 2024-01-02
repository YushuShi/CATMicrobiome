MiRKATR2<-function(metaData,outcomeVar,Ks){ 
  temp<-NA
  if(length(outcomeVar)>1){
    temp<-MiRKATS(obstime = metaData[,outcomeVar[1]], delta = metaData[,outcomeVar[2]], Ks = Ks, perm=TRUE,omnibus="permutation", nperm=1,returnR2 = TRUE)
  }else{
    outcome<-NULL
    if(is.vector(metaData)){
      outcome<-metaData
    }else{
      outcome<-metaData[,outcomeVar]
    }
    if(length(table(outcome))==2){
      temp<-MiRKAT(outcome, Ks = Ks, method="permutation",
                   out_type="D",omnibus="permutation", nperm=1,returnR2 = TRUE)   
    }else{
      temp<-MiRKAT(outcome, Ks = Ks,method="permutation",
                   out_type="C",omnibus="permutation", nperm=1,returnR2 = TRUE)  
    } 
  }
  max(temp$R2)
}