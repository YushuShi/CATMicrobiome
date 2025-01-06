MiRKATR2<-function(metaData,outcomeVar,adjVar,Ks){ 
  temp<-NA
  if(length(outcomeVar)>1){
    if(is.null(adjVar)){
    temp<-MiRKATS(obstime = metaData[,outcomeVar[1]], delta = metaData[,outcomeVar[2]], Ks = Ks, perm=TRUE,omnibus="permutation", nperm=1,returnR2 = TRUE)
    }else{
      X<-as.matrix(metaData[,adjVar],ncol=length(adjVar))
      temp<-MiRKATS(obstime = metaData[,outcomeVar[1]], delta = metaData[,outcomeVar[2]], X=X, Ks = Ks, perm=TRUE,omnibus="permutation", nperm=length(adjVar)+1,returnR2 = TRUE)
    }
  }else{
    outcome<-NULL
    if(is.vector(metaData)){
      outcome<-metaData
    }else{
      outcome<-metaData[,outcomeVar]
    }
    if(length(table(outcome))==2){
      if(is.null(adjVar)){
      temp<-MiRKAT(outcome, Ks = Ks, method="permutation",
                   out_type="D",omnibus="permutation", nperm=1,returnR2 = TRUE)
      }else{
        X<-matrix(metaData[,adjVar],ncol=length(adjVar))
        temp<-MiRKAT(outcome, X=X, Ks = Ks,method="permutation",
                     out_type="D",omnibus="permutation", nperm=length(adjVar)+1,returnR2 = TRUE)
      }
    }else{
      if(is.null(adjVar)){
      temp<-MiRKAT(outcome, Ks = Ks,method="permutation",
                   out_type="C",omnibus="permutation", nperm=1,returnR2 = TRUE)  
      }else{
        X<-matrix(metaData[,adjVar],ncol=length(adjVar))
        temp<-MiRKAT(outcome, X=X, Ks = Ks,method="permutation",
                     out_type="C",omnibus="permutation", nperm=length(adjVar)+1,returnR2 = TRUE)
      }
    } 
  }
  max(temp$R2)
}