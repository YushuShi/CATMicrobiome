# source("compDist.R")
# source("MiRKATR2.R")
# source("compDistList.R")
# source("paraCAT.R")
# source("paraCATMiRKAT.R")

CAT<-function(testList,otutable,taxonomy,metric="Weighted UniFrac",metaData,outcomeVar,tree=NULL,adjVar=NULL,method="PERMANOVA",numPerm=1000, parallel=TRUE,nCore=2){
  if(is.null(tree)){
    if(sum(metric%in%c("WeightUniFrac","Unweighted UniFrac","robust"))>0){
      stop("UniFrac distance needs a tree!")
    }
    if((sum(rownames(otutable)%in% rownames(taxonomy))!=nrow(taxonomy))|
       (sum(rownames(taxonomy) %in% rownames(otutable))!=nrow(otutable))){
      if(nrow(otutable)<nrow(taxonomy)){
        otutable<-otutable[rownames(otutable)%in% rownames(taxonomy),]
        print(paste0("Before padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
        padding<-matrix(0,ncol=ncol(otutable),nrow=nrow(taxonomy)-sum(rownames(taxonomy)%in% rownames(otutable)))
        rownames(padding)<-rownames(taxonomy)[!(rownames(taxonomy)%in% rownames(otutable))]
        colnames(padding)<-colnames(taxonomy)
        otutable<-rbind(otutable,padding)
        print(paste0("After padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
      }else{
        taxonomy<-taxonomy[rownames(taxonomy)%in% rownames(otutable),]
        print(paste0("Before padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
        padding<-matrix(0,ncol=ncol(taxonomy),nrow=nrow(otutable)-sum(rownames(otutable)%in% rownames(taxonomy)))
        rownames(padding)<-rownames(otutable)[!(rownames(otutable)%in% rownames(taxonomy))]
        colnames(padding)<-colnames(otutable)
        taxonomy<-rbind(taxonomy,padding)
        print(paste0("After padding: number of rows in OTU table: ",nrow(otutable),", and number of rows in taxonomy table: ",nrow(taxonomy)))
      }
    }
  }else{
    tree$root.edge<-0
    if((sum(rownames(otutable)%in% tree$tip.label)!=length(tree$tip.label))|
       (sum(tree$tip.label %in% rownames(otutable))!=nrow(otutable))){
      otutable<-otutable[rownames(otutable)%in% tree$tip.label,]
      print(paste0(nrow(otutable)," OTUs and ",length(tree$tip.label)," tip nodes."))
      if(nrow(otutable)<length(tree$tip.label)){
        padding<-matrix(0,ncol=ncol(otutable),nrow=length(tree$tip.label)-nrow(otutable))
        rownames(padding)<-tree$tip.label[!(tree$tip.label%in% rownames(otutable))]
        colnames(padding)<-colnames(otutable)
        otutable<-rbind(otutable,padding)
      }
    }
    if((sum(rownames(taxonomy)%in% tree$tip.label)!=length(tree$tip.label))|
       (sum(tree$tip.label %in% rownames(taxonomy))!=nrow(taxonomy))){
      taxonomy<-taxonomy[rownames(taxonomy)%in% tree$tip.label,]
      print(paste0(nrow(taxonomy)," rows in taxonomy table, and ",length(tree$tip.label)," tip nodes."))
      if(nrow(taxonomy)<length(tree$tip.label)){
        padding<-matrix(0,ncol=ncol(taxonomy),nrow=length(tree$tip.label)-sum(tree$tip.label%in% rownames(taxonomy)))
        rownames(padding)<-tree$tip.label[!(tree$tip.label%in% rownames(taxonomy))]
        colnames(padding)<-colnames(taxonomy)
        taxonomy<-rbind(taxonomy,padding)
      }
    }
  }
  otutable<-otutable[,rownames(metaData)]
  taxonomy<-taxonomy[rownames(otutable),]
  pvalue<-NA
  
  if(method=="PERMANOVA"){
    distMat<-compDist(otutable,metric,tree)
    if(is.null(adjVar)){
      origR2<-adonis2(distMat~metaData[,outcomeVar],permutations = 1)$R2[1]
    }else{
      newData<-metaData[,c(adjVar,outcomeVar)]
      origR2<-adonis2(distMat~.,data=newData,permutations = length(adjVar)+1)$R2[length(adjVar)+1]
    }

    if(parallel){
      cl <- makeCluster(nCore)
      registerDoParallel(cl)
      pvalue<-foreach(testTarget=testList,.combine='c',
                          .export =c("compDist","paraCAT"),
                          .packages=c("ape","vegan","GUniFrac")
      ) %dopar% paraCAT(testTarget,otutable,taxonomy,metric,metaData,outcomeVar,adjVar,numPerm,origR2,tree)
      stopImplicitCluster()
      stopCluster(cl)
      registerDoSEQ()
      rm(cl)
      gc()
    }else{
      pvalue<-NULL
      for(testIter in 1:length(testList)){
        pvalue<-c(pvalue,paraCAT(testList[testIter],otutable,taxonomy,metric,metaData,outcomeVar,adjVar,numPerm,origR2,tree))
      }
    }
  }else{
    Ks<-compDistList(otutable,metric,tree)
    origR2<-MiRKATR2(metaData,outcomeVar,adjVar,Ks)
    if(parallel){
      cl <- makeCluster(nCore)
      registerDoParallel(cl)
      pvalue<-foreach(testTarget=testList,.combine='c',
                      .export =c("compDistList","paraCATMiRKAT","MiRKATR2"),
                          .packages=c("ape","vegan","GUniFrac","MiRKAT")
      ) %dopar% paraCATMiRKAT(testTarget,otutable,taxonomy,metric,metaData,outcomeVar,numPerm,origR2,adjVar,tree)
      stopImplicitCluster()
      stopCluster(cl)
      registerDoSEQ()
      rm(cl)
      gc()
    }else{
      pvalue<-NULL
      for(testIter in 1:length(testList)){
        pvalue<-c(pvalue,paraCATMiRKAT(testList[testIter],otutable,taxonomy,metric,metaData,outcomeVar,numPerm,origR2,adjVar,tree))
      }
    }
  }
  names(pvalue)<-testList
  pvalue
}


