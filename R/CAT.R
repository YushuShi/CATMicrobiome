# source("compDist.R")
# source("MiRKATR2.R")
# source("compDistList.R")
# source("paraCAT.R")
# source("paraCATMiRKAT.R")

empiricalP<-function(val1000){
  val1000<-val1000[!is.na(val1000)]
  p.value<-sum(val1000>0)/length(val1000)
  p.value
}


CAT<-function(testList,otutable,taxonomy,metric="Weighted UniFrac",metaData,outcomeVar,tree=NULL,method="PERMANOVA",numBS=1000, parallel=TRUE,nCore=2){
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
  testResult<-NA
  
  if(method=="PERMANOVA"){
    distMat<-compDist(otutable,metric,tree)
    suppressMessages(distResult<-adonis(distMat~metaData[,outcomeVar],permutations = 1)$aov.tab)
    origR2<-distResult$R2[1]
    origBS<-rep(NA,numBS)
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
      if(length(unique(outcomeBS)) > 1){
          suppressMessages(resultBS<-adonis(distMatBS~outcomeBS,permutations = 1)$aov.tab)
          origBS[BSiter]<-resultBS$R2[1]
        }
    }
    if(parallel){
      registerDoParallel(nCore)
      testResult<-foreach(testIter=1:length(testList),.combine='rbind',
                          .packages=c("ape","vegan","GUniFrac")
      ) %dopar% paraCAT(testIter,testList,otutable,taxonomy,metric,metaData,outcomeVar,numBS,origR2,origBS,tree)
      stopImplicitCluster()
    }else{
      testResult<-NULL
      for(testIter in 1:length(testList)){
        testResult<-rbind(testResult,paraCAT(testIter,testList,otutable,taxonomy,metric,metaData,outcomeVar,numBS,origR2,origBS,tree))
      }
    }
  }else{
    Ks<-compDistList(otutable,metric,tree)
    origR2<-MiRKATR2(metaData,outcomeVar,Ks)
    origBS<-rep(NA,numBS)
    for(BSiter in 1:numBS){
      set.seed(BSiter)
      indi<-sample(1:nrow(metaData),nrow(metaData),replace = TRUE)
      outcomeBS<-metaData[indi,outcomeVar]
      KsBS<-Ks
      temp<-NA
      for(k in 1:length(Ks)){
        for(i in 1:length(outcomeBS)){
          for(j in 1:length(outcomeBS)){
            KsBS[[k]][i,j]<-Ks[[k]][indi[i],indi[j]]
          }
        }
      }
      origBS[BSiter]<-MiRKATR2(outcomeBS,outcomeVar,KsBS)
    }
    if(parallel){
      registerDoParallel(nCore)
      testResult<-foreach(seedNum=1:length(testList),.combine='rbind',
                          .packages=c("ape","vegan","GUniFrac","MiRKAT")
      ) %dopar% paraCATMiRKAT(seedNum,testList,otutable,taxonomy,metric,metaData,outcomeVar,numBS,origR2,origBS,tree)
      stopImplicitCluster()
    }else{
      testResult<-NULL
      for(testIter in 1:length(testList)){
        testResult<-rbind(testResult,paraCATMiRKAT(testIter,testList,otutable,taxonomy,metric,metaData,outcomeVar,numBS,origR2,origBS,tree))
      }
    }
  }
  rownames(testResult)<-testList
  colnames(testResult)<-c("original R2","After remove taxon R2","p-value")
  testResult
}


