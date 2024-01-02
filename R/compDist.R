compDist<-function(otutable,metric,tree=NULL){
  if(metric %in% c("Weighted UniFrac","Unweighted UniFrac","robust")){
    temp<-GUniFrac(t(otutable),tree)
    if(metric=="Weighted UniFrac"){
      distMat<-as.dist(temp$unifracs[,,"d_1"])         
    }else if(metric=="Unweighted UniFrac"){
      distMat<-as.dist(temp$unifracs[,,"d_UW"])   
    }else{
      distMat1<-as.matrix(as.dist(temp$unifracs[,,"d_UW"]))
      distMat1<-distMat1/max(distMat1)
      distMat2<-as.matrix(vegdist(t(otutable),"bray"))
      distMat2<-distMat2/max(distMat2) 
      distMat<-0.5*distMat1+0.5*distMat2
    }
  }else{
    distMat<-vegdist(t(otutable),metric)    
  }
}