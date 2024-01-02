compDistList<-function(otutable,metric,tree=NULL){
  Ks<-vector("list", length = length(metric))
  names(Ks)<-metric
  if(sum(metric %in% c("Weighted UniFrac","Unweighted UniFrac","robust"))>0){
    temp<-GUniFrac(t(otutable),tree)
    if(sum(metric%in% "Weighted UniFrac")>0){
      Ks[["Weighted UniFrac"]]<-D2K(temp$unifracs[,,"d_1"]) 
    }
    if(sum(metric%in% "Unweighted UniFrac")>0){
      Ks[["Unweighted UniFrac"]]<-D2K(temp$unifracs[,,"d_UW"])   
    }
    if(sum(metric%in% "robust")>0){
      distMat1<-temp$unifracs[,,"d_UW"]
      distMat1<-distMat1/max(distMat1)
      distMat2<-as.matrix(vegdist(t(otutable),"bray"))
      distMat2<-distMat2/max(distMat2) 
      Ks[["robust"]]<-D2K(0.5*distMat1+0.5*distMat2)
    } 
  }
  
  for(i in 1:length(metric)){
    if(!metric[i]%in%c("Weighted UniFrac","Unweighted UniFrac","robust")){
      Ks[[metric[i]]]<-D2K(as.matrix(vegdist(t(otutable),metric[i])))
    }
  }
  Ks
}