# roc functions
plot.roc=function(s2,ntrial,N_study=NULL,cex=1,pch=3){
  # plot roc curve of the experiment
  # and return the auc value
  par(mar=c(2,2,2,0)+1)
  T_total = s2$Total[1]
  F_total = ntrial - T_total
  plot(c(1,s2$count[,2]/F_total),c(1,s2$count[,1]/T_total),type="l",col="black",lwd = 2,main=paste('N_study=',as.character(N_study)),
       xlab='False positive rate',ylab='True positive rate',xlim=c(0,1),ylim=c(0,1),mgp=c(1.2,0.3,0),cex.lab=1.2, cex.axis=1.2)
  points(s2$count[,2]/F_total,s2$count[,1]/T_total,pch = pch,cex = cex)
  lines(c(0,1),c(0,1),lty=2)
  # return auc value
  l = length(s2$count[,1])+2
  sc1=c(T_total,s2$count[,1],0)/T_total
  height = (sc1[-1]+sc1[-l])/2
  width = -diff(c(F_total,s2$count[,2],0))/F_total
  auc = sum(height*width)
  return(auc)
}

lines.roc=function(s2,ntrial,color="black",lty=2,pch=2,cex=1){
  # plot roc curve of the experiment
  # and return the auc value
  T_total = s2$Total[1]
  F_total = ntrial - T_total
  lines(c(1,s2$count[,2]/F_total),c(1,s2$count[,1]/T_total),type="l",col=color,lwd=2,lty=lty)
  points(s2$count[,2]/F_total,s2$count[,1]/T_total,pch = pch,cex=cex)
  lines(c(0,1),c(0,1),lty=2)
  # return auc value
  l = length(s2$count[,1])+2
  sc1=c(T_total,s2$count[,1],0)/T_total
  height = (sc1[-1]+sc1[-l])/2
  width = -diff(c(F_total,s2$count[,2],0))/F_total
  auc = sum(height*width)
  return(auc)
}

roc.point.specify=function(s2,rfp0,rtp0=NULL,ntrial,levels){
  T_total = s2$Total[1]
  F_total = ntrial - T_total
  fp = s2$count[,2]/F_total
  tp = s2$count[,1]/T_total

  if(is.null(rtp0)){
    rseq = fp
    test = rfp0
  }else{
    rseq = tp
    test = rtp0
  }
  for(i in length(levels):1){
    if(rseq[i]>=test) break
  }
  if(rseq[i]<test){
    return("need more levels")
  }else{
    return(ifelse(is.null(rtp0),list(level=levels[i],tp=tp[i]),
                  list(level = levels[i],fp=fp[i])))
  }
}
