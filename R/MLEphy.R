#' Calculate the two-state CTMC likelihood
#'
#' @param par Positive numeric vector for the transition rates
#' @param data List phy for phylogenetic tree, tip01 by binary states encoded in 0,1
#' @return CTMC likelihood
#' @export
#'
MLEphy=function(par,data){
  data$phy=reorder(data$phy,order="pr")
  a=TwoStatePhyloLikelihood1(data$phy$edge,data$tip01,
                             data$phy$edge.length,par[1],par[2],c(1,0))
  return(-a)
}

#' Calculate the two-state CTMC likelihood with rate in log parametrization
#'
#' @param par Numeric vector for the log transition rates
#' @param data List phy for phylogenetic tree, tip01 by binary states encoded in 0,1
#' @return CTMC likelihood
#' @export

MLEphy.log=function(par,data){
  data$phy=reorder(data$phy,order="pr")
  a=TwoStatePhyloLikelihood1(data$phy$edge,data$tip01,
                             data$phy$edge.length,exp(par[1]),exp(par[2]),c(1,0))
  return(-a)
}

#'
#'
#'

getprs=function(phy,epsilon,root.prior,k=10){
  bb=brown_tree_prior2(phy,0,epsilon,root.prior,T,nsize=2000)
  cumpr=sapply(c(0:k),function(x) return(mean(bb<=x)) )
  pr=sapply(c(0:k),function(x) return(mean(bb==x)) )
  par(mfrow=c(2,1))
  plot(c(0:k),pr,main="probability mass plot",xlab="upcross",ylab="prob",type="h")
  plot(c(0:k),cumpr,main="cummulative probability plot",xlab="upcross",ylab="prob",type="h")
}



CTMCprob1=function(data,maxnum=20){
  p=optimx(c(0.1,0.1),MLEphy,lower=c(0,0),data=data)
  print(c(p$p1,p$p2))
  prior1=MLphy2(data,p$p1,p$p2,maxnum,md="ARD")
  return(prior1$prob1[,1])
}

CTMCprob2=function(data,maxnum=10){
  p=optim(c(0.1,0.1),MLEphy.log,data=data)$par
  p = exp(p)
  print(p)
  prior1=MLphy2(data,p[1],p[2],maxnum,md="ARD")
  return(prior1$prob1[,1])
}



ckprpostbf=function(phy,d,epsilon,root.prior,nsize=5000,k)
{
  bb = brown_tree_prior2(phy,0,epsilon,root.prior,up=T,nsize)
  bbp = brownian_crossing_post3(phy,d,0,epsilon,root.prior,ngen=30000,
                                burnin=10000,thin=5,k=0,less=F)
  return(list(prior=mean(bb<=k),post=mean(bbp<=k),bf=bfk(bb,bbp,k)))
}


#################################
