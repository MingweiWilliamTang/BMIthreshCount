#' Posterior probability for the number of crossing given tip binary observation
#'
#' @param phy A phylo type object.
#' @param d the discrete observation on the tip of the tree.
#' @param epsilon the radius of the interval for crossing
#' @param root.prior A numeric vector with two values, which is the mean and variance for the prior distribution (normal) for the root liability.
#' @return A vector associate with the posterior probabililty for the number of upcrossing

brownian_crossing_post2=function(phy,d,thed,epsilon,root.prior,ngen=20000,burnin=10000,thin=5,k=0,less=F){
  data=updatenode(phy,d,thed,root.prior,ngen,burnin,thin)
  cros=fills_in_mat(phy$edge,phy$edge.length,phy$Nnode,data,epsilon,0,0,F)
  return(cros)
}

#' Posterior probability for the number of crossing given tip binary observation
#'
#' @param phy A phylo type object.
#' @param d the discrete observation on the tip of the tree.
#' @param epsilon the radius of the interval for crossing
#' @param root.prior A numeric vector with two values, which is the mean and variance for the prior distribution (normal) for the root liability.
#' @return A vector associate with the posterior probabililty for the number of upcrossing

##################################
brownian_crossing_post3=function(phy,d,thed,epsilon,root.prior,ngen=30000,
                                 burnin=10000,thin=5,k=0,less=F){
  data=posterior_update(phy$edge,phy$edge.length,phy$Nnode,d,thed,epsilon,root.prior,ngen,burnin,thin)
  cros=fills_in_mat(phy$edge,phy$edge.length,phy$Nnode,data,epsilon,0,0,F)
  return(cros)
}


###############################
#' Posterior probability for the number of crossing given the laten value on each node
#'
#' @param d the discrete observation on the tip of the tree.
#' @param phy A phylo type object.
#' @param epsilon the radius of the interval for crossing
#' @param root.prior A numeric vector with two values, which is the mean and variance for the prior distribution (normal) for the root liability.
#' @return A vector associate with the posterior probabililty for the number of upcrossing

fastpost2=function(data,phy,epsilon,k=0){
  cros=fills_in_mat(phy$edge,phy$edge.length,phy$Nnode,data,epsilon,0,k,F)
  return(cros)
}
