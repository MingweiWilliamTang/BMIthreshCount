#' Prior probability for the number of crossing on phylogenetic tree
#'
#' @param phy A phylo type object.
#' @param thed the threshold value for Brownian motion threshold model.
#' @param epsilon the radius of the interval for crossing.
#' @param root.prior A numeric vector with two values, which is the mean and variance for the prior distribution (normal) for the root liability.
#' @param up Previous crossing direction. Set as TRUE for defaulted value.
#' @param nsize number of simulations for Monte-Carlo integral
#' @return A vector associate with the prior probabililty for the number of upcrossing

brown_tree_prior2=function(phy,thed=0,epsilon,root.prior,up=T,nsize=1000,k=0){
  data=brown_tree_prior_node_cpp(phy$edge,phy$edge.length,phy$Nnode,root.prior,nsize)
  res=fills_in_mat(phy$edge,phy$edge.length,phy$Nnode,data,epsilon,0,k,F)
  return(res)
}



brown_tree_prior_data=function(phy,up=T,nsize=2000,k=0){
  root.prior = c(0,0.000000001)
  data=brown_tree_prior_node_cpp(phy$edge,phy$edge.length,phy$Nnode,root.prior,nsize)
  return(data)
}

#' Generating libility value on each node of the phylogenetic tree
#'
#' @param phy A phylo type object.
#' @param up Previous crossing direction. Set as TRUE for defaulted value.
#' @param nsize number of simulations for Monte-Carlo integral
#' @return A vector associate with the prior probabililty for the number of upcrossing

brown_tree_prior_grid = function(phy,grid,Data,k=0){
  # grid is (epsilon,x0)
  #  print(Data)
  nrep = dim(Data)[1]
  p = dim(Data)[2]
  res_mat = matrix(nrow = length(grid[,1]),ncol = nrep)
  for(i in 1:length(grid[,1])){
    x0 = grid[i,2]
    epsilon = grid[i,1]
    root =rnorm(nrep,x0,x0/10)
    data2 = Data + matrix(rep(root,p),nrow=nrep,byrow = T)
    res=fills_in_mat(phy$edge,phy$edge.length,phy$Nnode,data2,epsilon,0,k,F)
    res_mat[i,] = res
  }
  return(res_mat)
}

#' Generating libility value on each node of the phylogenetic tree
#'
#' @param phy A phylo type object.
#' @param up Previous crossing direction. Set as TRUE for defaulted value.
#' @param nsize number of simulations for Monte-Carlo integral
#' @return A vector associate with the prior probabililty for the number of upcrossing

brown_tree_prior_node = function(phy,root.prior){
  p=length(phy$tip.label)
  N=p+phy$Nnode
  x0=rnorm(1,as.numeric(root.prior[1]),as.numeric(root.prior[2]))
  data=numeric(N)
  data[p+1]=x0
  start=phy$edge[,1]
  ts=phy$edge.length
  #  upd=logical(N)
  #  upd[p+1]=ifelse(x0>epsilon,1,0)
  #  cros=NULL
  for(i in 1:dim(phy$edge)[1]){
    j=phy$edge[i,2]
    data[j]=rnorm(1,data[start[i]],sqrt(ts[i]))
    #res=brown_crossing_sample(data[start[i]],data[j],ts[i],epsilon,up=upd[start[i]],nsize=1)
    #cros=c(cros,res$samples)
    #upd[j]=res$upnext
  }
  return(data)
}
#############

#' Generating the number of crossings given the nodes liabilties of the phylogenetic tree
#'
#' @param phy A phylo type object.
#' @param data The liabilities of the phylogenetic tree nodes.
#' @param epsilon the radius of the interval for crossing.
#' @param thed the threshold value for Brownian motion threshold model.
#' @param cut cutoff value: During the accumulating process, if the number of crossing is greater than the cutoff value, it will stop this iteration and yeild greater than cutoff.
#' @return A vector associate with the prior probabililty for the number of upcrossing

fills_in=function(phy,data,epsilon,thed=0,cut=0,less=F){
  N=2*phy$Nnode+1
  p=phy$Nnode+1
  x0=data[p+1]
  start=phy$edge[,1]
  ts=phy$edge.length
  upd=logical(N)
  upd[p+1]=ifelse(x0>-epsilon,1,0)
  cros=0
  rt=0
  data=as.numeric(data)
  #### fill in all the nodes
  if(less==T){
    for(i in 1:dim(phy$edge)[1]){
      j=phy$edge[i,2]
      res=brown_crossing_sample(data[start[i]],data[j],ts[i],epsilon,up=upd[start[i]],nsize=1)
      cros=cros+res$samples
      if(cros>cut){
        rt=1
        break
      }
      upd[j]=res$upnext
    }
    return(rt)
  }
  else{
    for(i in 1:dim(phy$edge)[1]){
      j=phy$edge[i,2]
      res=brown_crossing_sample(data[start[i]],data[j],ts[i],epsilon,up=upd[start[i]],nsize=1)
      cros=cros+res$samples
      upd[j]=res$upnext
    }
    return(cros)
  }
}
