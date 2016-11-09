#' Counting the number of Crossing of a process
#'
#'
#' @author Mingwei Tang
#' @param path  a vector of a stochastic process over a discrete grid
#' @param thed  the threshold value
#' @param epsilon  half length of the interval of the threshold, that is (thed-epsilon,thed+epsilon)
#' @return crossing: the number of crossing for a interval of the threshold


crossing=function(pth,thed=0,epsilon=0,prev=0,t=1){
  if(epsilon<0) stop('epsilon<0')
  n=length(pth$y)
  up=thed+epsilon
  down=thed-epsilon
  path = if(prev==1){
    path = c(down-0.5,pth$y)
  }
  else{
    path = c(up + 0.5,pth$y)
  }
  id=(path>up)-(path<down)
  id1=id>0
  id3=id<0
  # upc=sum(((id[2:n]-id[1:(n-1)])>0)&id1[2:n]>0)
  # downc=sum(((id[2:n]-id[1:(n-1)])<0)&id3[2:n]>0)
  nz_id = which(id!=0)
  id_non_zero=id[nz_id]
  l=length(id_non_zero)
  if(l<2){
    crossing=0
    upc=0
    downc=0
    px=NA
    py=NA
  }
  else{
    # print(l)
    crossing=sum(id_non_zero[2:l]*id_non_zero[1:(l-1)]<0)
    cr_id = which(id_non_zero[2:l]*id_non_zero[1:(l-1)]<0)
    if(length(cr_id)==0){
      px=NA
      py=NA
    }
    else{
      px=t*(nz_id[cr_id+1]-1)/(n-1)
      py=path[nz_id[cr_id+1]]
    }
    upc=sum((id_non_zero[2:l]*id_non_zero[1:(l-1)]<0)*(id_non_zero[1:l-1]<0))
    downc=sum((id_non_zero[2:l]*id_non_zero[1:(l-1)]<0)*(id_non_zero[1:l-1]>0))
  }

  return(list(crossing=crossing,upcrossing=upc,downcrossing=downc,xpoints=px,ypoints=py))
}


#' simulate the sample path of OU-process or Brownian motion
#'
#' @param x0 the starting point at t=0
#' @param TT the total time
#' @param n the number of the grid
#' @param theta Default as 0, which means brownian motion
#' @param sigma The diffusion parameter of OU process
#' @param pl Whether plot the path, default as false
#' @return A vector of OU-process path or Brownian motion path
OUpath=function(x0,TT,n,theta=0,mu,sigma=1,pl=F){
  Xt=rep(0,n)
  dt=TT/n
  grid=dt*c(1:n)
  if(theta!=0){
    mu1=(1-exps(theta,dt))*mu
    sigma_temp=sqrt(sigma^2/(2*theta)*(1-exps(2*theta,dt)))
    for(i in 1:n){
      if(i==1) mu_temp=mu1+x0*exps(theta,dt)
      else mu_temp=mu1+Xt[i-1]*exps(theta,dt)
      Xt[i]=mu_temp+sigma_temp*rnorm(1)
    }
  }
  else{
    sigma_temp=sigma*sqrt(dt)
  }
  for(i in 1:n){
    if(i==1) mu_temp=x0
    else mu_temp=Xt[i-1]
    Xt[i]=mu_temp+sigma_temp*rnorm(1)
  }
  if(pl==T){
    grid=c(0,grid)
    Xt=c(x0,Xt)
    plot(grid,Xt,type="l",xlab="time",ylab="liability value")
  }
  return(Xt)
}



exps=function(theta,x){
  return(exp(-theta*(x)))
}

#' simulating and counting the sample path of OU bridge or Brownian bridge
#'
#' @param x0 the starting point at t=0
#' @param XT the ending point at t=TT
#' @param TT the total time
#' @param n the number of the grid
#' @param theta Default as 0, which means brownian motion
#' @param sigma The diffusion parameter of OU process
#' @param pl Whether plot the path, default as false
#' @return a vector of OU-bridge
OUbridge=function(x0,xT,TT,n,theta,mu,sigma,pl=F){
  Xt=rep(0,n)
  dt=TT/n
  grid=dt*c(1:n)
  if(theta==0){
    sigma_temp=sigma*sqrt(dt)
    for(i in 1:n){
      if(i==1) mu_temp=x0-mu
      else mu_temp=Xt[i-1]
      Xt[i]=mu_temp+sigma_temp*rnorm(1)
    }
    Zt=Xt+(xT-mu-Xt[n])*grid/TT
  }
  else{
    sigma_temp=sqrt(sigma^2/(2*theta)*(1-exps(2*theta,dt)))
    for(i in 1:n){
      if(i==1) mu_temp=exps(theta,dt)*(x0-mu)
      else mu_temp=exps(theta,dt)*Xt[i-1]
      Xt[i]=mu_temp+sigma_temp*rnorm(1)
    }
    Zt=Xt+(xT-mu-Xt[n])*(exps(-theta,grid)-exps(theta,grid))/(exps(-theta,TT)-exps(theta,TT))
  }
  Zt=Zt+mu
  Zt=c(x0,Zt)
  if(pl==T){
    grid=c(0,grid)
    plot(grid,Zt,type="l",ylab="Liability Xt",xlab="time",ylim=c(min(Zt)-0.3,max(Zt)+0.1),mgp=c(2.5,1,0.5),cex.lab=2,las=-0.5,cex.axis=1.5)
    #abline(v=epsilon,col="red",lty=2)
    #abline(v=-epsilon,col="red",lty=2)
  }
  return(list(x=grid,y=Zt))
}
