#' Sample the number of gaining traits on a line based on CTMC
#'
#' @param x0 Start point from CTMC
#' @param t Total length of time
#' @param k12 transition rate from state 0 to state 1
#' @param k21 transition rate from state 1 to state 0
#' @return number of gain traits and the ending state

CTMC.line.sample=function(x0,t,k12,k21){
  # t total length of the time
  # return of dataframe has three colomns:
  # the 1st col is the jumping time, the 2nd col is the states

  # construct the transition rate
  ks = c(k12,k21)

  # generate the initial states
  state = numeric(1)
  state[1] = x0
  ts=numeric(1)
  i=1

  while( ts[i]<t ){
    k=ks[state[i]]
    t.temp=rexp(n = 1, rate = k)
    ts=c(ts,ts[i]+t.temp)
    state=c(state,ifelse(state[i]==1,2,1))
    i=i+1
  }
  n = length(state) -1
  nup=sum(diff(state[1:n])==1)
  return(list(nup=nup,end=state[n]))
}


CTMC.tree=function(phy,k12,k21){
  # on a given tree, simulate the number of upcrossings
  # returns the tip node and number of upcrossing
  p = phy$Nnode + 1
  n = 2*p - 1
  edge = phy$edge
  d1 = dim(edge)[1]
  edge.length = phy$edge.length
  # root node be 2
  states=numeric(n)
  states[p+1] = 2
  tup = 0
  for(i in 1:d1){
    temp = CTMC.line.sample(states[edge[i,1]],edge.length[i],k12,k21)
    states[edge[i,2]] = temp$end
    tup = tup + temp$nup
  }
  tip = states[1:p]-1
  tip01=tip
  tip = ifelse(tip ==1,"B","A")
  od=order(as.numeric(gsub("t", "", phy$tip.label)))
  names(tip01) = phy$tip.label[od]
  names(tip) = phy$tip.label[od]
  return(list(phy=phy,tup=tup,tip=tip,tip01=tip01))
}



tree.plot=function(phy.data){
  phy = phy.data$phy
  d = phy.data$tip
  n.simu = phy.data$tup
  par(mar=c(2,2,2,2))
  plot(phy,show.tip.label=FALSE,main=paste(n.simu," 0-1 crossings"))
  col <- c("#004165", "#eaab00")
  d=ifelse(d=="A",1,0)
  d=d[order(as.numeric(gsub("t", "", phy$tip.label)))]
  tiplabels(col=col[d+1], pch=20, adj=0.505)
  legend("topleft",col=c("#004165", "#eaab00"),pch=20,legend = c("1","0"))
}
#########################################



#brownian_crossing_post3=function(phy,d,thed,epsilon,root.prior,ngen=30000,
#                                 burnin=10000,thin=5,k=0,less=F){
#  data=posterior_update(phy$edge,phy$edge.length,phy$Nnode,d,thed,root.prior,ngen,burnin,thin)
#  cros=fills_in_mat(phy$edge,phy$edge.length,phy$Nnode,data,epsilon,0,0,F)
#  return(cros)
#}


MLphy=function(phy150.data1,max_num_jumps=10,md="ARD"){
  phy150=phy150.data1$phy
  d=phy150.data1$tip
  ER=fitDiscrete(phy150,d,model=md)

  MQ=matrix(c(ER$opt$q12*c(-1,1),ER$opt$q21*c(1,-1)),byrow = T,ncol=2)
  print(MQ)

  log_probabilities = treeConvolveTest(phy150, phy150.data1$tip01,
                                       ER$opt$q12, ER$opt$q21, root.dist.prob.0 = 0.01, max_num_jumps)

  # log_probabilities = treeConvolveTest(phy150, ifelse(d=="A",0,1), 0.01,10,0.1,max_num_jumps)

  rownames(log_probabilities) = c(paste(c(0:(max_num_jumps)),"0->1 jumps"))

  prob1=exp(log_probabilities)
  prob2=cbind(cumsum(exp(log_probabilities)[,1]),cumsum(exp(log_probabilities)[,2]))
  bf=bayesf(prob2)
  return(list(prob1=prob1,prob2=prob2,bayesfactor=bf))
}



MLphy2=function(phy150.data1,r1,r2,max_num_jumps=10,md="ARD"){
  #  ER=fitDiscrete(phy150,d,model=md)

  #  MQ=matrix(c(ER$opt$q12*c(-1,1),ER$opt$q21*c(1,-1)),byrow = T,ncol=2)
  # print(MQ)
  phy150=phy150.data1$phy
  d=phy150.data1$tip01
  log_probabilities = treeConvolveTest(phy150, d,
                                       r1, r2, root.dist.prob.0 = 0, max_num_jumps)

  # log_probabilities = treeConvolveTest(phy150, ifelse(d=="A",0,1), 0.01,10,0.1,max_num_jumps)

  rownames(log_probabilities) = c(paste(c(0:(max_num_jumps)),"0->1 jumps"))

  prob1=exp(log_probabilities)
  prob2=cbind(cumsum(exp(log_probabilities)[,1]),cumsum(exp(log_probabilities)[,2]))
  bf=bayesf(prob2)
  return(list(prob1=prob1,prob2=prob2,bayesfactor=bf))
}



qtile=function(x,type="quantile",a=0.25,b=0.75){
  if(type=="mode"){
    n1=which.max(x$prob1[,1]) - 1
    n2=which.max(xprob1[,2]) - 1
    res=c(n1,n2)
    names(res)=c("prior","posterior")
    return(res)
  }
  else
  {
    #   n1= (x$prob2[,1]>=a) & (x$prob2[,1]<=b)
    #   n2= (x$prob2[,2]>=a) & (x$prob2[,2]<=b)
    res1 = c(which.max((x$prob2[,1]>=a)) -1,which.min((x$prob2[,1]<=b))-1)
    #    res2 = c(which.max((x$prob2[,2]>=a)) -1,which.min((x$prob2[,2]<=b))-1)
    return(prior=res1)
  }
}

pr_test_epsilonx0=function(phy,qts,gds,reps=2000){
  epsilons=sqrt(quantile(phy$edge.length,qts))*0.125
  x0s=as.vector(sapply(epsilons,function(x) return(x*gds)))
  epsilons=rep(epsilons,each=length(gds))
  test.mtr.1=data.frame(x0s,epsilons)
  print(test.mtr.1)
  simu.2.res.1=apply(test.mtr.1,1,function(x) return(brown_tree_prior1(phy,0,x[2],c(x[1],x[1]/10),T,reps)))
  return(list(design=test.mtr.1,count.mx=simu.2.res.1))
}

pr_test_epsilonx0_2=function(phy,test.mtr.1,reps=2000){
  print(test.mtr.1)
  simu.2.res.1=apply(test.mtr.1,1,function(x) return(brown_tree_prior1(phy,0,x[2],c(x[1],x[1]/10),T,reps)))
  return(list(design=test.mtr.1,count.mx=simu.2.res.1))
}

post_test_epsilonx0=function(phy,d,test.mx){
  res=apply(test.mx,1,function(x){
    data=updatenode(phy,d,0,c(x[1],x[1]/10),20000,10000,5)
    return(fastpost(data,phy,x[2],0))
  }
  )
  return(res)
}


MLEphy=function(par,data){
  data$phy=reorder(data$phy,order="pr")
  a=TwoStatePhyloLikelihood1(data$phy$edge,data$tip01,
                             data$phy$edge.length,par[1],par[2],c(1,0))
  return(-a)
}



getprs=function(phy,epsilon,root.prior,k=10){
  bb=brown_tree_prior2(phy,0,epsilon,root.prior,T,nsize=2000)
  cumpr=sapply(c(0:k),function(x) return(mean(bb<=x)) )
  pr=sapply(c(0:k),function(x) return(mean(bb==x)) )
  par(mfrow=c(2,1))
  par(mar=c(2.3,2.5,2,0))
  par(oma=c(0,0,0,0))
  barplot(pr,width = rep(0.5,k),space=1.2,names.arg=as.character(c(0:k)),main="probability mass plot",xlab="upcross",
          ylab="prob",mgp=c(1.2,0.3,0),ylim=c(0,1),cex.axis=1.5,cex.names = 1.5,cex.lab=1.5)
  par(mar=c(2.3,2.5,1,0))
  barplot(c(cumpr),width = rep(0.5,k),space=1.2,names.arg=as.character(c(0:k)),main="cummulative probability plot",
          xlab="upcross",ylab="prob",mgp=c(1.2,0.3,0),ylim=c(0,1),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5)
}


CTMCprob1=function(data,maxnum=10){
  p=optimx(c(0.1,0.1),MLEphy,lower=c(0,0),data=data)
  prior1=MLphy2(data,p$p1,p$p2,maxnum,md="ARD")
  return(prior1$prob1[,1])
}

CTMCprob2=function(data,maxnum=10){
  p=optimx(c(0.1,0.1),MLEphy,lower = c(0,0),data=data)$par
  prior1=MLphy2(data$phy,data$tip01,p[1],p[2],maxnum,md="ARD")
  return(prior1$prob2[,1])
}


###
qtile2 = function(cum,p){
  for(i in (1:length(cum))){
    if(cum[i]>=p) break
  }
  if(i>1&&cum[i-1]>=0.95*p){
    return(i-2)
  }else{
    return(i-1)}
}

####

CTMC_BM_Connect = function(phy.data,pr_count_mx,grid,alpha,c1=0.75,c2=0.7,c3=1.2){
  # get information from CTMC model
  m1 = 40
  prior_CTMC=CTMCprob1(phy.data,m1)

  if(is.na(max(prior_CTMC))){
    N_test = 0
  }else{
    cum_CTMC = cumsum(prior_CTMC)
    #    while(max(cum_CTMC)<alpha && m1<=100){
    #      m1 = m1+10
    #     cum_CTMC=CTMCprob1(phy.data,m1)
    #    }
    N_test = qtile2(cum_CTMC,alpha)

    if(N_test == 39)
      alpha=max(cum_CTMC)
  }

  # link that with BM model
  lable = apply(pr_count_mx,1,function(x) {
    if(N_test==0){
      return(mean(x<=N_test)>=c1)
    }
    else{
      return(mean(x<=N_test)>alpha*c2 &&mean(x<=N_test) < alpha*c3)
    }
  })

  if(sum(lable)==0){return(0)
  }else{
    post_grid = grid[lable,]
    print(post_grid)
    BMpost_cross=apply(post_grid,1,function(x) return(brownian_crossing_post3(phy.data$phy,phy.data$tip01,
                                                                              0,x[1],c(x[2],x[2]/10),ngen=30000,burnin=10000,thin=5,less=F)))


    return(list(prmx=pr_count_mx[lable,],posmx=t(as.matrix(BMpost_cross))))
    # return of list of two matrices
    # prmx is the count related to prior, each row is 2000 countings at certain (epsilon,x0)
    # posmx is the count related to posterior, each row is 3000 countings at certain (epsilon,x0)
  }

}

##############33

#' link CTMC with brownian motion
#'
#'

mx_counts = function(phy.data,latent_data,alpha,grid)
{
  # phy.data = both tree and tip data
  # alpha = the probability to determine N that links CTMC with BM
  # epsilonx0.grid = a grid of (epsilon,x0) pair that use to analyze
  # the output is a grid of (epsilon,x0)

  m1 = 40
  prior_CTMC=CTMCprob1(phy.data,m1)

  if(is.na(max(prior_CTMC))){
    N_test = 0
  }else{
    cum_CTMC = cumsum(prior_CTMC)
    #    while(max(cum_CTMC)<alpha && m1<=100){
    #      m1 = m1+10
    #     cum_CTMC=CTMCprob1(phy.data,m1)
    #    }
    N_test = qtile2(cum_CTMC,alpha)

    if(N_test == 39)
      alpha=max(cum_CTMC)
  }
  #  brown_tree_prior_data(phy.data$phy,up=T,nsize = 2000,k=0)

  pr_count_mx=brown_tree_prior_grid(phy.data$phy,grid,latent_data,k=0)
  lable = apply(pr_count_mx,1,function(x) {
    if(N_test==0){
      return(mean(x<=N_test)>=0.75)
    }
    else{
      return(mean(x<=N_test)>alpha*0.7 &&mean(x<=N_test) < alpha*1.2)
    }
  })

  if(sum(lable)==0){return(0)
  }else{
    post_grid = grid[lable,]

    BMpost_cross=apply(post_grid,1,function(x) return(brownian_crossing_post3(phy.data$phy,phy.data$tip01,
                                                                              0,x[1],c(x[2],x[2]/10),ngen=30000,burnin=10000,thin=5,less=F)))


    return(list(prmx=pr_count_mx[lable,],posmx=t(as.matrix(BMpost_cross))))
    # return of list of two matrices
    # prmx is the count related to prior, each row is 2000 countings at certain (epsilon,x0)
    # posmx is the count related to posterior, each row is 3000 countings at certain (epsilon,x0)
  }
}


#function(pr,post,k){
#  p1=mean(pr<=k)
#  p2=mean(post<=k)
#  if(p2==0) {
#    return(0)
#  }else if(p2==1){
#    return(p2*(1-p1)/p1/(1-p2))
#  }
#}





########
bf_v=function(listmx,N_study){
  prmx = listmx$prmx
  n1 = dim(prmx)[2]
  posmx = listmx$posmx
  bf_vec = apply(cbind(prmx,posmx),1,function(x) return(bfk(x[1:n1],x[-c(1:n1)],N_study)))
  bf_vec=na.omit(bf_vec)
  return(bf_vec)
}


positive_negative=function(bf_vec,level,N_real,N_study){
  TP = (mean(bf_vec>=level)>0.5) && (N_study>=N_real)
  FP = (mean(bf_vec>=level)>0.5) && (N_study<N_real)
  return(c(TP,FP))
}


######

positive_negative_rate_parallel=function(phy,ntrial,x.root,
                                         eps,grid,alpha,level,
                                         pr_count_mx,
                                         c1 = 0.75,c2 = 0.7,c3 = 1.2){

  count = matrix(rep(0,length(level)*8),ncol=8) # a vector that counts # of TP and # of FP
  T_total = c(0,0,0,0)

  for(i in 1:ntrial){
    #   print(paste('finish ',i))
    # simulate ntrial sets of tip values
    # in each trail we simulate the binary values on the tip of the tree
    # The rate is calculated by averaging the results over trails
    listmx=0
    bf_vec = TRUE
    rrr = 0
    while(!(is.list(listmx)) || is.logical(bf_vec) ){
      phy.data=Browntip_simu(phy,eps,c(x.root,x.root/10))
      # listmx=mx_counts(phy.data,latent_data,alpha,grid)
      rrr = rrr + 1
      listmx = CTMC_BM_Connect(phy.data,pr_count_mx,grid,alpha,c1,c2,c3)
      if(!(is.list(listmx))){
        next
      }else{
        bf_vec = c(bf_v(listmx,0),bf_v(listmx,1),bf_v(listmx,2),bf_v(listmx,3))
      }
    }
    # print(bf_vec)
    N_real = phy.data$tup
    for(k in 1:4){
      T_total[k] = T_total[k] + (k-1>=N_real)
    }
    for(j in 1:length(level)){
      for(k in 1:4){
        count[j,c(2*k-1,2*k)] = count[j,c(2*k-1,2*k)] +
          positive_negative(bf_vec[k],level[j],N_real,k-1)
      }
    }
  }
  return(list(count=count,Total=T_total,reps=rrr))
}


get_prior_count_mx=function(phy,x.root,eps,grid,nsize=10000){
  # generate latent data
  latent_data = brown_tree_prior_data(phy,up=T,nsize,k=0)
  pr_count_mx=brown_tree_prior_grid(phy,grid,latent_data,k=0)
  return(pr_count_mx)
}


#####
positive_negative_rate=function(phy,ntrial,x.root,
                                eps,grid,alpha,level,
                                N_study,c1 = 0.75,c2 = 0.7,c3 = 1.2){

  #set grid
  # grid is a matrix with each row a vector of (epsilon,x0)

  # level is a vector of crossing values that creates different situations in ROC curve

  # argument x0 and epsilon are the original number that used to simulate the data on the tip
  # in BM model

  # generate latent data
  latent_data = brown_tree_prior_data(phy,up=T,nsize=10000,k=0)

  count = matrix(rep(0,length(level)*2),ncol=2) # a vector that counts # of TP and # of FP
  T_total = 0
  pr_count_mx=brown_tree_prior_grid(phy,grid,latent_data,k=0)

  for(i in 1:ntrial){
    #   print(paste('finish ',i))
    # simulate ntrial sets of tip values
    # in each trail we simulate the binary values on the tip of the tree
    # The rate is calculated by averaging the results over trails
    listmx=0
    bf_vec = TRUE
    while(!(is.list(listmx)) || is.logical(bf_vec) ){
      phy.data=Browntip_simu(phy,eps,c(x.root,x.root/10))
      # listmx=mx_counts(phy.data,latent_data,alpha,grid)
      listmx = CTMC_BM_Connect(phy.data,pr_count_mx,grid,alpha,c1,c2,c3)
      if(!(is.list(listmx))){
        next
      }else{
        bf_vec = bf_v(listmx,N_study)
      }
    }
    # print(bf_vec)
    N_real = phy.data$tup
    T_total = T_total + (N_study>=N_real)
    for(j in 1:length(level)){
      count[j,] = count[j,] + positive_negative(bf_vec,level[j],N_real,N_study)
    }
  }
  return(list(count=count,Total=T_total))
}

################
CTMC.positive.negative.rate = function(phy,ntrial,x.root,eps,level,N_study){

  T_total = 0
  count = matrix(rep(0,length(level)*2),ncol=2) # a vector that counts # of TP and # of FP

  # a forloop get the result in each trial, the averge of the result is the rate we want
  for(i in 1:ntrial){

    # Firstly simulate data from BM model with threshold
    phy.data=Browntip_simu(phy,eps,c(x.root,x.root/10))
    N_real=phy.data$tup
    # Secondly run indorigin to compute Bayes factors
    CTMC_test_1 = testIndOrigin(inputTrees=c(phy), traitData=phy.data$tip01,
                                initLambda01=.01, initLambda10=.01, priorAlpha01=1, priorBeta01=10,
                                priorAlpha10=1, priorBeta10=10, mcmcSize=20000, mcmcBurnin=10000,
                                mcmcSubsample=20, mcSize=10000,testThreshold = N_study)
    BF1 = as.numeric(getBF(CTMC_test_1)[1])
    #  print(paste('BF=',BF1))
    # compare N_study with N_real to True or False
    T_total = T_total + (N_real<=N_study)
    # another loop get positive or negative for each level
    for(j in 1:length(level)){

      # compare Bayes factor with levels to get positive or negetive
      TPFP=c((BF1 >= level[j])&&(N_study>=N_real),(BF1 >= level[j])&&(N_study<N_real))
      #print(TPFP)
      count[j,] = count[j,] + TPFP
      #    print(paste('finish',i,'th trial'))
    }
    #  if(i %% 10 == 0){ print(paste('finish ',i,' iterations'))}
    # output the count of TP and FP
  }
  return(list(count=count,Total=T_total))
}

###
# consistant model

BM_positive_negative_rate = function(phy,ntrial,x.root,eps,level,N_study){
  T_total = 0
  thed=0
  count = matrix(rep(0,length(level)*2),ncol=2) # a vector that counts # of TP and # of FP

  # first compute the prior probability of crossing
  pr_vec = brown_tree_prior2(phy,thed=0,eps,c(x.root,x.root/10),up=T,nsize=2000,k=0)
  pr1 = mean(pr_vec<=N_study)
  # a forloop get the result in each trial, the averge of the result is the rate we want
  for(i in 1:ntrial){

    # Firstly simulate data from BM model with threshold
    phy.data=Browntip_simu(phy,eps,c(x.root,x.root/10))
    postemp=brownian_crossing_post3(phy,phy.data$tip01,thed,eps,c(x.root,x.root/10),ngen=30000,
                                    burnin=10000,thin=5,k=0,less=F)
    postpr=mean(postemp<=N_study)
    bf=postpr*(1-pr1)/pr1/(1-postpr)
    N_real = phy.data$tup
    T_total = T_total + (N_study>=N_real)
    for(j in 1:length(level)){
      TPFP=c((bf >= level[j])&&(N_study>=N_real),(bf >= level[j])&&(N_study<N_real))
      count[j,] = count[j,] + TPFP
    }
  }
  return(list(count=count,Total=T_total))
}

#' Use calibrated Brownian motion method on testing the hypothesis of transition numbers
#'
#' @param phy a phylogenetic tree object
#' @param tipvalue the values of the tip node
#' @param grid the grid set for testing hypothesis
#' @param c1
#' @param c2
#' @param c3
#'
#' @return The bayes factor output from the grid

BM_calibrated_test = function(phy,tipvalue,grid,alpha,c1=0.75,c2=0.7,c3=1.2){
  #   print(paste('finish ',i))
  # simulate ntrial sets of tip values
  # in each trail we simulate the binary values on the tip of the tree
  # The rate is calculated by averaging the results over trails
  latent_data = brown_tree_prior_data(phy,up=T,nsize=10000,k=0)
  pr_count_mx=brown_tree_prior_grid(phy,grid,latent_data,k=0)

  listmx=0
  bf_vec = TRUE
  phy.data = list(phy=phy,tip01=tipvalue)
  while(!(is.list(listmx)) || is.logical(bf_vec) ){

    # listmx=mx_counts(phy.data,latent_data,alpha,grid)
    listmx = CTMC_BM_Connect(phy.data,pr_count_mx,grid,alpha,c1,c2,c3)
    if(!(is.list(listmx))){
      next
    }else{
      bf_vec = list(test0=bf_v(listmx,0),test1 = bf_v(listmx,1), test2 = bf_v(listmx,2), test3 = bf_v(listmx,3))
    }
  }
  # print(bf_vec)
  return(bf_vec)
}
print("source complete")

