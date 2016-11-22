mid.tree=function(n,scale,qtile,Mean=F,reps=100){
  tb=numeric(n)
  for(i in 1:reps){
    tree.temp=pbtree(n=n,scale=scale)
    if(Mean==T){
      tb[i]=mean(tree.temp$edge.length)
    }
    else{
      tb[i]=quantile(tree.temp$edge.length,qtile)
    }
  }
  return(tb)
}



bayesf=function(a)
{
  return(a[,2]*(1-a[,1])/a[,1]/(1-a[,2]))
}


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
    res2 = c(which.max((x$prob2[,2]>=a)) -1,which.min((x$prob2[,2]<=b))-1)
    return(list(prior=res1,posterior=res2))
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

##################################







getcrossingnum=function(test,n1,n2){
  return(apply(test,2,function(x) return(mean(x<=n2 & x>=n1))))
}



#####################################
# here, upcrossing mean from 1 to 2

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


################
Browntip_simu=function(phy,epsilon,root.prior){
  n=phy$Nnode+1
  data=brown_tree_prior_node(phy,root.prior)
  tip=ifelse(data[1:n]>epsilon,"B",ifelse(data[1:n]< -epsilon,"A","C"))
  tip01=ifelse(tip=="B",1,ifelse(tip == "A",0,0.5))
  od=order(as.numeric(gsub("t", "", phy$tip.label)))
  names(tip01) = phy$tip.label[od]
  names(tip) = phy$tip.label[od]
  res=fills_in(phy,data,epsilon,0,0,F)
  return(list(phy=phy,tup=res,tip=tip,tip01=tip01))
}



brpr=function(phy150,eps,x0s,rep=2000){
  k=1
  print(k)
  result=matrix(ncol=(rep+2),nrow=(length(eps)*length(x0s)))
  x0s=sort(x0s,decreasing = T)
  for(i in 1:length(eps))
  {
    for(j in 1:length(x0s)){
      x0=x0s[j]*eps[i]
      res=brown_tree_prior2(phy150,0,eps[i],c(x0,x0/10),nsize = rep)
      if(mean(res<=10)<0.75) break
      else{
        result[k,]=c(x0s[j],eps[i],res)
        k=k+1
        print(k)
        cat("i=",i,"j=",j)
      }
    }
  }
  return(result)
}




fitdisctmc=function(phy150,d){
  CTMC_test_br1 = testIndOrigin(inputTrees=c(phy150), traitData=d,
                                initLambda01=.01, initLambda10=.01, priorAlpha01=1, priorBeta01=10,
                                priorAlpha10=1, priorBeta10=10, mcmcSize=21000, mcmcBurnin=10000,
                                mcmcSubsample=20, mcSize=10000)
  colnames(CTMC_test_br1$mcmcOutput)
  mxid=which.max(CTMC_test_br1$mcmcOutput[,3])
  rts=CTMC_test_br1$mcmcOutput[mxid,c(4,5)]

  DS_res=MLphy2(phy150,dbrsim1$tip01,rts[1], rts[2],max_num_jumps = 20)
  plot(c(0:20),DS_res$prob1[,1],type="h")
  return(qtile(DS_res,a=0,b=0.9))
}

####
bfk=function(pr,post,k){
  p1=mean(pr<=k)
  p2=mean(post<=k)
  return(p2*(1-p1)/p1/(1-p2))
}

