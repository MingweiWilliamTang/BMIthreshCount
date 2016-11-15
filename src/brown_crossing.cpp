//@copyright Mingwei
/*
#include<RcppArmadilloExtensions/sample.h>
#include<RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include<Rcpp.h>
#include<R.h>
#include<math.h>
using namespace Rcpp;
using namespace arma;
 */

#include "twostates.h"
#define pi           3.14159265358979323846
//[[Rcpp::export()]]
double crossing_brownbridge(double x0,double xt, double t ,int ncross,
                            double epsilon ,bool up= true ,bool counditional = true)
{
  int k = ceil( ((double) ncross)/2 );
  double pro = R::dnorm(xt, x0, sqrt(t),0);
  double a1 = (4*k-3)*epsilon;
  double a2 = (4*k-1)*epsilon;
  if (x0 < -epsilon) {
    up = false;
  }
  if (up == true) {
    if (ncross == 0) {
      return((xt> -epsilon)/sqrt(2*pi*t) *(exp(-pow(x0-xt,2)/(2*t))-exp(-pow(2*epsilon+xt+x0,2)/(2*t)))/pro);
    }
    if (ncross % 2 == 1) {
      return((xt<epsilon)/sqrt(2*pi*t)*(exp(-pow(x0+a1+fabs(epsilon+xt),2)/(2*t))-exp(-pow(x0+a1+fabs(xt-3*epsilon),2)/(2*t)))/pro);
    }
    if ( ncross % 2 == 0 && ncross > 0 ) {
      return((xt > -epsilon)/sqrt(2*pi*t)*(exp(-pow(x0+a2+fabs(epsilon-xt),2)/(2*t))- exp(-pow(x0+a2+fabs(xt+3*epsilon),2)/(2*t)))/pro);
    }
  }
  else{
    if(ncross==0){
      return((xt<epsilon)/sqrt(2*pi*t)*(exp(-pow(xt-x0,2)/(2*t))-exp(-pow(x0-2*epsilon+xt,2)/(2*t)))/pro);
    }
    if (ncross % 2 == 1) {
      return((xt > -epsilon)/sqrt(2*pi*t)*(exp(-pow(-x0+a1+fabs(epsilon-xt),2)/(2*t))-exp(-pow(-x0+a1+fabs(xt+3*epsilon),2)/(2*t)))/pro) ;
    }
    if(ncross % 2==0 && ncross>0){
      return((xt<epsilon)/sqrt(2*pi*t)*(exp(-pow(-x0+a2+fabs(epsilon+xt),2)/(2*t))- exp(-pow(-x0+a2+fabs(-xt+3*epsilon),2)/(2*t)))/pro);
    }
  }
}

//[[Rcpp::export()]]
List brown_crossing_sample(double x0,double xt, double t, double epsilon,
                           bool up=true,int nsize=1)
{
  if(x0< - epsilon){ up = false;}
  if(x0 > epsilon) { up = true;}
  // Rcout<<up<<endl;
  int nup;
  if(up==false){ nup=0;}
  else {nup=1;}
  double cumsum=0,prob;
  int counts=51;
  double r = R::runif(0,1);
  for(int i = 0; i<50; i++)
  {
    prob = crossing_brownbridge(x0,xt,t,i,epsilon, up);
    //  Rcout<<"i="<<i<<"\t"<<prob(i)<<endl;
    cumsum += prob;
    if(r<cumsum)
    {
      counts = i;
      break;
    }
    //    cumsums(i) = cumsum;
    //   sp(i) = i;
  }

  int upnext =(counts + nup)%2;
  counts = floor((counts + 1-nup)/2);
  List res;
  res["samples"] = counts;
  res["upnext"] = upnext;
  return(res);
}

//[[Rcpp::export()]]
int fills_in_cpp(arma::mat edge,arma::vec edgelength,
                 int Nnode,arma::vec data,double epsilon,
                 double thed=0,int cut=0,bool less=false)
{
  int N,p,cros=0;
  N=2*Nnode+1;
  p=Nnode+1;
  double x0 = data(p);
  arma::vec start=edge.col(0);
  arma::vec upd(N);
  upd(p) = (x0>=-epsilon);
  int i,j;
  for(i=1;i<=edge.n_rows;i++)
  {
    j = edge(i-1,1);
    List res = brown_crossing_sample(data(start(i-1)-1),data(j-1),edgelength(i-1),\
                                     epsilon,upd(start(i-1)-1),1);
    cros=cros + as<int>(res["samples"]);
    upd[j-1] = as<int>(res["upnext"]);
  }
  return(cros);
}


//[[Rcpp::export()]]
arma::vec fills_in_mat(arma::mat edge,arma::vec edgelength,
                       int Nnode,arma::mat data,double epsilon,
                       double thed=0,int cut=0,bool less=false)
{
  int N,p,cros;
  arma::vec res(data.n_rows);
  N=2*Nnode+1;
  p=Nnode+1;
  arma::vec start=edge.col(0);
  arma::vec upd(N);
  int i,j,k;
  for(k=0;k<data.n_rows;k++){
    cros=0;
    double x0 = data(k,p);
    upd(p) = (x0>=-epsilon);
    for(i=1;i<=edge.n_rows;i++)
    {
      j = edge(i-1,1);
      List res = brown_crossing_sample(data(k,start(i-1)-1),data(k,j-1),edgelength(i-1),\
                                       epsilon,upd(start(i-1)-1),1);
      cros=cros + as<int>(res["samples"]);
      upd[j-1] = as<int>(res["upnext"]);
    }
    res(k) = cros;
  }
  return(res);
}


//[[Rcpp::export()]]
arma::mat brown_tree_prior_node_cpp(arma::mat edge,arma::vec edgelength,int Nnode,
                                    arma::vec rootprior,int nsize)
{
  int N,p;
  N = 2 * Nnode + 1;
  p = Nnode + 1;
  double x0;
  arma::mat res(nsize,N);
  arma::vec start = edge.col(0);
  int i,j;
  for(int k=0;k<nsize;k++)
  {
    x0 = R::rnorm(rootprior(0),rootprior(1));
    res(k,p) = x0;
    for(i=1;i<=edge.n_rows;i++)
    {
      j = edge(i-1,1);
      res(k,j-1) = R::rnorm(res(k,start(i-1)-1),sqrt(edgelength(i-1)));
    }
  }
  return(res);
}

//[[Rcpp::export()]]
arma::vec trnorm0(double mu,double sigma,int d,int nsize=1)
{
  //sample from truncated normal distribution on [a,b]
  arma::vec res(nsize);
  double u,alpha=R::pnorm((0-mu)/sigma,0,1,1,0);

  for(int i=0;i<nsize;i++)
  {
    u = R::runif(0,1);
    if(d ==1) res(i) = R::qnorm(alpha + (1-alpha)*u ,0,1,1,0)*sigma+mu;
    else res(i) = R::qnorm(alpha*u,0,1,1,0)*sigma+mu;
  }
  return(res);
}



//[[Rcpp::export()]]
arma::mat getneighbor_cpp(arma::mat edge,arma::vec edgelength,int Nnode)
{
  int N,p;
  N = 2 * Nnode + 1;
  p = Nnode + 1;
  arma::mat nbor,nedge;
  nbor.zeros(3,N);
  nedge.zeros(3,N);
  for(int i=1;i<=N;i++)
  {
    int t=0;
    for(int j=0;j<edge.n_rows;j++)
    {
      if(edge(j,0)==i)
      {
        nbor(t,i-1) = edge(j,1);
        nedge(t,i-1) = edgelength(j);
        t++;
      }
      if(edge(j,1)==i)
      {
        nbor(t,i-1) = edge(j,0);
        nedge(t,i-1) = edgelength(j);
        t++;
      }
    }
  }
  arma::mat res(6,N);
  res.submat(0,0,2,N-1) = nbor;
  res.submat(3,0,5,N-1) = nedge;
  return(res);
}

//[[Rcpp::export()]]

arma::mat posterior_update(arma::mat edge, arma::vec edgelength,int Nnode,arma::vec d,
                           double thed,arma::vec rootprior, int ngen, int burnin, int thin = 1)
{
  Rcpp::checkUserInterrupt();
  int N,p,j;
  double v1,v2,v3,z1,z2,z3,mu,sigma2;
  //  arma::mat nbor,nedge;
  N = 2 * Nnode + 1;
  p = Nnode + 1;
  arma::vec th(p),batch;
  arma::mat nb = getneighbor_cpp(edge,edgelength,Nnode);
  arma::mat nbor = nb.submat(0,0,2,N-1);
  arma::mat nedge = nb.submat(3,0,5,N-1);
  int k = ceil((ngen-burnin)/thin);
  arma::mat mcmc(k,N);
  batch = vectorise(brown_tree_prior_node_cpp(edge,edgelength,Nnode,rootprior,1));
  for(int i=0;i<p;i++)
  {
    batch(i) = 0;
  }
  //update node
  int l=0;
  for(int i=0;i<ngen;i++)
  {
    // root node
    v1 = nedge(0,p);
    v2 = nedge(1,p);
    v3 = rootprior(1)*rootprior(1);
    z1 = batch(nbor(0,p)-1);
    z2 = batch(nbor(1,p)-1);
    z3 = rootprior(0);
    mu = (z1/v1+z2/v2+z3/v3)/(1/v1+1/v2+1/v3);
    sigma2 = 1/(1/v1+1/v2+1/v3);
    batch(p) = R::rnorm(mu,sqrt(sigma2));
    // internal nodes
    for(j = p+1;j<N;j++)
    {
      v1 = nedge(0,j);
      v2 = nedge(1,j);
      v3 = nedge(2,j);
      z1 = batch(nbor(0,j)-1);
      z2 = batch(nbor(1,j)-1);
      z3 = batch(nbor(2,j)-1);
      mu = (z1/v1+z2/v2+z3/v3)/(1/v1+1/v2+1/v3);
      sigma2 = 1/(1/v1+1/v2+1/v3);
      batch(j) = R::rnorm(mu,sqrt(sigma2));
    }
    //tip nodes
    for(j=0;j<p;j++)
    {
      mu=batch(nbor(0,j)-1);
      sigma2=nedge(0,j);
      batch(j) = trnorm0(mu,sqrt(sigma2),d(j),1)(0);
    }
    if(!(i<burnin||i%thin!=0))
    {
      mcmc.row(l)=batch.t();
      l++;
    }

  }
  return(mcmc);

}

