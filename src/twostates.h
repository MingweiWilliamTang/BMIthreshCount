#ifndef BROWN_CROSSING_H_
#define BROWN_CROSSING_H_
#include<RcppArmadilloExtensions/sample.h>
#include<RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include<Rcpp.h>
#include<R.h>
#include<math.h>
using namespace Rcpp;
using namespace arma;

double crossing_brownbridge(double x0,double xt, double t ,int ncross,
                            double epsilon ,bool up ,bool counditional);

List brown_crossing_sample(double x0,double xt, double t, double epsilon,
                           bool up,int nsize);

int fills_in_cpp(arma::mat edge,arma::vec edgelength,
                 int Nnode,arma::vec data,double epsilon,
                 double thed,int cut,bool less);

arma::vec fills_in_mat(arma::mat edge,arma::vec edgelength,
                       int Nnode,arma::mat data,double epsilon,
                       double thed,int cut,bool less);

arma::mat brown_tree_prior_node_cpp(arma::mat edge,arma::vec edgelength,int Nnode,
                                    arma::vec rootprior,int nsize);

arma::mat getneighbor_cpp(arma::mat edge,arma::vec edgelength,int Nnode);


arma::mat posterior_update(arma::mat edge, arma::vec edgelength,int Nnode,arma::vec d,
                           double thed,arma::vec rootprior, int ngen, int burnin, int thin);


#endif
