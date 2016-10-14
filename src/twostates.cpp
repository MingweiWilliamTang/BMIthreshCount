#include<RcppArmadillo.h>
#include<Rcpp.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include<R.h>
#include<math.h>
using namespace Rcpp;
using namespace arma;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
arma::mat twoStateRateMatrix(double lambda_01, double lambda_10) {

  arma::mat rate_matrix = arma::mat(2,2);

  rate_matrix(0,0) = -lambda_01;
  rate_matrix(0,1) =  lambda_01;
  rate_matrix(1,0) =  lambda_10;
  rate_matrix(1,1) = -lambda_10;

  return rate_matrix;
}

arma::mat twoStateTransProb(double lambda_01, double lambda_10, double time) {

  double total_rate = lambda_01 + lambda_10;

  arma::mat prob_matrix = arma::mat(2,2);

  prob_matrix(0,0) = (lambda_10 + lambda_01*exp(-total_rate*time))/total_rate;
  prob_matrix(0,1) = (lambda_01 - lambda_01*exp(-total_rate*time))/total_rate;
  prob_matrix(1,0) = (lambda_10 - lambda_10*exp(-total_rate*time))/total_rate;
  prob_matrix(1,1) = (lambda_01 + lambda_10*exp(-total_rate*time))/total_rate;

  return prob_matrix;
}

int sampleOnce(arma::colvec weights, double rUnif) {

  double total = sum(weights);
  double cumProb=0;

  int i;

  for (i = 0; i < weights.n_elem; i++) {
    cumProb += weights(i) / total;
    if (rUnif < cumProb) break;
  }

  return(i);
}

// Samples 2-state continuous-time Markov chain conditional on starting and ending states
// Saves only sufficient statistics: number of 0->1 jumps, number of 1->0 jumps and dwell times in 0,1

NumericVector twoStateUnifSample(arma::mat& rateMatrix, int startState, int endState, double elapsedTime, double transProb)
{
  int numStates = rateMatrix.n_cols;

  // Set the rate of the dominating Poisson process
  double poissonRate = -1.0*arma::min(rateMatrix.diag());
  arma::mat dominTransProb(numStates,numStates);
  dominTransProb.eye();
  arma::cube powerTransProb = arma::cube(numStates, numStates,1);
  powerTransProb.slice(0) = dominTransProb;
  dominTransProb = dominTransProb + rateMatrix/poissonRate;


  // simulate the number of jumps, including fictitious ones

  double rU = as<double>(runif(1));
  double cum = 0;

  if (startState==endState){
    cum = ::Rf_dpois(0,poissonRate*elapsedTime,0)/transProb;
  }

  bool notExceed = true;

  if (cum > rU){
    notExceed = false;
  }

  int numJumps = 0;
  double nextProb;

  while(notExceed){

    numJumps++;
    powerTransProb.insert_slices(numJumps,1,false);
    powerTransProb.slice(numJumps) = powerTransProb.slice(numJumps-1)*dominTransProb;
    nextProb = ::Rf_dpois(numJumps,poissonRate*elapsedTime,0)*powerTransProb(startState,endState,numJumps)/transProb;
    cum += nextProb;

    if (cum > rU){
      notExceed = false;
    }
  }

  NumericVector sufStat(4);

  // if numJumps = 0: done
  if (numJumps == 0 || ((numJumps == 1) && (startState == endState))){

    sufStat[2+startState] = elapsedTime;

  }else{ // if one true jump: done
    if ((numJumps == 1) && (startState != endState)){

      sufStat[startState] = 1;
      sufStat[2+startState] = elapsedTime*(as<double>(runif(1)));
      sufStat[2+endState] = elapsedTime - sufStat[2+startState];

    }else{ // Case (nJmp >= 2)

      IntegerVector stateSpace = seq_len(numStates) - 1;

      // Simulate jumping times
      NumericVector dominJumpTimes = elapsedTime*runif(numJumps);
      std::sort(dominJumpTimes.begin(),dominJumpTimes.end());
      NumericVector dominStates(numJumps+1);
      dominStates[0] = startState;
      dominStates[numJumps] = endState;

      // Simulate states of the dominating Markov chain and records sufficient statistics
      for (int i = 0; i < numJumps; i++){
        dominStates[i+1] = sampleOnce(powerTransProb.slice(1).row(dominStates[i]).t()%powerTransProb.slice(numJumps-i-1).col(endState), as<double>(runif(1)));

        if (i == 0){
          sufStat[2+startState] = dominJumpTimes[0];
        }else{
          sufStat[2+dominStates[i]] += dominJumpTimes[i] - dominJumpTimes[i-1];
        }

        if (dominStates[i+1] != dominStates[i]){
          sufStat[dominStates[i]] += 1;
        }
      }

      sufStat[2+dominStates[numJumps-1]] += (elapsedTime - dominJumpTimes[numJumps-1]);
    }
  }

  return sufStat;
}

/* Arguments:
treeEdges: edge matrix of the ape object phylo
numIntNodes: number of internal nodes (not necessary, but convinient to get it from phylo)
tipStates: integer vector of tip states (-1=missing value)
cubeProbMat: array of probability matrices for each edge of the tree

Two important assumptions:
1. edges in the the edge matrix and probability matrices are in the "pruningwise" order;
see ?reorder.phylo for more details
2. tip state vector is ordered according to the tip numbering in the edge matrix
*/

arma::mat PartLikelihoods(const arma::Mat<int> & treeEdges, const IntegerVector & tipStates,
                          const arma::cube & cubeProbMat)
{
  /// get number of edges
  int numEdges = treeEdges.n_rows;

  // get number of internal notes
  int numIntNodes = numEdges/2;

  // get number of tips in the tree
  int numTips = tipStates.size();

  // prepare a matrix for storing regular (backward) partial likelihoods
  arma::mat partialLike = arma::zeros<arma::mat>(numTips + numIntNodes, 2);

  for (int i = 0; i < numTips; i++){
    if (tipStates(i) == -1) {// -1 denotes a missing value
      partialLike.row(i) = arma::ones<arma::rowvec>(2);
    } else {
      partialLike(i, tipStates(i)) = 1.0;
    }
  }

  // compute regular partial likelihoods for all internal nodes
  for (int i = 0; i < numEdges; i+=2){
    // parent1 = treeEdges[i,0] or treeEdges[i+1,0] also treeEdges indices should be shifted down by one
    partialLike.row(treeEdges(i,0)-1) = (partialLike.row(treeEdges(i,1)-1)*
      cubeProbMat.slice(i).t())%(partialLike.row(treeEdges(i+1,1)-1)*cubeProbMat.slice(i+1).t());
  }

  return partialLike;
}
// [[Rcpp::export]]
double TwoStatePhyloLikelihood1(arma::Mat<int>& treeEdges, IntegerVector& tipStates,
                                NumericVector& branchLengths, double lambda_01, double lambda_10,
                                NumericVector& rootDist) {

  // convert to rootDist to arma:rowvec
  arma::colvec armaRootDist = as<arma::colvec>(rootDist);

  // get number of edges
  int numEdges = treeEdges.n_rows;

  // get number of tips
  int numTips = tipStates.size();

  // Compute transition probabilities for each branch on the tree and store in arma:cube
  arma::cube cubeProbMat(2,2,numEdges);

  for (int i=0; i < numEdges; i++){
    cubeProbMat.slice(i) = twoStateTransProb(lambda_01, lambda_10, branchLengths(i));
  }

  // Compute partial likelihoods at all internal nodes
  arma::mat partLike = PartLikelihoods(treeEdges, tipStates, cubeProbMat);

  return log(sum(partLike.row(numTips).t()%armaRootDist));
}

// [[Rcpp::export]]
double TwoStatePhyloLikelihood2(const arma::Mat<int> & treeEdges, const IntegerVector & tipStates,
                                const arma::vec & branchLengths, const double & lambda_01,
                                const double & lambda_10, const arma::vec & armaRootDist)
{
  // get number of edges
  int numEdges = treeEdges.n_rows;

  // get number of tips
  int numTips = tipStates.size();

  // Compute transition probabilities for each branch on the tree and store in arma:cube
  arma::cube cubeProbMat(2, 2, numEdges);

  for (int i=0; i < numEdges; i++){
    cubeProbMat.slice(i) = twoStateTransProb(lambda_01, lambda_10, branchLengths(i));
  }

  // Compute partial likelihoods at all internal nodes
  arma::mat partLike = PartLikelihoods(treeEdges, tipStates, cubeProbMat);

  return sum(partLike.row(numTips).t()%armaRootDist);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
