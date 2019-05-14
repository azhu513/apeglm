// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
 
#include <RcppNumerical.h>
#include <Rcpp.h>
//required for the lgamma and digamma functions
#include <unsupported/Eigen/SpecialFunctions>
using namespace Numer;
using namespace Rcpp;

// betabinGLM is an extension of nbinomGLM to work for beta-binomial distribution
// written by Josh Zitovsky Spring semester 2019

typedef Eigen::Map<Eigen::MatrixXd> MapMat; 
typedef Eigen::Map<Eigen::VectorXd> MapVec; 

class optimFunBeta: public MFuncGrad
{
private:
  MapMat x;
  MapMat Y;
  MapMat sizes;
  MapVec thetas;
  MapMat weights;
  double sigma2;
  double S2;
  MapVec no_shrink;
  MapVec shrink;
  const MapVec cnst;
  int i;
  double lbd;
  double ubd;
public:
  optimFunBeta(MapMat x_, MapMat Y_, MapMat sizes_, MapVec thetas_, MapMat weights_, double sigma2_, double S2_, MapVec no_shrink_, MapVec shrink_, const MapVec cnst_, int i_, double lbd_, double ubd_) : x(x_), Y(Y_), sizes(sizes_), thetas(thetas_), weights(weights_), sigma2(sigma2_), S2(S2_), no_shrink(no_shrink_), shrink(shrink_), cnst(cnst_), i(i_) ,lbd(lbd_), ubd(ubd_) {}

  double f_grad(Constvec& beta, Refvec grad)
  {

    int B = x.cols(); // e.g. number of betas
    
    Eigen::ArrayXd xbeta = x * beta;
    Eigen::ArrayXd exp_neg_xbeta = xbeta.exp().inverse();
    exp_neg_xbeta = exp_neg_xbeta.max(lbd).min(ubd);
    Eigen::ArrayXd prob = (1+exp_neg_xbeta).inverse();
    double theta = thetas[i];
    Eigen::ArrayXd thetap = theta*prob;
    Eigen::ArrayXd size = sizes.col(i);
    Eigen::ArrayXd y = Y.col(i);
    Eigen::ArrayXd weight = weights.col(i);

    Eigen::ArrayXd e1 = size - y + theta - thetap;
    Eigen::ArrayXd e2 = y + thetap;
    Eigen::ArrayXd e3 = theta - thetap;
    Eigen::ArrayXd lik_parts = e1.lgamma() + 
                               e2.lgamma() -
                               e3.lgamma() -
                               thetap.lgamma();
    Eigen::ArrayXd weighted_parts = lik_parts * weight;
    double neg_obs = -1 * weighted_parts.sum(); //beta-binomial negative log-likelihood
  
    Eigen::ArrayXd e4 = theta*exp_neg_xbeta/pow(1+exp_neg_xbeta, 2);
    Eigen::VectorXd e5 = e4*(-e1.digamma()+e2.digamma()+e3.digamma()-thetap.digamma())*weight;
    Eigen::ArrayXd d_neg_obs = -1 * x.transpose() * e5; //beta-binomial gradient

    double neg_prior = 0.0;
    Eigen::ArrayXd d_neg_prior(B);

    // the coefficients to "not shrink" (very wide Normal prior)
    for (int j = 0; j < no_shrink.size(); j++) {
      int k = no_shrink[j] - 1; // change to 0-based indexing
      neg_prior += pow(beta[k], 2.0)/(2.0 * sigma2);
      d_neg_prior[k] = beta[k]/sigma2;
    }

    // the coefficients to shrink (Cauchy prior)
    for (int j = 0; j < shrink.size(); j++) {
      int k = shrink[j] - 1; // change to 0-based indexing
      neg_prior += log1p(pow(beta[k], 2.0)/S2);
      d_neg_prior[k] = 2.0 * beta[k] / (S2 + pow(beta[k], 2.0));
    }

    // this is the negative log posterior, scaled,
    // plus a constant to keep the values somewhat small (and not too close to zero)
    double cn = cnst[i];
    const double f = (neg_obs + neg_prior)/cn + 10.0;

    // this is the gradient of the negative log posterior, scaled
    grad = (d_neg_obs + d_neg_prior)/cn;
    
    return f;
  }
};

// [[Rcpp::export]]
Rcpp::List betabinGLM(Rcpp::NumericMatrix x, Rcpp::NumericMatrix Y,
		      Rcpp::NumericMatrix sizes, Rcpp::NumericVector thetas, 
		      Rcpp::NumericMatrix weights,double sigma2, double S2,
		      Rcpp::NumericVector no_shrink, Rcpp::NumericVector shrink,
		      Rcpp::NumericVector init, Rcpp::NumericVector cnst, double tol, double lbd, double ubd)
{

  // this optimization code is modeled after the RcppNumerical example:
  // https://github.com/yixuan/RcppNumerical/blob/master/README.md

  // see nbinomGLM for details on variables
  
  const MapMat mx = Rcpp::as<MapMat>(x);
  const MapMat mY = Rcpp::as<MapMat>(Y);
  const MapMat msizes = Rcpp::as<MapMat>(sizes);
  const MapVec mthetas = Rcpp::as<MapVec>(thetas);
  const MapMat mweights = Rcpp::as<MapMat>(weights);
  const MapVec mno_shrink = Rcpp::as<MapVec>(no_shrink);
  const MapVec mshrink = Rcpp::as<MapVec>(shrink);
  const MapVec minit = Rcpp::as<MapVec>(init);
  const MapVec mcnst = Rcpp::as<MapVec>(cnst);

  int G = Y.ncol(); // e.g. number of genes
  int B = x.ncol(); // e.g. number of betas

  Eigen::MatrixXd betas(B, G);
  Eigen::VectorXd beta(B);

  Rcpp::NumericVector value(G);
  Rcpp::IntegerVector convergence(G);

  double fopt;
  for (int i = 0; i < G; i++) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    optimFunBeta nll(mx, mY, msizes, mthetas, mweights, sigma2, S2, mno_shrink, mshrink, mcnst, i, lbd, ubd);
    beta = minit;
    int status = optim_lbfgs(nll, beta, fopt, 300, tol, 1e-8);
    betas.col(i) = beta;
    value[i] = fopt;
    convergence[i] = status;
  }

  return Rcpp::List::create(Rcpp::Named("betas") = betas,
			    Rcpp::Named("value") = value,
			    Rcpp::Named("convergence") = convergence);
}
