// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
using namespace Numer;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class optimFun: public MFuncGrad
{
private:
  const MapMat x;
  const MapMat Y;
  const MapVec size;
  const MapMat weights;
  const MapMat offset;
  double sigma2;
  double S2;
  const MapVec no_shrink;
  const MapVec shrink;
  const MapVec cnst;
  int i;
public:
  optimFun(const MapMat x_, const MapMat Y_, const MapVec size_, const MapMat weights_, const MapMat offset_, double sigma2_, double S2_, const MapVec no_shrink_, const MapVec shrink_, const MapVec cnst_, int i_) : x(x_), Y(Y_), size(size_), weights(weights_), offset(offset_), sigma2(sigma2_), S2(S2_), no_shrink(no_shrink_), shrink(shrink_), cnst(cnst_), i(i_) {}

  double f_grad(Constvec& beta, Refvec grad)
  {

    int B = x.cols(); // e.g. number of betas
    
    Eigen::ArrayXd xbeta = x * beta;
    Eigen::ArrayXd exp_xbeta_off = xbeta.exp().array() + offset.col(i).array();

    Eigen::ArrayXd a = Y.col(i).array() + size[i];
    Eigen::ArrayXd b = exp_xbeta_off + size[i];

    Eigen::VectorXd c = a * exp_xbeta_off * b.inverse() - Y.col(i).array();
    Eigen::VectorXd cw = c.array() * weights.col(i).array();
    
    Eigen::ArrayXd d = a * b.log() - Y.col(i).array() * xbeta;
    Eigen::ArrayXd dw = d * weights.col(i).array();

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

    // this is the negative log posterior plus a constant to keep it above 0
    const double f = dw.sum() + neg_prior + cnst[i];

    // this is the gradient of the negative log posterior
    Eigen::ArrayXd d_neg_lik = x.transpose() * cw;
    grad = d_neg_lik + d_neg_prior;
    
    return f;
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix nbinomGLM(Rcpp::NumericMatrix x, Rcpp::NumericMatrix Y,
			      Rcpp::NumericVector size, Rcpp::NumericMatrix weights,
			      Rcpp::NumericMatrix offset, double sigma2, double S2,
			      Rcpp::NumericVector no_shrink, Rcpp::NumericVector shrink,
			      Rcpp::NumericVector intercept, Rcpp::NumericVector cnst)
{

  // this optimization code is modeled after the RcppNumerical example:
  // https://github.com/yixuan/RcppNumerical/blob/master/README.md

  // here, x is the design matrix, Y is num samples x num genes
  // size is 1/dispersion
  // weights and offset are same dimension as Y
  // no_shrink and shrink are vectors of indices of which coefs to shrink (1-based)
  // intercept is the log of the basemean for initializing beta[0]
  // cnst is a constant to keep the negative log posterior above 0
  
  const MapMat mx = Rcpp::as<MapMat>(x);
  const MapMat mY = Rcpp::as<MapMat>(Y);
  const MapVec msize = Rcpp::as<MapVec>(size);
  const MapMat mweights = Rcpp::as<MapMat>(weights);
  const MapMat moffset = Rcpp::as<MapMat>(offset);
  const MapVec mno_shrink = Rcpp::as<MapVec>(no_shrink);
  const MapVec mshrink = Rcpp::as<MapVec>(shrink);
  const MapVec mintercept = Rcpp::as<MapVec>(intercept);
  const MapVec mcnst = Rcpp::as<MapVec>(cnst);

  int G = Y.ncol(); // e.g. number of genes
  int B = x.ncol(); // e.g. number of betas
  Eigen::MatrixXd betas(B, G);
  Eigen::VectorXd beta(B);

  double fopt;
  for (int i = 0; i < G; i++) {
    optimFun nll(mx, mY, msize, mweights, moffset, sigma2, S2, mno_shrink, mshrink, mcnst, i);
    beta.setZero();
    beta[0] = mintercept[i];
    int status = optim_lbfgs(nll, beta, fopt);
    betas.col(i) = beta;
  }

  return Rcpp::wrap(betas);
}
