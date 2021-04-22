#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::plugins("cpp11")]]

// Define some constants

#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

// Define helper functions
namespace help{

// The following codes are from PGBVS.cpp by Dr. Matt Koslovsky (https://github.com/mkoslovsky/PGBVS)

// Generate exponential distribution random variates
double exprnd(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }
  return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {
    out = mu*mu / out;
  }
  return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = help::exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss(double z, double t)
{
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / help::truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = help::randinvg(mu);
    }
  }
  return X;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q);
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1)
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + help::exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = help::tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = help::aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1)
    {
      Sn = Sn + asgn * help::aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// Simulate MVT normal data
arma::mat mvrnormArma( int n, arma::vec mu, arma::mat sigma ) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn( n, ncols );
  return arma::repmat( mu, 1, n ).t()  + Y * arma::chol( sigma );
}

// Code above are from PGBVS.cpp





























// Function :: propose new gamma
arma::vec new_gamma(arma::vec gamma, int k){
  
  int tmp=0;
  if (gamma(k)==1){
    tmp=0;
  }else if(gamma(k)==0){
    tmp=1;
  }
  gamma(k)=tmp;
  
  return gamma;
}

// Function :: propose new beta after propose new_gamma
arma::vec new_beta( arma::vec beta, arma::vec new_gamma, int k, double v2){
  
  double tmp=0;
  if (new_gamma(k)==1){
    tmp = beta(k) + rnorm( 1,0,sqrt(v2))[ 0 ];
    beta(k)=tmp;
  }else if(new_gamma(k)==0){
    tmp=0;
    beta(k)=tmp;
  }
  return beta;
}

// Function :: Return the difference between the log of prior using old and new beta
double update_log_beta(int k, arma::vec beta, double v1){
  double diff = -0.50*log( 2*atan(1)*4*v1 ) - 1/( 2*v1 )*pow( beta( k ), 2 );
  return diff;
}

// Function :: Return the difference between old and new probability of log gamma
double update_log_gamma(int k, arma::vec gamma, double a, double b){
  double g = gamma(k);
  double diff = lgamma( g+a ) + lgamma( 1-g+b ) + lgamma(a+b) - lgamma( a ) - lgamma( b ) - lgamma( 1+a+b );
  return diff;
}

// Function :: Call sample function from R. Sample one number from a vector.
int sample_cpp( IntegerVector x ){
  // Calling sample()
  Function f( "sample" );
  IntegerVector sampled = f( x, Named( "size" ) = 1 );
  return sampled[0];
}

// Function :: Initial probability
double iprob(int h1,int y1){
  double pi1h1 = 0;
  if (h1==y1){
    pi1h1=1;
  }
  return pi1h1;
}

// Function :: Emission probability
double eprob(int y, int h, double p0, double p1){
  double fyh = 0;
  if (h==1){
    fyh = pow(p1,y)*pow((1-p1),(1-y));
  }else if(h==0){
    fyh = pow(p0,(1-y))*pow((1-p0),y);
  }
  return fyh;
}

// Function :: Transition probability
double tprob(int hj_prev, int hj_curr, double lambda, double mu, double delta){
  
  if (delta == NA){
    std::cout << "error!" << std::endl;
  }
  
  double q = 0;
  double pns = 0;
  double psn = 0;
  double pnn = 0;
  double pss = 0;
  
  pns = lambda/(lambda+mu)*(1-exp(-(lambda+mu)*delta));
  psn = mu/(lambda+mu)*(1-exp(-(lambda+mu)*delta));
  pss = 1 - psn;
  pnn = 1 - pns;
  
  if (hj_prev==1){
    if (hj_curr==1){
      q=pss;
    }else if(hj_curr==0){
      q=psn;
    }
  }else if (hj_prev==0){
    if (hj_curr==1){
      q=pns;
    }else if(hj_curr==0){
      q=pnn;
    }
  }
  return q;
}

// Function :: Scaled Forward-backward algorithm
arma::vec fb(arma::mat datmat, arma::vec beta, int pt, int pe, int subject_index, arma::vec sp, arma::vec nobs){
  
  int spi = sp[subject_index];
  int ni = nobs[subject_index];
  
  arma::mat xti = datmat(span(spi,spi+ni-1),span(1,pt));
  arma::mat xei = datmat(span(spi,spi+ni-1),span(pt+1,pt+pe));
  arma::vec y = datmat(span(spi,spi+ni-1),span(pt+pe+1));
  arma::vec deltai = datmat(span(spi,spi+ni-1),pt+pe+2);
  
  arma::mat fwd( 2, ni, fill::zeros ); // forward matrix 2 states * ni times
  
  arma::vec A(2, fill::zeros); // To store temporary variable for two states
  double b = 0; // to store the max between two states
  arma::vec bvec = A; // to make b a vector
  
  arma::vec beta_lambda = beta(span(0,pt-1));
  arma::vec beta_mu = beta(span(pt,2*pt-1));
  arma::vec beta_0 = beta(span(2*pt,2*pt+pe-1));
  arma::vec beta_1 = beta(span(2*pt+pe,2*pt+2*pe-1));
  
  // Forward
  
  // Initialization
  
  double p0 = 1/(1 + 1/exp( xei.row(0) * beta_0 )[0]);
  double p1 = 1/(1 + 1/exp( xei.row(0) * beta_1 )[0]);
  
  fwd(0,0) = log(help::iprob(0,y[0])) + log(help::eprob(y[0],0,p0,p1));
  fwd(1,0) = log(help::iprob(1,y[0])) + log(help::eprob(y[0],1,p0,p1));
  
  // Iterations
  
  for( int t = 1; t < ni; ++t ){
    
    double p0 = 1/(1 + 1/exp( xei.row(t) * beta_0 )[0]);
    double p1 = 1/(1 + 1/exp( xei.row(t) * beta_1 )[0]);
    
    double lambda = exp( xti.row(t-1) * beta_lambda)[0];
    double mu = exp( xti.row(t-1) * beta_mu)[0];
    
    // transition from t-1 to t
    
    for (int j=0; j<2; j++){
      for (int i=0; i<2; i++){
        A[i] = fwd(i,t-1) + log(help::tprob(i,j,lambda,mu,deltai[t])) + log(help::eprob(y[t],j,p0,p1));
      }
      b = max(A);
      bvec.fill(b);
      fwd(j,t) = b + log(sum(exp(A - bvec)));
      if (Rcpp::traits::is_nan<REALSXP>(fwd(j,t))){
        fwd(j,t) = -999;
      }
    }
    
  } // end for t
  
  // Backward
  
  arma::mat bwd(2,ni,fill::zeros);
  
  // Initialization
  
  bwd(0,ni-1) = 0;
  bwd(1,ni-1) = 0;
  
  // Iterations
  
  for (int t = ni-2; t>-1; t--){
    // transition from t to t+1
    
    double p0 = 1/(1 + 1/exp( xei.row(t+1) * beta_0 )[0]);
    double p1 = 1/(1 + 1/exp( xei.row(t+1) * beta_1 )[0]);
    
    double lambda = exp( xti.row(t) * beta_lambda)[0];
    double mu = exp( xti.row(t) * beta_mu)[0];
    
    
    for (int i=0; i<2; i++){
      for (int j=0; j<2; j++){
        A[j] = bwd(j,t+1) + log(help::tprob(i,j,lambda,mu,deltai[t+1])) + log(help::eprob(y[t+1],j,p0,p1));
      }
      b = max(A);
      bvec.fill(b);
      bwd(i,t) = b + log(sum(exp(A-bvec)));
      if (Rcpp::traits::is_nan<REALSXP>(bwd(i,t))){
        bwd(i,t) = -999;
      }
    }
    
  } // end for t
  
  arma::mat gamma(2,ni,fill::zeros);
  gamma = fwd + bwd;
  
  arma::vec h(ni, fill::zeros);
  for (int t=0;t<ni;t++){
    double ga = exp(gamma(1,t))/(exp(gamma(1,t))+exp(gamma(0,t)));
    h[t] = R::rbinom(1,ga);
  }
  
  return h;
}

// Function :: simulate the whole hidden chain for all subjects
arma::vec simulate_hidden(arma::mat datmat, arma::vec beta, int pt, int pe, arma::vec sp, arma::vec nobs){
  arma::vec h(datmat.n_rows);
  int n = sp.size(); // number of subjects
  for (int i=0; i<n; i++){
    arma::vec hi = fb(datmat,beta,pt,pe,i,sp,nobs);
    h(span(sp[i],sp[i]+nobs[i]-1)) = hi;
  }
  return h;
}

// Function :: Log likelihood
double log_lik_cpp(arma::mat datmat, arma::vec h, arma::vec beta, int pt, int pe, arma::vec nobs, arma::vec sp, int k){
  double loglik = 0;
  int m = sp.size();
  int rc = 0;
  
  arma::vec beta_lambda = beta(span(0,pt-1));
  arma::vec beta_mu = beta(span(pt,2*pt-1));
  arma::vec beta_0 = beta(span(2*pt,2*pt+pe-1));
  arma::vec beta_1 = beta(span(2*pt+pe,2*pt+2*pe-1));
  
  for (int i=0; i<m; i++){
    int ni = nobs[i];
    for (int j=0; j<ni; j++){
      
      // transition part of the log-likelihood
      if ( k > (2*pt-1) ){
        double p00 = exp(datmat(rc,span(pt+1,pt+pe))*beta_0)[0];
        double p10 = exp(datmat(rc,span(pt+1,pt+pe))*beta_1)[0];
        double p0 = p00/(1+p00);
        double p1 = p10/(1+p10);
        loglik += log(eprob(datmat(rc,pt+pe+1),h[rc],p0,p1)); // emission prob
      }
      
      // emission part of the log-likelihood
      if ( j > 0 && k < 2*pt ){
        double lambda = exp(datmat(rc-1,span(1,pt))*beta_lambda)[0];
        double mu = exp(datmat(rc-1,span(1,pt))*beta_mu)[0];
        loglik += log(tprob(h[rc-1],h[rc],lambda,mu,datmat(rc,pt+pe+2))); // transition prob
      }
      
      rc += 1;
      
    }
  }
  return loglik;
}

// Function :: sampling posterior coefficients assuming prior beta~N(0,diag(v1)), given beta from last step
arma::vec sample_beta_PG(arma::mat x, arma::vec y, arma::vec beta, double v1, int nrep){
  int p = x.n_cols; // dimension of the vector of regressors
  int N = x.n_rows; // number of observations
  arma::vec lil_omega_vec(N);
  double z = 0;
  
  for (int r=0;r<nrep;r++){
    for (int k=0;k<N;k++){
      z = ( x.row(k) * beta ).eval()(0,0);
      lil_omega_vec[k] = samplepg(z);
    }
    
    arma::mat omega_mat(N,N); // big omega diagonal matrix
    omega_mat.zeros();
    for (int k=0;k<N;k++){
      omega_mat(k,k)=lil_omega_vec[k];
    }
    
    arma::mat v_omega(p,p);
    v_omega.zeros();
    v_omega = inv( x.t() * omega_mat * x + (1/v1) * eye(p,p));
    
    arma::vec kappa(N);
    kappa.zeros();
    kappa = y - 0.5;
    
    arma::vec m_omega(p);
    m_omega.zeros();
    m_omega = v_omega * ( x.t() * kappa );
    
    beta = mvrnormArma(1,m_omega,v_omega).row(0).t();
  }
  
  return beta;
}

// Function :: polya-gamma update
arma::vec propose_beta_PG(arma::vec h, arma::vec beta, arma::vec gamma, arma::mat datmat, int pt, int pe, int nu, double v1, int nrep){
  
  arma::uvec index0 = find(h==0);
  arma::uvec index1 = find(h==1);
  arma::vec gamma0 = gamma(span(2*pt,2*pt+pe-1));
  arma::uvec indexg0 = find(gamma0==1);
  arma::vec gamma1 = gamma(span(2*pt+pe,2*pt+2*pe-1));
  arma::uvec indexg1 = find(gamma1==1);
  arma::mat xe = datmat.cols(pt+1,pt+pe);
  arma::vec y = datmat.col(pt+pe+1);
  
  arma::vec beta_lambda = beta(span(0,pt-1));
  arma::vec beta_mu = beta(span(pt,2*pt-1));
  arma::vec beta_0 = beta(span(2*pt,2*pt+pe-1));
  arma::vec beta_1 = beta(span(2*pt+pe,2*pt+2*pe-1));
  
  int nj = h.size();
  arma::vec tellTruth(nj);
  tellTruth.zeros();
  for (int tmp=0; tmp<nj; tmp++){
    if (h[tmp]==y[tmp]){
      tellTruth[tmp]=1;
    }
  }
  
  if (nu==0){
    arma::vec beta_nu_nonzero = sample_beta_PG(xe(index0,indexg0), tellTruth(index0), beta_0(indexg0), v1, nrep);
    beta_0(indexg0) = beta_nu_nonzero;
    int betatmp = 2*pt; // starting point of the coefficients of interest
    beta(span(betatmp,betatmp+pe-1)) = beta_0;
  }else if (nu==1){
    arma::vec beta_nu_nonzero = sample_beta_PG(xe(index1,indexg1), tellTruth(index1), beta_1(indexg1), v1, nrep);
    beta_1(indexg1) = beta_nu_nonzero;
    int betatmp = 2*pt+pe; // starting point of the coefficients of interest
    beta(span(betatmp,betatmp+pe-1)) = beta_1;
  }
  return beta;
}

// between step for the kth covariate
List between_step(int k,
                  arma::mat datmat,
                  arma::vec h,
                  arma::vec beta,
                  arma::vec gamma,
                  arma::vec nobs,
                  arma::vec sp,
                  double v1,
                  double v2,
                  int pt,
                  int pe,
                  double a,
                  double b
){
  
  bool add = (gamma(k)==0);
  
  // propose gamma and beta
  arma::vec gamma_proposal = new_gamma(gamma, k);
  arma::vec beta_proposal = new_beta(beta, gamma_proposal, k, v2);
  
  // calculate ratio
  double r = 0;
  double loglik = log_lik_cpp(datmat, h, beta, pt, pe, nobs, sp, k);
  double loglik_proposal = log_lik_cpp(datmat, h, beta_proposal, pt, pe, nobs, sp, k);
  
  r = loglik_proposal - loglik  +
    update_log_gamma(k, gamma_proposal, a, b) - update_log_gamma(k, gamma, a, b);
  if (add){
    r = r + update_log_beta(k, beta_proposal, v1);
  }else{
    r = r - update_log_beta(k, beta, v1);
  }
  
  // Calculate acceptance probability
  double tmp = log(runif(1)[0]);
  
  // Determine acceptance
  if(tmp < r){
    beta = beta_proposal;
    gamma = gamma_proposal;
  }
  
  // Return output
  List between_return( 2 );
  between_return[ 0 ] = beta;
  between_return[ 1 ] = gamma;
  return between_return;
}

// Function :: within step
arma::vec within_step(bool hmm,
                      arma::mat datmat,
                      arma::vec h,
                      arma::vec beta,
                      arma::vec gamma,
                      arma::vec nobs,
                      arma::vec sp,
                      double v1,
                      double v2,
                      int pt,
                      int pe
){
  
  // update transition parameters
  for (int k=0; k<2*pt; ++k){
    if (gamma(k) == 1){
      arma::vec beta_proposal = new_beta(beta, gamma, k, v2);
      // calculate ratio
      double r = 0;
      r = log_lik_cpp(datmat, h, beta_proposal, pt, pe, nobs, sp, k)  + update_log_beta(k, beta_proposal, v1) -
        log_lik_cpp(datmat, h, beta, pt, pe, nobs, sp, k) - update_log_beta(k, beta, v1);
      // Calculate acceptance probability
      double tmp = log(runif(1)[0]);
      // Determine acceptance
      if(tmp < r){
        beta = beta_proposal;
      }
    }
  }
  
  // propose new betas for emission
  if (hmm){
    beta = propose_beta_PG(h, beta, gamma, datmat, pt, pe, 0, v1, 1);
    beta = propose_beta_PG(h, beta, gamma, datmat, pt, pe, 1, v1, 1);
  } // if not hmm, avoid this part of the calculation for the emission covariates
  
  // Return output
  return beta;
}

} // end help namespace

// Function :: The MCMC algorithm
// [[Rcpp::export]]
List HMMbvs(
    arma::mat datmat,
    arma::vec beta,
    arma::vec gamma,
    arma::vec nobs,
    arma::vec sp,
    IntegerVector klist,
    bool hmm=true,
    int iterations=5000,
    double v1=5,
    double v2=1,
    int pt=30,
    int pe=20,
    double a=1,
    double b=9,
    int thin=10,
    int thin_hidden=10,
    int base01larger=-1,
    int base00larger=-1
){
  
  // INITIALIZATION
  
  int hidden_int = thin*thin_hidden; // save the hidden states every (hidden_int) steps
  
  arma::mat beta_list(iterations/thin, beta.size(), fill::zeros); // initialize space to store beta for each iteration (50000*104 is around 38Mb)
  arma::mat gamma_list(iterations/thin, gamma.size(), fill::zeros); // initialize space to store gamma for each iteration (50000*104)
  int ncol = datmat.n_rows; // number of observations in total (4650)
  arma::vec h(ncol, fill::zeros); // initialize a temporary space to store the hidden states
  
  // initialize space to store hidden states.
  // Every hidden_int=25 iterations we save one set of simulated hidden states after the between and within step
  arma::mat h_mat(iterations/hidden_int, ncol, fill::zeros);
  // (original dataset 50000*4650 is around 1.8 Gb, need to thin it down to 2000*4650, which is around 68 Mb)
  int bcount = 0; // count number of beta/gamma saved
  int hcount = 0; // count number of hidden states saved
  arma::vec loglik_list(iterations/thin, fill::zeros); // initialize space to store the total log likelihood for each iteration
  arma::vec obs=datmat.col(pt+pe+1);
  
  // std::cout << "good" << std::endl;
  
  for (int iter=0; iter<iterations; ++iter){
    
    if(hmm){
      h = help::simulate_hidden(datmat, beta, pt, pe, sp, nobs);
    } else {
      h = obs;
    }
    
    if (iter%hidden_int==0){
      h_mat.row(hcount) = h.t();
      hcount += 1;
    }
    
    ///////////////////
    // between steps //
    ///////////////////
    
    
    // transition indices: 0 (intercept), 1, ..., pt-1, pt (intercept), pt+1, ..., 2*pt-1
    // emission indices: 2*pt (intercept), 2*pt+1, ..., 2*pt+pe-1, 2*pt+pe (intercept), 2*pt+pe+1, ..., 2*pt+2*pe-1
    
    // sample k
    int k = help::sample_cpp(klist);
    List between_return = help::between_step(k, datmat, h, beta, gamma, nobs, sp, v1, v2, pt, pe, a, b);
    beta = as<arma::vec>(between_return[0]);
    gamma = as<arma::vec>(between_return[1]);
    
    
    //////////////////
    // within steps //
    //////////////////
    
    beta = help::within_step(hmm, datmat, h, beta, gamma, nobs, sp, v1, v2, pt, pe);
    
    double tmp_con=0;
    
    if (base01larger == 1){
      tmp_con = beta[0];
      if(tmp_con < beta[pt]){
        beta[0] = beta[pt];
        beta[pt] = tmp_con;
      }
    }else if(base01larger == 0){
      tmp_con = beta[0];
      if(tmp_con > beta[pt]){
        beta[0] = beta[pt];
        beta[pt] = tmp_con;
      }
    }
    
    if (base00larger == 1){
      tmp_con = beta[2*pt];
      if(tmp_con < beta[2*pt+pe]){
        beta[2*pt] = beta[2*pt+pe];
        beta[2*pt+pe] = tmp_con;
      }
    }else if(base00larger == 0){
      tmp_con = beta[2*pt];
      if(tmp_con > beta[2*pt+pe]){
        beta[2*pt] = beta[2*pt+pe];
        beta[2*pt+pe] = tmp_con;
      }
    }
      

    
    double loglik = 0;
    
    if (iter%thin==0){
      beta_list.row(bcount) = beta.t();
      gamma_list.row(bcount) = gamma.t();
      loglik = help::log_lik_cpp(datmat, h, beta, pt, pe, nobs, sp, pt); // transition log-lik
      if(hmm){
        loglik += help::log_lik_cpp(datmat, h, beta, pt, pe, nobs, sp, 2*pt); // emission log-lik, hmm only
      }
      loglik_list[bcount] = loglik;
      bcount += 1;
    }
    
    
    double printer = iter % 500;
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
    
    
    
  } // end iter
  
  // output
  List output( 4 );
  output[ 0 ] = beta_list;
  output[ 1 ] = gamma_list;
  output[ 2 ] = h_mat;
  output[ 3 ] = loglik_list;
  
  return output;
}
