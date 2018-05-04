// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;

vec invlogit(vec x){
  vec y = exp(x)/(1 + exp(x));
  return(y);
}


// [[Rcpp::export]]
double log_post(vec b, vec Y, mat X, vec b_s, double jump_v) {
  int n = X.n_rows;
  int p = X.n_cols;
  vec xb;
  vec pi;
  vec mu;
  mat sd_prior(p, p);
  mat sd(p, p);
  double likelihood = 0;
  double prior = 0;
  double proposal;
  double log_posterior = 0;
  
  // linear predictor
  xb =  X*b;
  
  // correct for numerical issues
  for( int i=0; i<n; i++){
    if( xb[i] > 10 ){
      xb[i] = 10;
    }else if( xb[i] < -10 ){
      xb[i] = -10;
    }else{
      xb[i] = xb[i];
    }
  } 
  
  // apply inverse link function
  pi = invlogit(xb);
  
  // compute log likelihood 
  for(int i=0; i<n; i++){
    likelihood = likelihood + R::dbinom(Y[i], 1, pi[i], 1);
  }
  
  // compute log prior contribution for each parameter
  for (int i = 0; i < p; i++)
  {
	  prior = prior + R::dnorm4(b[i], 0, 10, 1);
  }

  // compute the proposal density
  for (int i = 0; i < p; i++)
  {
	  proposal = proposal + R::dnorm4(b_s[i], b[i], sqrt(jump_v), 1);
  }
  
  // evaluate log posterior as sum of likelihood and prior
  log_posterior = likelihood + prior;
  return log_posterior;
}

// [[Rcpp::export]]
List for_log_post(mat X, vec Y, int iteration, double jump_v){
    // create shell
    int p = X.n_cols;
    mat BETA(iteration, p);
    vec accept_rate(iteration);
    
	// starting value
    BETA.row(0) = ones<vec>(p).t();
    vec beta_0(p);
    mat sd(p, p);
    vec beta_c(p);
    double r_denum = 0;
    double r_num = 0;
    double r = 0;
    double rmin = 0;
    double uni = 0;

    for(int i = 1; i < iteration; i++){
        beta_0 = trans(BETA.row(i-1));

		sd = diagmat(jump_v * ones<vec>(p));
        
		// sample from proposal distribution
        beta_c = mvnrnd(beta_0, sd);
        
		// ratio computation
        r_num = log_post(beta_c, Y, X, beta_0, jump_v);
        
		r_denum = log_post(beta_0, Y, X, beta_c, jump_v);
        
		// acceptance ratio
        r = exp(r_num - r_denum);
        rmin = fmin(r, 1);
        
        //decision
        uni = randu<double>();
        if (rmin > uni){
            BETA.row(i) = beta_c.t();
        } else {
            BETA.row(i) = beta_0.t();
        }
        accept_rate(i) = rmin;
    }
    List result;
    result["BETA"] = BETA;
    result["accept_rate"] = accept_rate;
    return result;
}
