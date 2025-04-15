/// stan code for weibull model
functions {
  real weibull2_lpdf(real y, real alpha, real lambda) { // log density
    real prob;
    real lprob;
    prob = alpha * lambda * y^(alpha - 1) * exp(- lambda*y^alpha);
    lprob = log(prob);
    return lprob;
  }
  real weibull2_lcdf(real y, real alpha, real lambda) { // log cdf
    real prob;
    real lprob;
    prob = 1-exp(-lambda * y^alpha);
    lprob = log(prob);
    return lprob;
  }
  real weibull2_lccdf(real y, real alpha, real lambda) { // log survival
    real prob;
    real lprob;
    prob = exp(-lambda * y^alpha);
    lprob = log(prob);
    return lprob;
  }
  real weibull2_rng(real alpha, real lambda) { // generate random data points from survival
    real u;
    real z;
    u = uniform_rng(0, 1);
    z = (-log(u)/lambda)^(1/alpha);
    return z;
  }
  real weibull2_lb_rng(real alpha, real lambda, real lb) { // generate data points from weibull truncated below at lb 
    real p = exp(weibull2_lcdf(lb | alpha, lambda));
    real u;
    real z;
    if(is_nan(p)) p = 0.999;
    else if (p >= 1) p = 0.999;
    else if (p <= 0) p = 0;
    u = uniform_rng(p, 1);      
    z = (-log(1-u)/lambda)^(1/alpha);
    return z;
  }
}

 
data {
  int<lower=0> N;   //num obs current
  int<lower=0> Nh;  // num obs historical
  int<lower=0> K;   // num covariates
  real y[N];       // outcome
  real yh[Nh];      // outcome historical
  int v[N];       // event indicator
  int vh[Nh];     // event indicator historical
  matrix[N, K] X; // design matrix
  matrix[Nh,K] Xh; // design matrix historical
  real<lower=0> alpha_a0;  //hyperparameters
  real<lower=0> alpha_b0;
  real<lower=0> beta_sig;
  real<lower=0,upper=1> power;
}

parameters {
  real<lower=0> alpha;
  vector[K] beta;
}

transformed parameters {
  vector[N] lp;
  vector[Nh] lph;
  real shape;
  real scale;
  lp = exp(X*beta);
  lph = exp(Xh*beta);
  shape = alpha;
  scale = exp(beta[1]);
}

model {
  alpha ~ gamma(alpha_a0, alpha_b0);
  beta ~ normal(0, beta_sig);
  for(n in 1:N) { // current
    if(v[n] == 0) // right-censored
      target += weibull2_lccdf(y[n] | alpha, lp[n]);
    else if(v[n] == 1) // observed event time
      target += weibull2_lpdf(y[n] | alpha, lp[n]);
  }
  for (nh in 1:Nh) { //historical
    if(vh[nh] == 0) // right-censored
      target += weibull2_lccdf(yh[nh] | alpha, lph[nh])*power;
    else if(vh[nh] == 1) // observed event time
      target += weibull2_lpdf(yh[nh] | alpha, lph[nh])*power;
  }
}

      
              

generated quantities { // take draws from the posterior predictive distribution
} 
