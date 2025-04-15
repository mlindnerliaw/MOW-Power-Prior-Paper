 data { 
     //Outcome data current
     int<lower=0> N; 
     real<lower=0,upper=1> y[N]; 
     /* X1 matrix for the mean parameter. Include a constant column. */ 
     int<lower=1> X1_dim; 
     matrix[N, X1_dim] X1; 
     
     //Outcome data historical
     int<lower=0> Nh; 
     real<lower=0,upper=1> yh[Nh]; 
     /* X1 matrix for the mean parameter. Include a constant column. */ 
     int<lower=1> X1_dimh; 
     matrix[Nh, X1_dimh] X1h; 
     
     // alpha power
     real<lower=0, upper=1>alpha;
     
 } 
  
 transformed data { 
 } 
  
 parameters { 
     vector[X1_dim] beta_x1; 
     real<lower=0> phi;
 } 
  
 transformed parameters { 
     vector[N] eta_x1 = X1 * beta_x1; 
     vector[Nh] eta_x1h = X1h * beta_x1; 
     
  
     /* logit for mean. expit is the inverse. */ 
     vector[N] mu = inv_logit(eta_x1); 
     vector[Nh] muh = inv_logit(eta_x1h); 
} 
  
 model { 
   
  
     /* Priors */ 
     beta_x1~normal(0,1);
     
  
     /* Mean model */ 
     for (i in 1:N) { 
         target += beta_lpdf(y[i] | (mu * phi)[i], ((1-mu) * phi)[i]); 
     } 
     
     for (ii in 1:Nh) { 
         target += beta_lpdf(yh[ii] | (muh * phi)[ii], ((1-muh) * phi)[ii])*alpha; 
     }      
  
 } 
  
 generated quantities { 
  
 } 
 
