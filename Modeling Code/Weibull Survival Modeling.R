########################
#############Weibull regression (survival outcome)

library(rstan)
library(flexsurv)
library(survival)


# load file name for stan
fileName<-"Weibull.stan"

# covariate matrix for current data
X.mat<-model.matrix(object= ~grp+hkage1+factor(race)+hksex+prior_nh+prior_inp+
                      adrd+ami+anem+atrfib+chf+cancer+chrkid+copd+
                      depr+diab+hyperten+hyprlip+ischhd+osteo+stroke+smi, 
                    data=dat[dat$study==1,])

# covariate matrix for historical data
Xh.mat<-model.matrix(object= ~grp+hkage1+factor(race)+hksex+prior_nh+prior_inp+
                       adrd+ami+anem+atrfib+chf+cancer+chrkid+copd+
                       depr+diab+hyperten+hyprlip+ischhd+osteo+stroke+smi, 
                     data=dat[dat$study==0,])


# data list for stan
dataList<-list(# current data
               # number of current observations
               N=nrow(dat[dat$study==1,]),
               # number of columns in model matrix
               K=ncol(X.mat),
               # vector of survival times (inpatient or nursing home)
               y=dat[dat$study==1,]$time,
               # vector of event indicators (inpatient or nursing home)
               v=dat[dat$study==1,]$event,
               # covariate matrix
               X=X.mat,
               
               # historical data
               # number of historical observations
               Nh=nrow(dat[dat$study==0,]),
               # vector of survival times (inpatient or nursing home)
               yh=dat[dat$study==0,]$time,
               # vector of event indicators (inpatient or nursing home)
               vh=dat[dat$study==0,]$event,
               # covariate matrix
               Xh=Xh.mat,
               
               # "on-trial" score weights for DAW
               #ph=dat[dat$study==0,]$iow.weight,
               
               #hyperparameters
               # gamma(a,b) prior for Weibull shape parameter a
               alpha_a0=1,
               alpha_b0=2,
               # sigma for normal prior on beta coefficients
               beta_sig=1,
               
               # power prior alpha
               power=0.4
)

# fit Weibull model in stan
stan_fit = rstan::stan(fileName, data = dataList, chains = 1, iter = 5000,
                       thin = 10)

# extract samples from model
post<-rstan::extract(stan_fit)

# initialize matrices to store values
cox<-matrix(0, nrow=length(post$alpha), ncol=2)
colnames(cox)<-c("grp_coef", "SE")

# for each posterior sample
for (i in 1:length(post$alpha)){
  # impute the data from the modeled Weibull distribution in PH form
  dat$impute<-rweibullPH(nrow(dat), 
                              shape = post$alpha[i], 
                              scale = exp(as.matrix(rbind(Xh.mat, X.mat))%*%
                                            post$beta[i,]))
  # new event indicator with 180 day follow up
  event1<-ifelse(dat$impute>180, 0,1)
  # censor imputed times based on follow-up period
  dat$impute<-ifelse(dat$impute>180, 180, dat$impute)
  # fit cox PH model on imputed data
  fit<-coxph(Surv(impute, event1)~grp+hkage1+factor(race)+hksex+prior_nh+
               prior_inp+adrd+ami+anem+atrfib+chf+cancer+chrkid+copd+
               depr+diab+hyperten+hyprlip+ischhd+osteo+stroke+smi, data=dat)
  s<-summary(fit)$coefficients
  
  # store grp coefficient estimate and standard error
  cox[i,]<-c(s[1,1], s[1,3])
}

# summarize treatment estimate from posterior samples using mean estimate, 
# mean variance of estimate, and variance of estimates
# the two variance measured required for two-stage multiple imputation
cox_est<-c(mean(cox[,1]), mean(cox[,2]^2), var(cox[,1]))  

