########################
#############Beta regression (proportional outcome)
#############treatment and control models fit separately for efficiency

library(rstan)

# load file names for stan
# for DAW method
#file<-"Daw Beta Constant Phi.stan"
# for power prior method
filename<-"Beta Reg Constant Phi.stan"


# Treated ######################################################################
# covariate matrices for treated current and historical data
X1A<-model.matrix(object= ~hkage1+factor(race)+hksex+prior_nh+prior_inp+adrd+
                    ami+anem+atrfib+chf+cancer+chrkid+copd+
                    depr+diab+hyperten+hyprlip+ischhd+osteo+stroke+smi, 
                  data=dat[dat$study==1 & dat$grp==1,])

# historical
X1hA<-model.matrix(object= ~hkage1+factor(race)+hksex+prior_nh+prior_inp+adrd+
                     ami+anem+atrfib+chf+cancer+chrkid+copd+
                     depr+diab+hyperten+hyprlip+ischhd+osteo+stroke+smi,
                   data=dat[dat$study==0& dat$grp==1,])


# create data list for stan 
dataListA=list(# current data
              # number of observations
              N=nrow(dat[dat$study==1& dat$grp==1,]),
              # proportion outcome vector
              y=dat[dat$study==1& dat$grp==1,]$prop,
              # number of columns in model matrix
              X1_dim=dim(X1A)[2],
              # covariate matrix
              X1=X1A,
              
              # historical data
              # number of observations
              Nh=nrow(dat[dat$study==0& dat$grp==1,]),
              # proportion outcome vector
              yh=dat[dat$study==0& dat$grp==1,]$y,
              # number of columns in model matrix
              X1_dimh=dim(X1hA)[2],
              # covariate matrix
              X1h=X1hA,
              
              # "on-trial" score weights for DAW
              # weight=dat[dat$study==0& dat$grp==1,]$iow.weight,
              
              # power prior alpha
              alpha=0.8
)

# fit treated model
trt_fit<-rstan::stan(file=filename, data = dataListA, chains = 1, iter = 5000,
                     thin = 10)
# extract beta coefficient samples
trt_sample_vals=t(extract(trt_fit)$beta_x1)
################################################################################



# Control ######################################################################
# covariate matrices for historical data
X1B<-model.matrix(object= ~hkage1+factor(race)+hksex+prior_nh+prior_inp+adrd+
                    ami+anem+atrfib+chf+cancer+chrkid+copd+
                    depr+diab+hyperten+hyprlip+ischhd+osteo+stroke+smi, 
                  data=dat[dat$study==1 & dat$grp==0,])

# historical
X1hB<-model.matrix(object= ~hkage1+factor(race)+hksex+prior_nh+prior_inp+adrd+
                     ami+anem+atrfib+chf+cancer+chrkid+copd+
                     depr+diab+hyperten+hyprlip+ischhd+osteo+stroke+smi, 
                   data=dat[dat$study==0& dat$grp==0,])

# data list for stan
dataListB=list(# current data
              # number of observations
              N=nrow(dat[dat$study==1& dat$grp==0,]),
              # proportion outcome vector
              y=dat[dat$study==1& dat$grp==0,]$y,
              # number of columns in model matrix
              X1_dim=dim(X1B)[2],
              # covariate matrix
              X1=X1B,
               
              # historical data
              # number of observations
              Nh=nrow(dat[dat$study==0& dat$grp==0,]),
              # proportion outcome vector
              yh=dat[dat$study==0& dat$grp==0,]$y,
              # number of columns in model matrix
              X1_dimh=dim(X1hB)[2],
              # covariate matrix
              X1h=X1hB,
              
              # "on-trial" score weights for DAW
              #weight=dat[dat$study==0& dat$grp==0,]$iow.weight,
              
              # power prior alpha
              alpha=0.8
)

# fit control model
ctrl_fit<-rstan::stan(file=filename, data = dataListB, chains = 1, iter = 5000,
                      thin = 10)
# extract beta coefficient samples
ctrl_sample_vals=t(extract(ctrl_fit)$beta_x1)

# for each posterior sample
est<-c()

# impute values for Y1 for control group and Y0 for treated group
for(i in 1:ncol(trt_sample_vals)){
  # estimate trt value for control group
  trt_impute<-rbind(X1hB, X1B)%*%as.matrix(trt_sample_vals[,i])
  # estimate ctrl value for trt group
  ctrl_impute<-rbind(X1hA, X1A)%*%as.matrix(ctrl_sample_vals[,i])
  
  # Y(1) for treated and controls (imputed)
  Y1<-c(dat[dat$grp==1,]$prop, expit(trt_impute))
  # Y(0) for treated (imputed) and controls
  Y0<-c(expit(ctrl_impute),dat[dat$grp==0,]$prop)
  # calculate ate
  ate<-mean(Y1)-mean(Y0)
  est<-c(est, ate)
}

# summarize ATE estimate by taking mean and variance of samples
ate_est<-c(mean(est), var(est))
