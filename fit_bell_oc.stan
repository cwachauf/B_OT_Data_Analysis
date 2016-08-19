functions
{
  real LogProbMeanExp(real mean,int num,real lt)
  {
    int n;
    real t;
    real tau;
    real log_lik;
    n <- num;
    t <- mean;
    tau <- lt;
    log_lik <- log(n) + n*log(1.0/tau) +(n-1)*log(n*t) - n*t/tau - lgamma(n);
    return(log_lik);
  }
}

data 
{
  int<lower=1> N_open;
  int<lower=1> N_closed;
  
  real forces_open[N_open];
  real forces_closed[N_closed];
  
  real lts_open[N_open];
  real lts_closed[N_closed];
  
  int num_ltsteps_open[N_open];
  int num_ltsteps_closed[N_closed];
  
  real<lower=0> rel_uncertainty;
}

parameters
{
  real dx_off; // persistence lenth[nm]
  real dx_on;
  
  real<lower=0> tau0_open;
  real<lower=0> tau0_closed;
  
  real<lower=0> factor;
}

transformed parameters
{
  real<lower=0> Temp;
  real<lower=0> kB;
  kB <- 1.38e-02;
  Temp <- 298.15;
}

model 
{
  real fopen_calib[N_open];
  real fclosed_calib[N_closed];
  
  real taus_open[N_open];
  real taus_closed[N_closed];
  
  dx_on ~ normal(0,20);
  dx_off ~ normal(0,20);
  
  tau0_open ~ cauchy(0,2.5);
  tau0_closed ~ cauchy(0,2.5);
  
  factor ~ normal(1.0,rel_uncertainty); // allow for a small correction in the force calibration
  
  for(i in 1:N_open)
  {
    fopen_calib[i] <- forces_open[i]*factor;
    taus_open[i] <- tau0_open*exp(dx_on*fopen_calib[i]/(kB*Temp));  //(dx_on,tau_open,n_forces_open,f_min,f_max)
    increment_log_prob(LogProbMeanExp(lts_open[i],num_ltsteps_open[i],taus_open[i]));
   
  }
  for(i in 1:N_closed)
  {
    fclosed_calib[i] <- forces_closed[i]*factor;
    taus_closed[i] <- tau0_closed*exp(dx_off*fclosed_calib[i]/(kB*Temp));
    increment_log_prob(LogProbMeanExp(lts_closed[i],num_ltsteps_closed[i],taus_closed[i]));
  }
}

generated quantities
{
  real dG;
  dG <- -(kB*Temp)*log(tau0_closed/tau0_open);
}