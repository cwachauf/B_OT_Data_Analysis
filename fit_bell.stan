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
  int<lower=1> N; 
  real forces[N]; 
  real<lower=0> lifetimes[N];
  int num_lt_steps[N];
  real<lower=0> rel_uncertainty;
}

parameters
{
  real dx; // persistence lenth[nm]
  real<lower=0> tau0;
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
  real forces_calibrated[N];
  real taus[N];
  
  dx ~ normal(0,5); // flat normal prior
  tau0 ~ cauchy(0,2.5);
  factor ~ normal(1.0,rel_uncertainty); // allow for a small correction in the force calibration
  for(i in 1:N)
  {
    forces_calibrated[i] <- forces[i]*factor;
    taus[i] <- tau0*exp(-(dx*forces_calibrated[i])/(kB*Temp));
    increment_log_prob(LogProbMeanExp(lifetimes[i],num_lt_steps[i],taus[i]));
  }
}