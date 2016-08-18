source("polym_mechanics.R")

## Create some synthetic data to test the bell fit functions
Simulate_Bell_Model <- function(dx,tau0,nforces,f_min,f_max)
{
  ## for the simulation:
  ## assume Temp = 25°C
  kBT <- 1.38e-02*298.15
  forces <- runif(n=nforces,min=f_min,max=f_max)
  num_lts <- sample(20:50,length(forces),replace=T)
  taus <- array(0,dim=c(length(forces)))
  lts <- array(0,dim=c(length(forces)))
  for(i in 1:length(forces))
  {
    taus[i] <- tau0*exp(-forces[i]*dx/(kBT))
    lts[i] <- rgamma(n=1,shape=num_lts[i],scale=taus[i])/num_lts[i]
  }
  df_bell_data <- data.frame(forces=forces,lts=lts,ns=num_lts)
  return(df_bell_data)
}

## Fit_Bell_Stan(forces,lifetimes,num_steps,dx0,tau0_0)
Fit_Bell_Stan <- function(forces,lifetimes,num_steps,dx0,tau0_0,num_iter,num_chains)
{
  init_stan_bell <- list(list(dx=dx0,tau0=tau0_0))
  data_stan_bell <- list(N=length(forces),forces=forces,lifetimes=lifetimes,num_lt_steps=num_steps)
  fit_stan_bell <-  stan(file="fit_bell.stan",data=data_stan_bell,init=init_stan_bell,iter=num_iter,chains=num_chains)
  return(fit_stan_bell)
}


## Extension_Gain_From_Force_Open_WLC(f_open,L,P,Temp=298.15,dz)
## Takes as an input the force in the open conformation "f_open" (in pN),
## the contour length "L" (in nm), the persistence length "P" (in nm),
## the temperature "Temp" (in K), default 298.15 and the distance from
## the structure end to the beginning of the leash on one side "dz" (in nm)
## and calculates the extension gain
Extension_Gain_From_Force_Open_WLC <- function(f_open,L,P,Temp=298.15,dz)
{
  x_leash_open <- Extension_From_Force_WLC(f_open,L,P,Temp=298.15)
  dx <- x_leash_open - 2*dz;
  return(dx)
}

## Extension_Gain_From_Force_Closed_WLC(f_closed,L,P,Temp=298.15,dz,k_eff)
## calculates the extension gain from a given force in the closed conformation
## for a system of characterized by the contour length "L" (in nm), the persistence length "P" (in nm),
## the temperature "Temp" (in K), default 298.15 and the distance from
## the structure end to the beginning of the leash on one side "dz" (in nm)
## as well as the effective spring constant "k_eff" (in pN/nm)
Extension_Gain_From_Force_Closed_WLC <- function(f_closed,L,P,Temp=298.15,dz,k_eff)
{
  F_max <- 40
  npnts_array <- 400;
  x_open <- f_closed/k_eff + 2*dz
  
  ## interpolate the force in the open state
  forces <- seq(from=0,to=F_max,length=npnts_array)
  extensions <- Create_Extension_Array_From_Forces_WLC_Hook <- function(forces,L,P,k_eff,Temp=298.15)
  df_res <- approx(x=extensions,y=forces,xout=x_open)
  f_open <- df_res$y
  
  x_leash_open <- Extension_From_Force_WLC(f_open,L,P,Temp=298.15)
  dx <- x_leash_open - 2*dz;
  return(dx)
}

TestStuff <- function()
{
  df_sim <- Simulate_Bell_Model(-1.0,1e-03,12,2.0,8.0)
  stan_fit <- Fit_Bell_Stan(df_sim$forces,df_sim$lts,df_sim$ns,0,1e-02,20000,1)
  return(stan_fit)
}



