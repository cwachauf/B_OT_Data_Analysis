source("polym_mechanics.R")

Occupation_Probability_From_Force_Open <- function(F_open,L,P,k_eff,G_int,dz,Temp=298.15)
{
  kB <- 1.38e-02
  G_open_trap <- (0.5*F_open^2)/(k_eff)
  x_open_leash <- Extension_From_Force_WLC(F_open,L,P,Temp)
  G_open_leash <- Integrated_WLC(x_open_leash,L,P,Temp)
  G_open <- G_open_leash + G_open_trap
  
  x_closed_trap <- x_open_leash-2*dz + F_open/k_eff
  F_closed <- k_eff*x_closed_trap
  G_closed_trap <- (0.5*F_closed^2)/(k_eff)
  G_closed_leash <- Integrated_WLC(2*dz,L,P,Temp)
  G_closed <- G_closed_leash + G_closed_trap + G_int
  exp_open <- exp(-G_open/(kB*Temp))
  exp_closed <- exp(-G_closed/(kB*Temp))
  p_open <- exp_open/(exp_open+exp_closed)
  return(p_open)
}

Occupation_Probability_From_Force_Closed <- function(F_closed,L,P,k_eff,G_int,dz,Temp=298.15)
{
  kB <- 1.38e-02
  G_closed_trap <- (0.5*F_closed^2)/(k_eff)
  G_closed_leash <- Integrated_WLC(2*dz,L,P,Temp)
  G_closed <- G_closed_trap+G_closed_leash+G_int
  
  x_open_total <- F_closed/k_eff + 2*dz
  F_open <- Force_From_Extension_WLC_Hook(x_open_total,L,P,k_eff,Temp)
  x_open_leash <- Extension_From_Force_WLC(F_open,L,P,Temp)
  G_open_trap <- (0.5*F_open^2)/(k_eff)
  G_open_leash <- Integrated_WLC(x_open_leash,L,P,Temp)
  G_open <- G_open_trap + G_open_leash
  exp_open <- exp(-G_open/(kB*Temp))
  exp_closed <- exp(-G_closed/(kB*Temp))
  p_closed <- exp_closed/(exp_open+exp_closed)
  return(p_closed)
}

## Implementation of "Berkemeier-Schlierf"-model in R
Lifetime_From_Force_BS_Open <- function(F_open,tau0_open,L_trans,L_open,P,k_eff,Temp=298.15)
{
  kB <- 1.38e-02
  G_open_trap <- (0.5*F_open^2)/(k_eff)
  x_open_leash <- Extension_From_Force_WLC(F_open,L_open,P,Temp)
  G_open_leash <- Integrated_WLC(x_open_leash,L_open,P,Temp)
  G_open <- G_open_leash + G_open_trap
  
  x_trans_total <- x_open_leash + F_open/k_eff
  F_trans <- Force_From_Extension_WLC_Hook(x_trans_total,L_trans,P,k_eff,Temp)
  x_trans_leash <- Extension_From_Force_WLC(F_trans,L_trans,P,Temp)
  
  G_trans_leash <- Integrated_WLC(x_trans_leash,L_trans,P,Temp)
  G_trans_trap <- (0.5*F_trans^2)/(k_eff)
  G_trans <- G_trans_leash + G_trans_trap
    
  tau_open <- tau0_open*exp(-(G_open-G_trans)/(kB*Temp))
  return(tau_open)
}

Lifetime_From_Force_BS_Closed <- function(F_closed,tau0_closed,L_trans,P,k_eff,Temp=298.15)
{
  kB <- 1.38e-02
  G_closed_trap <- (0.5*F_closed^2)/(k_eff)
  G_closed <- G_closed_trap
  
  x_trans_total <- F_closed/k_eff
  F_trans <- Force_From_Extension_WLC_Hook(x_trans_total,L_trans,P,k_eff,Temp)
  G_trans_trap <- (0.5*F_trans^2)/(k_eff)
  x_trans_leash <- Extension_From_Force_WLC(F_trans,L_trans,P,Temp)
  G_trans_leash <- Integrated_WLC(x_trans_leash,L_trans,P,Temp)
  G_trans <- G_trans_leash + G_trans_trap
  tau_closed <- tau0_closed*exp(-(G_closed-G_trans)/(kB*Temp))
  return(tau_closed)
}

## Fit_Bell_Stan(forces,lifetimes,num_steps,dx0,tau0_0)
Fit_Bell_Stan <- function(forces,lifetimes,num_steps,dx0,tau0_0,factor0,num_iter,num_chains)
{
  init_stan_bell <- list(list(dx=dx0,tau0=tau0_0,factor=factor0))
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

TestStuff <- function(df_sim)
{
 ## df_sim <- Simulate_Bell_Model(-1.0,1e-03,12,2.0,8.0)
  ##stan_fit <- Fit_Bell_Stan(df_sim$forces,df_sim$lts,df_sim$ns,0,1e-02,1.0,20000,1)
  ##return(stan_fit)
  dx0 <- -1.0
  tau0_0 <- 1e-03
  factor0 <- 1.0
  num_steps <- df_sim$ns
  forces <- df_sim$forces
  lifetimes <- df_sim$lts
  num_iter <- 400000
  num_chains <- 1
  init_stan_bell <- list(list(dx=dx0,tau0=tau0_0,factor=factor0))
  
  data_stan_bell <- list(N=length(forces),forces=forces,lifetimes=lifetimes,num_lt_steps=num_steps,rel_uncertainty=0.001)
  fit_stan_bell1 <-  stan(file="fit_bell.stan",data=data_stan_bell,init=init_stan_bell,iter=num_iter,chains=num_chains)
  mat1 <- as.matrix(fit_stan_bell1)
  
  data_stan_bell <- list(N=length(forces),forces=forces,lifetimes=lifetimes,num_lt_steps=num_steps,rel_uncertainty=0.1)
  fit_stan_bell2 <-  stan(file="fit_bell.stan",data=data_stan_bell,init=init_stan_bell,iter=num_iter,chains=num_chains)
  mat2 <- as.matrix(fit_stan_bell2)
  
  data_stan_bell <- list(N=length(forces),forces=forces,lifetimes=lifetimes,num_lt_steps=num_steps,rel_uncertainty=0.2)
  fit_stan_bell3 <-  stan(file="fit_bell.stan",data=data_stan_bell,init=init_stan_bell,iter=num_iter,chains=num_chains)
  mat3 <- as.matrix(fit_stan_bell3)
  
  breaks <- seq(from=-3.8,to=0.2,by=0.02)
  
  ## create histograms and plot them...
##  hist1 <- hist(mat1[,1],breaks=breaks)
  ##hist2 <- hist(mat2[,1],breaks=breaks)
  ##hist3 <- hist(mat3[,1],breaks=breaks)
  ##plot(hist1,col=rgb(1,0,0,0.25))
  ##plot(hist2,col=rgb(0,0,1,0.25),add=T)
  ##plot(hist3,col=rgb(0,1,0,0.25),add=T)
  
  df_res <- data.frame(mat1,mat2,mat3)
  
}

Stan_Fit_Bell <- function(df_data,rel_unc=0.0001,num_iter=20000)
{
  require("rstan")
  num_chains <- 1
  data_stan_bell <- list(N_open=length(df_data$forces_open),N_closed=length(df_data$forces_closed),forces_open=df_data$forces_open,forces_closed=df_data$forces_closed,
                        lts_open=df_data$lts_open,lts_closed=df_data$lts_closed,num_ltsteps_open=df_data$num_ltsteps_open,num_ltsteps_closed=df_data$num_ltsteps_closed,rel_uncertainty=rel_unc)
  init_stan_bell <- list(list(dx_off=-0.5,dx_on=6.0,tau0_open=1e-3,tau0_closed=1e02,factor=1.00))
  fit_stan_bell <- stan(file="fit_bell_oc.stan",data=data_stan_bell,init=init_stan_bell,iter=num_iter,chains=num_chains)
  return(fit_stan_bell)
}

Make_NLS_Fits <- function(df_data)
{
  require("minpack.lm")
  forces_open <- df_data$forces_open
  lts_open <- df_data$lts_open
  wts_open <- df_data$sem_open
  data <- list(x=forces_open,y=lts_open)
  kBT <- 1.38e-02*298.15;
  wts_open <- as.vector(wts_open)
  ##print(wts)
  
  nls_fit_open <- nls(formula <- lts_open ~ tau0*exp(dx*forces_open/kBT),start=c(tau0=1e-03,dx=5),weights=1/wts_open^2)#c(0.1,0.2,0.2,0.1,0.2,0.2,0.1,0.2,0.2,0.1,0.2,0.2))#1/wts^2)
  coefs_open <- coefficients(nls_fit_open)
  
  forces_closed <- df_data$forces_closed
  lts_closed <- df_data$lts_closed
  wts_closed <- df_data$sem_closed
  wts_closed <- as.vector(wts_closed)
  nls_fit_closed <- nls(formula <- lts_closed ~ tau0*exp(dx*forces_closed/kBT),start=c(tau0=3e01,dx=-1.5),weights=1/wts_closed^2)
  coefs_closed <- coefficients(nls_fit_closed)
  
  df_nls_results <- data.frame(coefs_open,coefs_closed)
  summary(nls_fit_open)
  summary(nls_fit_closed)
  print(summary(nls_fit_open))
  print(summary(nls_fit_closed))
  return(df_nls_results)
}

#TestAll <- function()
#{
#  
#  ## on -Lebenszeit nimmt ab:
#  
#  print(head(data_closed))
#  plot(data_closed$forces,data_closed$lts,log="y",xlim=c(0,10),ylim=c(1e-4,1e02),xlab="force [pN]",ylab="lifetime [s]")
#  points(data_open$forces,data_open$lts,col="red")
#  
#  ## test: do inference on the data:
#  data_stan_bell <- list(N_open=length(data_open$forces),N_closed=length(data_closed$forces),forces_open=data_open$forces,forces_closed=data_closed$forces,
#                         lts_open=data_open$lts,lts_closed=data_closed$lts,num_ltsteps_open=data_open$ns,num_ltsteps_closed=data_closed$ns,rel_uncertainty=0.0001)
#  init_stan_bell <- list(list(dx_off=-0.5,dx_on=6.0,tau0_open=1e-3,tau0_closed=1e02,factor=1.00))
#                         
#  fit_stan_bell <- stan(file="fit_bell_oc.stan",data=data_stan_bell,init=init_stan_bell,iter=num_iter,chains=num_chahttp://127.0.0.1:22908/graphics/plot_zoom_png?width=563&height=643ins)
#  return(fit_stan_bell)
#}

TestBSPlot <- function()
{
  L_trans <- 20
  L_open <- 40
  P <- 0.9
  k_eff <- 0.15
  tau0_open <- 1e-04
  tau0_closed <- 1e02
  forces <- seq(from=1.5,to=10,by=0.1)

  G_int <- -18
  popen <- array(0,dim=c(length(forces)))
  pclosed <- array(0,dim=c(length(forces)))
  for(i in 1:length(forces))
  {
    popen[i] <- Occupation_Probability_From_Force_Open(forces[i],L_open,P,k_eff,G_int,dz=3.5,Temp=298.15)
    pclosed[i] <- Occupation_Probability_From_Force_Closed(forces[i],L_open,P,k_eff,G_int,dz=3.5,Temp=298.15)
    
  }
  plot(forces,popen,type="l")
  points(forces,pclosed,type="l")
}
## PlotBell(dx,tau0,fmin,fmax,npnts,Temp=298.15)
PlotBell <- function(dx,tau0,fmin,fmax,npnts,Temp=298.15)
{
  kB <- 1.38e-02;
  forces <- seq(from=fmin,to=fmax,length=npnts)
  lts <- array(0,dim=c(npnts))
  for(i in 1:length(forces))
  {
    lts[i] <- tau0*exp(forces[i]*dx/(kB*Temp))
  }
  
  df_plot <- data.frame(forces,lts);
  return(df_plot);
}

Plot_Normal_Distribution <- function(mu,sigma,x_min,x_max,num_points)
{
  x_values <- seq(from=x_min,to=x_max,length=num_points)
  y_values <- matrix(nrow=num_points,ncol=1)
  for(i in 1:num_points)
  {
    y_values[i] <- dnorm(x_values[i],mean=mu,sd=sigma)
  }
  df_norm_plot <- data.frame(x_norm_dist=x_values,y_norm_dist=y_values)
  return(df_norm_plot)
}

## HMM_Calculate_Lifetime_LogLikelihood(t_m,tau,n)
## calculates the log likelihood for the mean value (t_m)
## out of n exponentially distributed RVs with underlying 
## lifetime tau, see description above
HMM_Calculate_Lifetime_LogLikelihood <- function(t_m,tau,n)
{
  log_lik <- log(n) - n*log(tau)+(n-1)*log(n*t_m)-n*t_m/tau - lgamma(n)
  return(log_lik)
}

