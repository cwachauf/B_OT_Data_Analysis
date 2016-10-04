source("polym_mechanics.R")

Load_Extension_Gain_Data <- function(index,num_iters)
{
  full_path <- "/Users/christianwachauf/Documents/Daten/data_stack_faki/Ext_Gains_vs_Force/"
  filename <- paste0("ext_gain_A0",toString(index),".txt")
  full_path <- paste0(full_path,filename)
  df_ext_gain_data <- read.table(full_path,sep="\t",header=TRUE)
  head(df_ext_gain_data)
  
  plot(df_ext_gain_data[,3],df_ext_gain_data[,4],xlim=c(0,8),ylim=c(0,11),xlab="force [pN]",ylab="extension gain [nm]")
##  df_ext_gain_theor <- Plot_Extension_Gain_From_Force_Open(L=40,P=0.88,Temp=298.15,dz=3.5,F_min=2,F_max=6,npnts=200)
 ## points(df_ext_gain_theor$forces_open,df_ext_gain_theor$ext_gains_open,type="l")
  sf_ext_gain_data <- Stan_Fit_Extension_Gain_Data(forces_open=df_ext_gain_data[,3],extension_gains = df_ext_gain_data[,4],num_iters=num_iters)
  df_prediction_bands <- Plot_Stan_Extension_Gain_From_Force_Open(sf_ext_gain_data,F_min=2,F_max=7,npnts=100,nsamples=num_iters/2)
  points(df_prediction_bands$forces,df_prediction_bands$ext_gain_mean,type="l")
  points(df_prediction_bands$forces,df_prediction_bands$ext_gain_lp,type="l",lty=2)
  points(df_prediction_bands$forces,df_prediction_bands$ext_gain_up,type="l",lty=2)
  return(sf_ext_gain_data)
}

Stan_Fit_Extension_Gain_Data <- function(forces_open,extension_gains,num_iters=10000)
{
  require("rstan")
  stan_data <- list(N=length(forces_open),forces_open=forces_open,ext_gains_fopen=extension_gains)
  stan_init <- list(list(P=1.0,CL=40.0,sigma=0.5))
  stan_fit_ext_gain_data <- stan(file="/Users/christianwachauf/Documents/Scripts/GithubRepo/B_OT_Data_Analysis/Stan/fit_ext_gain.stan",data=stan_data,init=stan_init,iter=num_iters,chain=1)
  return(stan_fit_ext_gain_data)
}

Plot_Stan_Extension_Gain_From_Force_Open <- function(stan_fit_ext_gain_data,F_min=2,F_max=7,npnts,nsamples)
{
  require("coda")
  mat <- as.matrix(stan_fit_ext_gain_data)
  
  forces <- seq(from=F_min,to=F_max,length=npnts)
  
  ext_gain_mean <- array(0,dim=c(npnts))
  ext_gain_lp <- array(0,dim=c(npnts))
  ext_gain_up <- array(0,dim=c(npnts))
  
  predictions_ext_gain <- array(0,dim=c(nsamples,npnts))
  
  for(i in 1:nsamples)
  {
    curr_P <- mat[i,1]
    curr_CL <- mat[i,2]
    curr_sigma <- mat[i,3]
    curr_Temp <- mat[i,4]
    curr_dz <- mat[i,5]
    curr_k_eff <- mat[i,6]
    
    
    for(j in 1:npnts)
    {
      predictions_ext_gain[i,j] <- Extension_Gain_From_Force_Open(forces[j],curr_CL,curr_P,curr_Temp,curr_dz)
      predictions_ext_gain[i,j] <- predictions_ext_gain[i,j] + rnorm(n=1,mean=0,sd=curr_sigma)
    }
  }
  for(j in 1:npnts)
  {
    ext_gain_mean[j] <- mean(predictions_ext_gain[,j])
    hdi_values <- HPDinterval(as.mcmc(predictions_ext_gain[,j]))
    ext_gain_lp[j] <- hdi_values[1];
    ext_gain_up[j] <- hdi_values[2];
  }
  df_preds <- data.frame(forces,ext_gain_mean,ext_gain_lp,ext_gain_up)
  return(df_preds)
}

Extension_Gain_From_Force_Open <- function(F_open,L,P,Temp=298.15,dz=3.5)
{
  x_leash_open <- Extension_From_Force_WLC(F_open,L,P,Temp)
  dx <- x_leash_open - 2*dz
  return(dx)
}

Plot_Extension_Gain_From_Force_Open <- function(L,P,Temp=298.15,dz=3.5,F_min,F_max,npnts)
{
  forces_open <- seq(from=F_min,to=F_max,length=npnts)
  ext_gains_open <- array(0,dim=c(npnts))
  for(i in 1:npnts)
  {
    ext_gains_open[i] <- Extension_Gain_From_Force_Open(forces_open[i],L,P,Temp)
  }
  df_ext_gain <- data.frame(forces_open,ext_gains_open)
  return(df_ext_gain)
}
