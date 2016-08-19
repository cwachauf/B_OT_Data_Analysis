source("B_OT_data_analysis.R")
## Create some synthetic data to test the bell fit functions
Simulate_Bell_Model <- function(dx,tau0,nforces,f_min,f_max)
{
  kBT <- 1.38e-02*298.15 ## assume Temp = 25°C
  forces <- runif(n=nforces,min=f_min,max=f_max)
  num_lts <- sample(5:50,length(forces),replace=T)
  taus <- array(0,dim=c(length(forces)))
  lts <- array(0,dim=c(length(forces)))
  sem <- array(0,dim=c(length(forces)))
  for(i in 1:length(forces))
  {
    taus[i] <- tau0*exp(forces[i]*dx/(kBT))
    rexp_vs <- rexp(n=num_lts[i],rate=(1.0/taus[i]))
    lts[i] <- mean(rexp_vs)
    sem[i] <- sd(rexp_vs)/sqrt(num_lts[i])
  }
  df_bell_data <- data.frame(forces=forces,lts=lts,ns=num_lts,sem=sem)
  return(df_bell_data)
}

## Simulate_Data_Set() creates synthetic data for a whole system
## (i.e. the open and the closed branch in a Chevron plot)...
Simulate_Data_Set <- function()
{
  dx_on <- 8.0
  dx_off <- -1.0
  tau_closed <- 1e01
  tau_open <- 1e-04
  n_forces_open <- 12
  n_forces_closed <- 12
  f_min <- 2.0
  f_max <- 8.0
  
  num_iter <- 20000;
  num_chains <- 1;
  data_closed <- Simulate_Bell_Model(dx_off,tau_closed,n_forces_closed,f_min,f_max)
  data_open <- Simulate_Bell_Model(dx_on,tau_open,n_forces_open,f_min,f_max)
  df_data <- data.frame(forces_open=data_open$forces,lts_open=data_open$lts,num_ltsteps_open=data_open$ns,sem_open=data_open$sem,
                        forces_closed=data_closed$forces,lts_closed=data_closed$lts,num_ltsteps_closed=data_closed$ns,sem_closed=data_closed$sem)
  return(df_data)
}

## PlotDataSet(df_data)
## Plots the entire data set (i.e. open and closed branch)
## on a logarithmic scale
## adds error bars (the standard errors of the mean)
Plot_Data_Set <- function(df_data)
{
  epsilon <- 0.02
  plot(df_data$forces_open,df_data$lts_open,xlim=c(0,10),log="y",ylim=c(1e-04,1e03),xlab="force [pN]",ylab="lifetime [s]")
  
  points(df_data$forces_closed,df_data$lts_closed,col="red")
  segments(df_data$forces_open,df_data$lts_open-df_data$sem_open,df_data$forces_open,df_data$lts_open+df_data$sem_open)
  segments(df_data$forces_open-epsilon,df_data$lts_open-df_data$sem_open,df_data$forces_open+epsilon,df_data$lts_open-df_data$sem_open)
  segments(df_data$forces_open-epsilon,df_data$lts_open+df_data$sem_open,df_data$forces_open+epsilon,df_data$lts_open+df_data$sem_open)
  
  segments(df_data$forces_closed,df_data$lts_closed-df_data$sem_closed,df_data$forces_closed,df_data$lts_closed+df_data$sem_closed,col="red")
  segments(df_data$forces_closed-epsilon,df_data$lts_closed-df_data$sem_closed,df_data$forces_closed+epsilon,df_data$lts_closed-df_data$sem_closed,col="red")
  segments(df_data$forces_closed-epsilon,df_data$lts_closed+df_data$sem_closed,df_data$forces_closed+epsilon,df_data$lts_closed+df_data$sem_closed,col="red")
}

## Plot_All(df_data)
## plots the data set as 
## well as the nonlinear least squares fits
Plot_All <- function(df_data)
{
  Plot_Data_Set(df_data)
  nls_fit_coefs <- Make_NLS_Fits(df_data)
  df_plot_open <- PlotBell(nls_fit_coefs[2,1],nls_fit_coefs[1,1],0.0,10.0,200)
  df_plot_closed <- PlotBell(nls_fit_coefs[2,2],nls_fit_coefs[1,2],0.0,10.0,200)
  
  points(df_plot_open$forces,df_plot_open$lts,type="l")
  points(df_plot_closed$forces,df_plot_closed$lts,type="l")
  
  ## print out the coefficient values,....
  print(nls_fit_coefs)
  kB <- 1.38e-02
  Temp <- 298.15
  dG <- -(kB*Temp)*log(nls_fit_coefs[1,2]/nls_fit_coefs[1,1]);
  print("dG_ML: ")
  print(dG)
  
  dx_1 <- nls_fit_coefs[2,2]
  
  
  sf_bell1 <- Stan_Fit_Bell(df_data,rel_unc=0.0001,num_iter=400000)
  sf_bell2 <- Stan_Fit_Bell(df_data,rel_unc=0.15,num_iter=400000)
  ##sf_bell3 <- Stan_Fit_Bell(df_data,rel_unc=0.15,num_iter=40000)
  
  breaks = seq(from=-5,to=+1,by=0.02)
  breaks2 = seq(from=-65,to=-20,by=0.2)
  
  mat1 <- as.matrix(sf_bell1)
  mat2 <- as.matrix(sf_bell2)
  print(summary(sf_bell1))
  print(summary(sf_bell2))
##  mat3 <- as.matrix(sf_bell3)
  
  hist1 <- hist(mat1[,1],breaks=breaks)
  hist2 <- hist(mat2[,1],breaks=breaks)
  
  hist1_1 <- hist(mat1[,8],breaks=breaks2)
  hist2_1 <- hist(mat2[,8],breaks=breaks2)
  
 ## hist3 <- hist(mat3[,8],breaks=breaks)
  
  mw1 <- mean(mat1[,1])
  mw2 <- mean(mat2[,1])
  sd1 <- sd(mat1[,1])
  sd2 <- sd(mat2[,1])
  print(mw1)
  print(mw2)
  print(sd1)
  print(sd2)
  mw1 <- mean(mat1[,8])
  mw2 <- mean(mat2[,8])
  sd1 <- sd(mat1[,8])
  sd2 <- sd(mat2[,8])
  print(mw1)
  print(mw2)
  print(sd1)
  print(sd2)
  
  plot(hist1,col=rgb(1,0,0,0.25),xlab="dx [nm]",ylab="frequency",main="uncertainty in k_eff",xlim=c(-3,0))
  plot(hist2,col=rgb(0,0,1,0.25),add=T)
  xs_dx1 <- c(dx_1,dx_1)
  vs <- par("usr")
  ys_dx1 <- c(0,vs[4])
  points(xs_dx1,ys_dx1,type="l",lwd=2,lty=2)
  plot(hist1_1,col=rgb(1,0,0,0.25),xlab="dG [pN/nm]",ylab="frequency",main="uncertainty in k_eff",xlim=c(-55,-40))
  plot(hist2_1,col=rgb(0,0,1,0.25),add=T)
  
  
  
 ## plot(hist3,col=rgb(0,1,0,0.25),add=T)
  
  ## add line 
  
  vs <- par("usr")
  
  xs_dg_ml <- c(dG,dG)
  ys_dg_ml <- c(0,vs[4])
  points(xs_dg_ml,ys_dg_ml,type="l",lwd=2,lty=2)
  ## plot te 
}