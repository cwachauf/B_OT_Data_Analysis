## Force_From_Extension_WLC(x,L,P,Temp=298.15)
## calculates the force at a given extension for a
## wormlike-chain model (interpolation formula from
## J.F. Marko and E.D. Siggia "Stretching DNA", Macromolecules 1995
## F(x) = (kBT/P)*(1/(4*(1-x/L)^2) - 1/4 +x/L)
## where "x" is the extension (in nm), "L" is the contour length
## (also in nm), "P" is the persistence length (in nm), where
## P has typically values around 1nm for single-stranded DNA
## The default temperature "Temp" is 25°C (298.15 K)
Force_From_Extension_WLC <- function(x,L,P,Temp=298.15)
{
 kBT <- 1.38e-02*Temp
 force <- (kBT/P)*(1/(4*(1-x/L)^2) - 1/4 +x/L)
 return(force)
}

## Extension_From_Force_WLC(f,L,P,Temp=298.15)
## code from E22 (Igor-Functions), 
## Felix Berkemeier, Johannes Stigler,..
Extension_From_Force_WLC <- function(f,L,P,Temp=298.15)
{
  kBT <- 1.38e-02*Temp;
  d <- (f*P)/(kBT);
  a <- -2.25 - d
  b <- 1.5 + 2.0*d
  c <- -d
  p <- b - (a^2)/3.0
  q <- c + (2.0*a^3-9*a*b)/27.0
  Dis <- complex(real=(q^2)/4.0+(p^3)/27.0,imaginary=0)
  u <- complex(real=0,imaginary=0)
  u <- (q/2.0 + sqrt(Dis))^(1.0/3.0)
  R <- complex(real=0,imaginary=0)
  R <- p/(3.0*u) - u - a/3.0
  return(Re(L*R))
}

## Integrated_WLC(x,L,P,Temp=298.15)
## calculates the elastic energy stored in a 
## wormlike chain with the parameters: contour length "L" (in nm),
## persistence length "P" (in nm) at a given Temperature "Temp" (in K,
## by default 25°C = 298.15K) and a given extension "x" (in nm)
Integrated_WLC <- function(x,L,P,Temp=298.15)
{
 kBT <- 1.38e-02*Temp
 a <- kBT/P
 res <- (a*x^2*(3*L-2*x))/(4*L*(L-x))
 return(res)
}

## Integrated_WLC2(x,L,P,Temp=298.15)
## see function "Energy_from_Extension_WLC(x,L,P,Temp=298.15)
## slightly different version, gives the identical result,
## haven't done any speed-testing so far
Integrated_WLC2 <- function(x,L,P,Temp=298.15)
{
 kBT <- 1.38e-02*Temp
 a <- kBT/P
 res <- a*(L/(4*(1-x/L)) - L/4 - ((1/4)*x)+(x^2)/(2*L))
 return(res)
}

## Force_From_Extension_EWLC(x,L,P,K,Temp=298.15)
## calculates the force at a given extension "x" for
## an extensible wormlike-chain with a contour lnegth "L" (in nm)
## a persistence length "P" (in nm), and the stretch modulus
## "K" (in nm) and a Temperature "Temp" (in K)
Force_From_Extension_EWLC <- function(x,L,P,K,Temp=298.15)
{
  kBT <- 1.38e-02*Temp;
  R <- x/L
  D <- (kBT/(P*K))
  a <- (8.0-8.0*R + 9.0*D - 12.0*R*D) / (4.0 + 4.0*D)
  b <- (4.0 - 8.0*R + 4.0*R^2 + 6.0*D -18.0*D*R +12.0*D*R^2) / (4.0 + 4.0*D)
  c <- D * (-6.0*R + 9.0*R^2 - 4.0*R^3) / (4.0 + 4.0*D)
  p <- b - a^2 / 3.0
  q = c + (2.0 * a^3 - 9.0 * a * b) / 27.0
  Dis <- complex(real=(q^2)/4.0 + (p^3)/27.0,imaginary=0)
  u <- complex(real=0,imaginary=0)
  u <- (-q/2.0 + sqrt(Dis))^(1.0/3.0)
  Force <- complex(real=0,imaginary=0)
  Force <- -p/(3.0*u) + u - a/3.0
  return(Re(Force*K))
}

## Extension_From_Force_EWLC(f,L,P,K,Temp=298.15)
## calculates the extension as a function of the force
## for an extensible wormlike chain (with a given contour length
## "L" (in nm), persistence length "P" (in nm), stretch modulus "K"
## (in nm) and Temperature "Temp" (in K), default value is 298.15
## (corresponding to 25°C)
Extension_From_Force_EWLC <- function(f,L,P,K,Temp=298.15)
{
  kBT <- 1.38e-02*Temp
  Z <- f/K
  D <- (kBT/(P*K))
  a <- (4.0*Z + 9.0*D + 12.0*D*Z) / (-4.0*D)
  b <- (-8.0*Z - 8.0*Z^2 - 6.0*D - 18.0*D*Z - 12.0*D*Z^2) / (-4.0*D)
  c <- (4.0*Z + 8.0*Z^2 + 4.0*Z^3 + 6.0*D*Z + 9.0*D*Z^2 + 4.0*D* Z^3) / (-4.0*D)
  p <- b - (a^2)/3.0
  q <- c + (2.0 * a^3 - 9.0 * a * b)/27.0
  Dis <- complex(real=(q^2)/4.0+(p^3)/27.0,imaginary=0)
  u <- complex(real=0,imaginary=0)
  u <- ((q/2.0) + sqrt(Dis))^(1.0/3.0)
  R <- complex(real=0,imaginary=0)
  R <- p/(3.0*u) - u - a/3.0
  return(Re(R*L))
}

## Extension_From_Force_WLC_Hook(f,L,P,k_eff,Temp=298.15)
## Returns the extension of a system of a series of a
## wormlike (characterized by "L","P" and "k_eff") chain 
## and a hookean spring ("k_eff")
Extension_From_Force_WLC_Hook <- function(f,L,P,k_eff,Temp=298.15)
{
  ext_hook <- f/k_eff
  ext_wlc <- Extension_From_Force_WLC(f,L,P,Temp=298.15)
  return(ext_hook+ext_wlc)
}

## Create_Extension_Array_From_Forces_WLC_Hook(forces,L,P,k_eff,Temp=298.15)
## creates and returns an array of extensions for a given array of forces
## "forces", for a system of a series of a wormlike chain (characterized by
## "L","P") and a hookean spring (characterized by "k_eff")
Create_Extension_Array_From_Forces_WLC_Hook <- function(forces,L,P,k_eff,Temp=298.15)
{
  extensions <- array(0,dim=c(length(forces)))
  for(i in 1:length(forces))
  {
    extensions[i] <- Extension_From_Force_WLC_Hook(forces[i],L,P,k_eff,Temp)
  }
  df_wlc_hook <- data.frame(forces=forces,extensions=extensions)
  return(df_wlc_hook)
}

## Test_Polym_Mechanics()
##
Test_Polym_Mechanics <- function()
{
  ## plot force from extension
  
  P <- 1.0
  L <- 40
  
  npnts_array <- 300
  
  
  forces1 <- array(0,dim=c(npnts_array))
  forces2 <- seq(from=0,to=100,length=npnts_array)
  forces3 <- array(0,dim=c(npnts_array))
  forces4 <- seq(from=0,to=100,length=npnts_array)
  
  extensions1 <- seq(from=0,to=0.9*L,length=npnts_array)
  extensions2 <- array(0,dim=c(npnts_array))
  extensions3 <- seq(from=0,to=0.9*L,length=npnts_array)
  extensions4 <- array(0,dim=c(npnts_array))
  K <- 1200
  k_eff <- 0.15
  for(i in 1:npnts_array)
  {
    forces1[i] <- Force_From_Extension_WLC(extensions1[i],L,P)
    extensions2[i] <- Extension_From_Force_WLC(forces2[i],L,P)
    forces3[i] <- Force_From_Extension_EWLC(extensions1[i],L,P,K)
    extensions4[i] <- Extension_From_Force_EWLC(forces2[i],L,P,K)
  }
  
  ## Plot for testing the WLC-functions
  plot(extensions1,forces1,type="l")
  points(extensions2,forces2,col="red")
  
  ## Plot for testing the EWLC-functions
  plot(extensions3,forces3,type="l")
  points(extensions4,forces4,col="red")
  
  ## test the WLC+hook (series) system
  forces_wlc_hook <- seq(from=0,to=30,length=npnts_array)
  df_wlc_hook <- Create_Extension_Array_From_Forces_WLC_Hook(forces_wlc_hook,L,P,k_eff,Temp=298.15)
  extensions5 <- seq(from=0,to=60,by=0.2)
  result <- approx(x=df_wlc_hook$extensions,y=df_wlc_hook$forces,xout=extensions5)
  plot(result$x,result$y)
  forces6 <- seq(from=0,to=6,by=0.05)
  extensions6 <- array(0,dim=c(length(forces6)))
  for(i in 1:length(forces6))
  {
    extensions6[i] <- forces6[i]/k_eff +   Extension_From_Force_WLC(forces6[i],L,P)
  }
  points(extensions6,forces6,col="red")
}
