
functions {
  // real ForceWLC(real P,real L,real Temp,real x)
  // calculates the force of a WLC-model (Marko & Siggia interpolation formula)
  // for an extension of "x", the WLC is characterized by
  // its persistence lengh P (nm)
  // its contour lenght L (nm)
  // and the Temperature (K)
  real ForceWLC(real P,real L,real Temp,real x) // calculate the force-extension for a WLC-model
  {
    real kB;
    real a;
    real force;
    kB <- 1.38e-02;
    a <- (kB*Temp/P);
    force <- a*(1.0/(4.0*(1-x/L)^2) - 0.25 + x/L);
    return force;
  }
  
  // real IntegratedWLC(real P,real L, real Temp,real x)
  // calculates the Integral of ForceWLC(...x') from x'=0 to
  // x' = x
  real IntegratedWLC(real P,real L, real Temp,real x)
  {
    real kB; 
    real a;
    real res;
    kB <- 1.38e-02;
    a <- (kB*Temp)/P;
    res <- a*(L/(4.0*(1-x/L))-L/4.0 - (0.25*x)+(x^2)/(2*L));
    return res;
  }
  
  // real[] NormalToPolar(real[] normal)
  // expects a complex number in normal form (z = [a,b], with z=a+ib)
  // and returns the polar form: z=[r,phi],with z= r*exp(i*phi)
  real[] NormalToPolar(real[] normal)
  {
    real polar[2];
    polar[1] <-sqrt(normal[1]^2+normal[2]^2);
    polar[2] <- atan2(normal[2],normal[1]);
    return polar;
  }
  
  // real[] PolarToNormal(real[] polar)
  // input: a complex number in polar form
  // returns the normal form
  real[] PolarToNormal(real[] polar)
  {
    real normal[2];
    normal[1] <- polar[1]*cos(polar[2]);
    normal[2] <- polar[1]*sin(polar[2]);
    return normal;
  }
  
  // real[] ReturnFirstNthRoot(real[] polar,real n)
  // input: a complex number in polar form z [r,phi], with z = r*exp(i*phi)
  // and n (determining which root shall be returned)
  // output/return value: "first" complex root: z^(1/n)_1 = r^(1/n)*exp(i*phi/n)
  // in polar form
  real[] ReturnFirstNthRoot(real[] polar,real n)
  {
    real polar_root[2];
    polar_root[1] <- (polar[1])^(1/n);
    polar_root[2] <- (polar[2]/n);
    return polar_root;
  }
  
  // real[] InverseFromNormalForm(real[] normal)
  // input: a complex number in normal form (z [a,b], with z=a+b*i)
  // output: the inverse z^(-1), in normal form...
  real[] InverseFromNormalForm(real[] normal)
  {
    real normal_inverse[2];
    real r;
    real theta;
    r <- sqrt(normal[1]^2+normal[2]^2);
    theta <- atan2(normal[2],normal[1]);
    normal_inverse[1] <- (1/r)*cos(theta);
    normal_inverse[2] <- -(1/r)*sin(theta);
    return normal_inverse;
  }
  
  // real ExtensionFromForceWLC(real P,real L, real Temp,real Force)
  // input: persistence length P(in nm), contour length L(in nm),
  // temperature Temp (in K) and the force Force (in pN)
  // output: returns the extension of a WLC-model described by the
  // above parameters at the given force
  real ExtensionFromForceWLC(real P,real L, real Temp, real Force)
  {
    real kB;
    real a;
    real b;
    real c;
    real d;
    real p;
    real q;
    real Dis;
    real result;
    real polar[2]; // Dis in polar form
    real normal[2];
    real Dis_root[2];
    real second_root[2];
    real norm_inv[2];
    
    kB <- 1.38e-02;
    d <- (Force*P)/(kB*Temp);
    
    a <- -2.25 - d;
    b <- 1.5+2.0*d;
    c <- -d;
    
    p <- b - (a^2)/(3.0);
    q <- c+ (2*a^3 - 9*a*b)/27.0;
    Dis <- (q^2)/4.0 + (p^3)/27.0;
    normal[1] <- Dis;
    normal[2] <- 0;
    polar <- NormalToPolar(normal);
    Dis_root <- ReturnFirstNthRoot(polar,2.0);
    
    normal <- PolarToNormal(Dis_root);
    normal[1] <- normal[1] + q/2.0;
    polar <- NormalToPolar(normal);
    second_root <- ReturnFirstNthRoot(polar,3.0);
    normal <- PolarToNormal(second_root);
    norm_inv <- InverseFromNormalForm(normal);
    
    result <- (p/3.0)*norm_inv[1] - normal[1] - a /3.0;
    return (result*L);
  }
  
  // real ExtensionFromForceWLCHook(real P, real L, real Temp, real k_eff, real Force)
  // input: persistence-length P (nm), contour-length L (nm), Temperature Temp (K),
  // effective spring constant of both traps (in pN/nm) and a force value "Force"
  // returns the extension of a series of WLC and hookean spring at the given force
  real ExtensionFromForceWLCHook(real P,real L, real Temp,real k_eff, real Force)
  {
    real extension_hook;
    real extension_wlc;
    real extension_total;
    extension_hook <- Force/k_eff;
    extension_wlc <- ExtensionFromForceWLC(P,L,Temp,Force);
    extension_total <- extension_wlc + extension_hook;
    return extension_total;
  }
  
  
  // int BinarySearch(real ext_value,real[] extension_array)
  /*int BinarySearch(real ext_value,real[] extensions_array,int low,int high);
  int BinarySearch(real ext_value,real[] extensions_array,int low,int high)
  {
    int mid;
    
    mid <- low+(high-low)/2;
    
    if((extensions_array[mid]<=ext_value)&&(extensions_array[mid+1]>ext_value))
      return mid;
    if(ext_value>extensions_array[mid])
      return BinarySearch(ext_value,extensions_array,mid+1,high);
    else if(ext_value<extensions_array[mid])
    {
      return BinarySearch(ext_value,extensions_array,low,mid);
    }
    return mid;
    
  }*/
    
  // int Search(real ext_value,real[] extensions_array,int max_pos)
  // very simple Search-routine that looks for the place where
  // "ext_value" would be inserted in "extensions_array"
  // important: "extensions_array[] " has to be sorted, 
  // max_pos is the length of the vector
  int Search(real ext_value,real[] extensions_array,int max_pos)
  {
    //   print("first value in array: ");
    //  print(extensions_array[1]);
    //  print("last value in array: ");
    //  print(extensions_array[max_pos]);
    
    for(j in 1:(max_pos-1))
    {
      if((ext_value>extensions_array[j])&&(ext_value<=extensions_array[j+1]))
      {
        //    print("Position found: ");
        //    print(j);
        return j;
      }
    }
    //}
  return +1;
}

// real LinearInterpolate(int pos,real ext_value,real[] extensions,real[] forces)
// input: int pos (positions, where the value "ext_value" would fit in the ordered
// "extensions" - vector), "ext_value" is the extension, "extensions[]" is a vector
// containing extensions that belong to the corresponding forces in "forces[]"
// output: returns the (linearly) interpolated force at the extension "ext_value"
real LinearInterpolate(int pos,real ext_value,real[] extensions,real[] forces)
{
  real dF;
  real dx_part;
  real dx_total;
  real interp_force;
  dF <- forces[pos+1] - forces[pos];
  dx_part <- ext_value - extensions[pos];
  dx_total <- extensions[pos+1] - extensions[pos];
  interp_force <- forces[pos] + dF*dx_part/dx_total;
  return interp_force;
}

// real ForceFromExtensionWLCHook(real P,real L,real Temp,real k_eff,real extension)
// input: persistence-length P (nm), contour-length L (nm), Temperature Temp (K),
// effective spring constant of the traps (pN/nm) and an extension value
// computes the force for the given extension for a system of WLC and hookean spring
// in series
real ForceFromExtensionWLCHook(real P,real L, real Temp, real k_eff, real extension)
{
  // create array of forces and extensions for this system...
  int npnts;
  real forces_local[500];
  real extensions_local[500];
  real dF;
  real F_min;
  real F_max;
  int pos;
  real interp_force;
  
  F_min <- 0.01;
  F_max <- 9.0;
  dF <- (F_max-F_min)/(500);
  
  for(i in 1:500)
  {
    forces_local[i] <- F_min+(i-1)*dF;
    extensions_local[i] <- ExtensionFromForceWLCHook(P,L,Temp,k_eff,forces_local[i]);
  }
  // now, search the correction "insertion position" for the extension
  print("Extension: ");
  print(extension);
  pos <- Search(extension,extensions_local,500);
  print("Position: ");
  print(pos);
  // finally, interpolate the force....
  interp_force <- LinearInterpolate(pos,extension,extensions_local,forces_local);
  return interp_force;
}

// real OccupationProbabilityFromForceOpen(real P,real L,real Temp,real k_eff, real G_int,real dz,real F_open)
// input: persistence-length P (nm), contour-length L (nm), Temperature Temp (K), effectives spring constant
// of both traps (in pN/nm), the free energy of the interaction (in pN nm) and half of the end-to-end distance
// in the closed state and a force value "F_open"
// output: returns the occupation probability of the open state at the force "F_open"
real OccupationProbabilityFromForceOpen(real P,real L,real Temp,real k_eff, real G_int,real dz,real F_open)
{
  real G_open_trap;
  real G_open_leash;
  real G_closed_trap;
  real G_closed_leash;
  real G_open;
  real G_closed;
  real kB;
  real x_open_leash;
  real x_closed_trap;
  real F_closed;
  real P_open;
  real exp_open;
  real exp_closed;
  kB <- 1.38e-02;
  G_open_trap <- (0.5*F_open^2)/(k_eff);
  x_open_leash <- ExtensionFromForceWLC(P,L,Temp,F_open);
  G_open_leash <- IntegratedWLC(P,L,Temp,x_open_leash);
  G_open <- G_open_leash+G_open_trap;
  
  x_closed_trap <- x_open_leash + F_open/k_eff-2*dz; // distance in the closed conformations is shorter by 2*dz...
  F_closed <- k_eff*x_closed_trap;
  G_closed_trap <- (0.5*F_closed^2)/(k_eff);
  G_closed_leash <- IntegratedWLC(P,L,Temp,2*dz);
  G_closed <- G_closed_trap + G_closed_leash + G_int;
  exp_open <- exp(-G_open/(kB*Temp));
  exp_closed <- exp(-G_closed/(kB*Temp));
  P_open <- exp_open/(exp_open+exp_closed);
  return P_open;
}

// real OccupationProbabilityFromForceClosed(real P,real L,real Temp,real k_eff, real G_int,real dz,real F_closed)
// input: persistence-length P (nm), contour-length L (nm), Temperature Temp (K), effectives spring constant
// of both traps (in pN/nm), the free energy of the interaction (in pN nm) and half of the end-to-end distance
// in the closed state and a force value "F_closed"
// output: returns the occupation probability of the closed state at the closed state equilibrium force "F_closed"
real OccupationProbabilityFromForceClosed(real P,real L,real Temp,real k_eff,real G_int,real dz, real F_closed)
{
  real G_closed_trap;
  real G_closed_leash;
  
  real G_closed;
  real G_open_trap;
  real G_open_leash;
  real x_open;
  real x_open_leash;
  real F_open;
  real G_open;
  real exp_open;
  real exp_closed;
  real P_closed;
  real kB;
  kB <- 1.38e-02;
  
  G_closed_trap <- (0.5*F_closed^2)/(k_eff);
  G_closed_leash <- IntegratedWLC(P,L,Temp,2*dz);
  
  G_closed <- G_closed_trap + G_closed_leash + G_int;
  
  x_open <- F_closed/k_eff+2*dz; // effective distance of the WLC+ spring system is largen than only the spring extension in the closed version....
  F_open <- ForceFromExtensionWLCHook(P,L,Temp,k_eff,x_open);
  print("F_open: ");
  print(F_open);
  x_open_leash <- ExtensionFromForceWLC(P,L,Temp,F_open);
  
  G_open_trap <- (0.5*F_open^2)/(k_eff);
  G_open_leash <- IntegratedWLC(P,L,Temp,x_open_leash);
  G_open <- G_open_trap + G_open_leash;
  exp_open <- exp(-G_open/(kB*Temp));
  exp_closed <- exp(-G_closed/(kB*Temp));
  P_closed <- exp_closed/(exp_open+exp_closed);
  return P_closed;
}
  
  real LifetimeFromForce_BS_Open(real tau0_open, real P,real L_trans,real L_open,real Temp,real k_eff,real F_open)
  {
    // calculate the energy in the open state:
    real G_open_leash;
    real G_open_trap;
    real G_open;
    real x_open_leash;
    real x_total_trans;
  
    real G_trans_leash;
    real G_trans_trap;
    real G_trans;
    real x_trans_leash;
    real F_trans;
    real tau_open;
    real kB;
    kB <- 1.38e-02;
  
    G_open_trap <- (0.5*F_open^2)/(k_eff);
  
    // equilibrium extension of the leash:
    x_open_leash <- ExtensionFromForceWLC(P,L_open,Temp,F_open);
    G_open_leash <- IntegratedWLC(P,L_open,Temp,x_open_leash);
    G_open <- G_open_leash+G_open_trap;
  
    // calculate the corresponding force
    // in the transition state:
    x_total_trans <- x_open_leash + F_open/k_eff;
    F_trans <- ForceFromExtensionWLCHook(P,L_trans,Temp,k_eff,x_total_trans);
    x_trans_leash <- ExtensionFromForceWLC(P,L_trans,Temp,F_trans);
    G_trans_leash <- IntegratedWLC(P,L_trans,Temp,x_trans_leash);
  
    G_trans_trap <- (0.5*F_trans^2)/(k_eff);
    G_trans <- G_trans_leash + G_trans_trap;
  
    tau_open <- tau0_open*exp(-(G_open - G_trans)/(kB*Temp));
    return tau_open;
  }
  
  real LifetimeFromForce_BS_Closed(real tau0_closed, real P,real L_trans,real L_open,real Temp,real k_eff,real F_closed)
  {
    //real G_closed_leash;
    real G_closed_trap;
    real G_closed;
    
    real x_trans_total;
    real F_trans;
    
    real x_trans_leash;
    real G_trans_trap;
    real G_trans_leash;
    real G_trans;
    
    real tau_closed;
    real kB;
    kB <- 1.38e-02;
    
    G_closed_trap <- (0.5*F_closed^2)/(kB*Temp);
    G_closed <- G_closed_trap;
    
    // 
    x_trans_total <- F_closed/k_eff;
    F_trans <- ForceFromExtensionWLCHook(P,L_trans,Temp,k_eff,x_trans_total);
    G_trans_trap <- (0.5*F_trans^2)/(kB*Temp);
    x_trans_leash <- ExtensionFromForceWLC(P,L_trans,Temp,F_trans);
    G_trans_leash <- IntegratedWLC(P,L_trans,Temp,x_trans_leash);
    G_trans <- G_trans_trap + G_trans_leash;
    
    tau_closed <- tau0_closed*exp(-(G_closed-G_trans)/(kB*Temp));
    return tau_closed;
  }
  
  // real ExtensionGainFromForceOpen(real P, real L, real Temp, real dz, real F_open)
  // input: persistence-length [nm], contour-length L [nm], Temperature Temp [K], half of the
  // end-to-end distance of the leash in the closed state dz [nm] and the force in the
  // open state F_open [nm]
  // returns the extension gain of our system at this particular value of the force
  // in the open state
  real ExtensionGainFromForceOpen(real P,real L,real Temp,real dz,real F_open)
  {
    real x_leash_open;
    real dx;
    // end-to-end distance of the WLC...
    x_leash_open <- ExtensionFromForceWLC(P,L,Temp,F_open);
    dx <-x_leash_open - 2*dz;
    return(dx);
  }
  
  // real ExtensionGainFromForceClosed(real P, real L, real Temp, real dz, real k_eff, real F_closed)
  // input: persistence length P [nm], contour-length L [nm], Temperature Temp [K], half of the
  // end-to-end distance of the leash in the closed state dz [nm], the effective trap-stiffness
  // [pN/nm] and the force in the closed state
  // returns the extension gain of our system at this particular force-value
  real ExtensionGainFromForceClosed(real P,real L,real Temp,real dz,real k_eff,real F_closed)
  {
     real x_open;
     real F_open;
     real x_leash_open;
     real dx;
     x_open <- F_closed/k_eff+2*dz;
     F_open <- ForceFromExtensionWLCHook(P,L,Temp,k_eff,x_open);
     x_leash_open <- ExtensionFromForceWLC(P,L,Temp,F_open);
     dx <- x_leash_open - 2*dz;
     return(dx);
  }
}

data 
{
  int<lower=1> N_open; // number of data points (forces_open)
  real forces_open[N_open];
  real<lower=0> P;
  real<lower=0> CL; // contour length[nm]
  real<lower=0> k_eff; // the effective spring constant
  real G_int; // interaction energy (in pN*nm)
  real<lower=0> dz;
  real<lower=0> Temp;
  real<lower=0> sigma; // noise level (has to be improved, noise not normally distributed)
}


model 
{
}

generated quantities
{
  real occ_probs_open_nf[N_open];
  real occ_probs_open[N_open];
  real mu;
  real alpha_par;
  real beta_par;
  for(i in 1:N_open)
  {
    occ_probs_open_nf[i] <- OccupationProbabilityFromForceOpen(P,CL,Temp,k_eff,G_int,dz,forces_open[i]);
    mu <- occ_probs_open_nf[i];
    // transform mu and sigma to alpha and beta
    alpha_par <- ((1.0-mu)/(sigma^2)-1.0/mu)*mu^2;
    beta_par <- alpha_par*(1.0/mu - 1.0);
    occ_probs_open[i] <- beta_rng(alpha_par,beta_par);
  }
}