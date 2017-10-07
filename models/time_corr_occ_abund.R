model{
  #  Naming 
  #  Parameter names begin with a capitalized letter
  #  Data are all lower case
  
  #  Priors - try different priors
  #  Mean abundance
  for(i in 1:nspecies){
    Mean_n[i] ~ dnorm(n_mu, n_tau)T(0,)
  }
  #Mean_n ~ dgamma(0.001, 0.001)
  
  #  Mean detection
  Mean_p ~ dnorm(0, 0.35)T(-10, 10)
  
  #trap covariate
  n_trap_eff ~ dnorm(0, 0.01)
  
  #covariate - wetland type prior
  # type_eff ~ dunif(-10,10)
  # Type_eff ~ dnorm(0,0.01)
  #for(i in 1:3){
  #Type_eff[i] ~ dnorm(0,0.01)
  #}
  for(i in 1:nspecies){
    created_eff[i] ~ dnorm(0,0.01)
    impacted_eff[i] ~ dnorm(0,0.01)
  }
  
  occ_eff ~ dnorm(0, 0.35)
  
  #  Random effect prior
  for(sp in 1:nspecies){
    Tau_site[sp] <- 1/(Sd_site[sp]^2)
    Sd_site[sp] ~ dunif(0, 200)
    for(s in 1:nsite){
      Site_eff[sp,s] ~ dnorm(0, Tau_site[sp])
    }
  }    
  
  #  Linear predictors
  for(y in 1:nyear){
   for(sp in 1:nspecies){
    for(s in 1:nsite){
      Lambda[sp, s] <- exp(Mean_n[sp] + Site_eff[sp,s] + impacted_eff[sp] * impacted[s] + created_eff[sp] * created[s]) 
      for(t in 1:nprim[s]){
        for(occ in 1:nsec){
          logit(P[sp, s, t, occ]) <- Mean_p# + n_trap_eff * n_trap[s, t]
          }
        }
      }
    }
  }
  #  Should never have to change code below this point
  #  True latent abundance
  for(y in 1:nyear){
    for(sp in 1:nspecies){
     for(s in 1:nsite){
      N[sp,s,1] ~ dpois(Lambda[sp,s])
    }
   }
  }
  
  for(s in 1:nsite){
    gamma[s] ~ dbeta(1, 1)
  }
  
 for(y in 1:nyear) {
    for(sp in 1:nspecies){
      for(s in 1:nsite){
       for(t in 2:nprim[s]){
          N[sp,s,t] ~ dbinom(gamma[s], N[sp, s, t-1])
      }
    }
  }
}
  #  Observation process
  for(i in 1:nobs){
    y[i] ~ dbinom(P[species[i], site[i], time[i], occs[i]], N[species[i], site[i], time[i]])
  }
  
  
  
}

##add loop for year before species
##might need to change things for 26 vs 28 sites
##add occupancy line
