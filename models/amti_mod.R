#  Workflow for Nmix models
#  Josh Nowak and Leah Swartz
#  02/2017
################################################################################
#  Packages
require(R2jags)
require(readr)
require(tidyr)
require(dplyr)
################################################################################
#  Set working directory
setwd("~/Leah_frogs")


#  Source helper functions
source("~/Leah_frogs/helpers/Nmix_utility_funs.R")

#  Load observation data
raw_dat <- read_csv("~/Blackrock/Data/TadpolePaper/2015AND2016TrappingData.csv")

#  Load covariate data
cov_dat <- read.csv("~/Blackrock/Data/Invert Paper/BR_site_covariates.csv") %>%
  mutate(
    site = site_dic$site_num[match(SiteName, site_dic$site_nm)]
  )%>%
  filter(site!="NA")%>%
  mutate(type=as.numeric(WetlandType))%>%
  mutate(reference=if_else(type=="3",1,1))%>%
  mutate(created=if_else(type=="1",1,0)) %>%
  mutate(impacted=if_else(type=="2",1,0))


#  Load dics
load("data/site_dic.RData")
################################################################################
#  Morph observation data
y_obs <- morph_data(raw_dat, site_dic) %>%
  left_join( .,cov_dat, by="site") 

AMTI <- filter(y_obs,sp=="1")%>%
  filter(year=="1")%>%
  select(sp,site,year, prim,sec,cnt)%>%
  spread(sec,cnt)%>%
  ungroup()%>%
  select(-sp)%>%
  select(-year)
write.csv(AMTI,file="AMTI.csv") ##had to add NA's in manually to make it the right size
amti <- read.csv("AMTI.csv",header=TRUE)

yy <- array(NA,dim=c(27,2,4))
for(k in 1:4){
  sel.rows <- amti$prim == k
  yy[,,k] <- as.matrix(amti)[sel.rows, 3:4]
}

sink("Nmix1.jags")
cat("
    model {
    
    # Priors
    omega ~ dunif(0, 1)
    for (k in 1:4){
    alpha.lam[k] ~ dnorm(0, 0.01)
    p[k] ~ dunif(0, 1)
    }
    
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:R){                          # Loop over R sites (95)
    z[i] ~ dbern(omega)                   # Latent suitability state (coin-flip for site suitability)
    for (k in 1:4){                       # Loop over primary periods (seasons)
    N[i,k] ~ dpois(lam.eff[i,k])       # Latent abundance state
    lam.eff[i,k] <- z[i] * lambda[i,k]
    log(lambda[i,k]) <- alpha.lam[k]
    # Observation model for replicated counts
    for (j in 1:T){                    # Loop over temporal reps (2)
    yy[i,j,k] ~ dbin(p[k], N[i,k])   # Detection
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j,k] <- p[k] * N[i,k]
    E[i,j,k] <- pow((yy[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j,k] ~ dbin(p[k], N[i,k])
    E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
    } #j
    } #k
    } #i
    
    # Derived and other quantities
    for (k in 1:4){
    totalN[k] <- sum(N[,k])	# Estimate total pop. size across all sites
    mean.abundance[k] <- exp(alpha.lam[k])
    }
    fit <- sum(E[,,])
    fit.new <- sum(E.new[,,])
    }
    ",fill = TRUE)
sink()

# Bundle data
R = nrow(yy)
T = ncol(yy)
win.data <- list(yy = yy, R = R, T = T)

# Initial values
# It is advisable to give initial values for the latent state z, the best option is to provide a vector of 1
Nst <- apply(yy, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(4, -1, 1), z = rep(1, 27))}

# Parameters monitored
params <- c("omega", "totalN", "alpha.lam", "p", "mean.abundance", "fit", "fit.new")

# MCMC settings
ni <- 3000
nt <- 15
nb <- 1500
nc <- 3

# Call JAGS from R (BRT 3 min)
out1 <- jags(win.data, inits, params, "Nmix1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3)

# Evaluation of fit
plot(out1$BUGSoutput$sims.list$fit, out1$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out1$BUGSoutput$sims.list$fit.new > out1$BUGSoutput$sims.list$fit)
mean(out1$BUGSoutput$mean$fit) / mean(out1$BUGSoutput$mean$fit.new)




# 12.3.3. Binomial-mixture model with overdispersion in both abundance and detection
# Specify model in BUGS language
sink("Nmix2.jags")
cat("
    model{
    # Priors
    for (k in 1:4){
    alpha.lam[k] ~ dnorm(0, 0.1)
    beta[k] ~ dnorm(0, 0.1)
    }
    
    # Abundance site and detection site-by-day random effects
    for (i in 1:R){
    eps[i] ~ dnorm(0, tau.lam)                    # Abundance noise
    }
    tau.lam <- 1 / (sd.lam * sd.lam)
    sd.lam ~ dunif(0, 3)
    tau.p <- 1 / (sd.p * sd.p)
    sd.p ~ dunif(0, 3)
    
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:R){                                 # Loop over R sites (95)
    for (k in 1:4){                              # Loop over days (7)
    N[i,k] ~ dpois(lambda[i,k])               # Abundance
    log(lambda[i,k]) <- alpha.lam[k] + eps[i]
    
    # Observation model for replicated counts
    for (j in 1:T){                           # Loop over temporal reps (2)
    yy[i,j,k] ~ dbin(p[i,j,k], N[i,k])      # Detection
    p[i,j,k] <- 1 / (1 + exp(-lp[i,j,k])) 
    lp[i,j,k] ~ dnorm(beta[k], tau.p) # random delta defined implicitly
    
    # Assess model fit using Chi-squared discrepancy
    # Compute fit statistic for observed data
    eval[i,j,k] <- p[i,j,k] * N[i,k]
    E[i,j,k] <- pow((yy[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
    # Generate replicate data and compute fit stats for them
    y.new[i,j,k] ~ dbin(p[i,j,k], N[i,k])
    E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k]+0.5)
    } #j
    ik.p[i,k] <- mean(p[i,,k])
    } #k
    } #i
    
    # Derived and other quantities
    for (k in 1:4){
    totalN[k] <- sum(N[,k])   # Estimate total pop. size across all sites
    mean.abundance[k] <- mean(lambda[,k])
    mean.N[k] <- mean(N[,k])
    mean.detection[k] <- mean(ik.p[,k])
    }
    fit <- sum(E[,,])
    fit.new <- sum(E.new[,,])
    }
    ",fill = TRUE)
sink()

# Bundle data
R = nrow(yy)
T = ncol(yy)
win.data <- list(yy = yy, R = R, T = T)

# Initial values
Nst <- apply(yy, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(4, -3, 3), beta = runif(4, -3, 3), sd.lam = runif(1, 0, 1), sd.p = runif(1, 0, 1))}

# Parameters monitored
params <- c("totalN", "alpha.lam", "beta", "sd.lam", "sd.p", "mean.abundance", "mean.N", "mean.detection", "fit", "fit.new")

# MCMC settings
ni <- 35000
nt <- 300
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 215 min)
out2 <- jags(win.data, inits, params, "Nmix2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Evaluation of fit
plot(out2$BUGSoutput$sims.list$fit, out2$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(50, 200), yylim = c(50, 200))
abline(0, 1, lwd = 2, col = "black")
mean(out2$BUGSoutput$sims.list$fit.new > out2$BUGSoutput$sims.list$fit)
mean(out2$BUGSoutput$mean$fit) / mean(out2$BUGSoutput$mean$fit.new)

# Summarize posteriors
print(out2, dig = 2)

max.day.count <- apply(yy, c(1, 3), max, na.rm = TRUE)
max.day.count[max.day.count == "-Inf"] <- NA
mean.max.count <- apply(max.day.count, 2, mean, na.rm = TRUE)
mean.max.count

par(mfrow = c(2, 1))
plot(1:4, mean.max.count, xlab = "Day", ylab = "Mean daily abundance", las = 1, ylim = c(0, 10), type = "b", main = "", frame.plot = FALSE, pch = 16, lwd = 2)
lines(1:4, out2$BUGSoutput$summary[25:31,5], type = "b", pch = 16, col = "blue", lwd = 2)
segments(1:4, out2$BUGSoutput$summary[25:31,3], 1:4, out2$BUGSoutput$summary[25:31,4], col = "blue")

plot(1:4, out2$BUGSoutput$summary[17:20,1], xlab = "Day", ylab = "Detection probability ", las = 1, ylim = c(0, 1), type = "b", col = "blue", pch = 16, frame.plot = FALSE, lwd = 2)
segments(1:4, out2$BUGSoutput$summary[17:20,3], 1:4, out2$BUGSoutput$summary[17:20,4], col = "blue")

out2$BUGSoutput$summary










yy <- array(NA,dim=c(28,2,4,2,4)) ##28 sites, 2 secondary occ., 4 prim. occ., 2 years, 4 species 

for(s in 1:4){ #for species in 1:4
  sel.rows <- y$sp==s
  for(k in 1:2){ #for year in 1:2
    sel.rows <- y$year==k
    for(p in 1:4){ #for prim. occ. in 1:4
    sel.rows <- y$prim==p
    yy[,,p,y,s] <- as.matrix(y)[sel.rows,5:6]
    }}}
##fill array with NA's in the right place
ar <- array(NA, dim=c(2,3,4))
dat <- data.frame(
  r=c(1,2,3),
  col=c(1,3,2),
  dept=c(4,2,3),
  val=1:3
)

for(i in 1:3){
  ar[dat[i,1],dat[i,2],dat[i,3]] <- dat[i,4]
}


  sel.rows <- y$prim==p
  yy[,,p,,] <- as.matrix(y)[sel.rows,5:6]
}#i


for(s in 1:2){
  sel.rows <- y_obs$sec==s
  yy[,s,,,]
}                                                
                                            (colnames(dyn.cov[,3:11])), 
                                                 (unique(as.character(dyn.cov$year)))))) ##28 sites, 2 secondary occ., 4 prim. occ., 2 years, 4 species 

Byy <- filter(y_obs,sp=="1")

