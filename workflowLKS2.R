#  Workflow for Nmix models
#  Josh Nowak and Leah Swartz
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

#  Load dics
load("data/site_dic.RData")

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
################################################################################
#  Morph observation data
y_obs <- morph_data(raw_dat, site_dic) %>%
  left_join( .,cov_dat, by="site")# %>%
# filter(year == 2)

type_cov <- y_obs %>%
  ungroup() %>%
  select(site, impacted, created) %>%
  distinct %>%
  arrange(site) %>%
  select(impacted, created)

n_traps <- y_obs %>%
  ungroup() %>%
  select(site, year, prim, n_trap) %>%
  spread(prim, n_trap, -site, -year)

new_site <- as.numeric(as.factor(y_obs$site))

####################################################################################
##summarize data
####################################################################################
AMTI <- filter(y_obs, sp=="1") %>%
  group_by(SiteName, year)%>%
  summarise(nprim=max(prim),sum = sum(cnt),WetlandType=first(WetlandType))%>%
  group_by(WetlandType)%>%
  arrange(WetlandType,SiteName,sum)

BUBO <- filter(y_obs, sp=="2") %>%
  group_by(SiteName, year)%>%
  summarise(nprim=max(prim),sum = sum(cnt),WetlandType=first(WetlandType))%>%
  group_by(WetlandType)%>%
  arrange(WetlandType,SiteName,sum)

PSMA <- filter(y_obs, sp=="3") %>%
  group_by(SiteName, year)%>%
  summarise(nprim=max(prim),sum = sum(cnt),WetlandType=first(WetlandType))%>%
  group_by(WetlandType)%>%
  arrange(WetlandType,SiteName,sum)

RALU <- filter(y_obs, sp=="4") %>%
  group_by(SiteName, year)%>%
  summarise(nprim=max(prim),sum = sum(cnt),WetlandType=first(WetlandType))%>%
  group_by(WetlandType)%>%
  arrange(WetlandType,SiteName,sum)



#  Data as a list
jdat <- list(
  "nyear" = 2,
  "year"= y_obs$year,
  "species" = y_obs$sp,
  "nspecies" = length(unique(y_obs$sp)),
  "nsite" = length(unique(y_obs$site)),
  "nprim" = y_obs %>% 
    group_by(site) %>% 
    summarise(ntimes = max(prim)) %>% 
    .$ntimes,
  #"nsec" = 2,
  "nobs" = nrow(y_obs),
  "n_mu" = 1000,
  "n_tau" = 1/1000,
  "site" = new_site,
  "time" = y_obs$prim,
  #"occs" = y_obs$sec,
  "y" = y_obs$cnt,
  "impacted" = as.numeric(type_cov$impacted),
  "created" = as.numeric(type_cov$created)
)

N_init <- y_obs %>%
  group_by(sp, site) %>% 
  summarise(
    N = sum(cnt, na.rm = T)
  ) %>%
  .$N

n_mat <- matrix(N_init, nrow = jdat$nspecies, byrow = T)


nn <- array(NA, dim = c(jdat$nspecies, jdat$nsite, jdat$nyear, max(jdat$nprim)))
for(y in 1:2){
  for(i in 1:4){
    for(j in 1:jdat$nsite){
      for(k in 1:jdat$nprim[j]){
        nn[i,j,y,k] <- n_mat[i, j]
      }
    }
  }
}  


#  Initial Values
inits <- function(){
  list(
    Mean_n = rpois(jdat$nspecies, mean(y_obs$cnt)),
    Mean_p = runif(1, 0.02, 0.5),
    N = nn,
    w = array(1,dim=c(jdat$nyear,jdat$nspecies,jdat$nsite))
  )
}

#  Parameters to monitor
parms <- c("Mean_n", "Mean_p", "created_eff", "created_eff2","impacted_eff", "impacted_eff2", "gamma", "Phi", "N")

fit <- jags(
  jdat,
  inits, 
  parms,
  model.file = "models/time_corr_occ_abund2.txt",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 3000, 
  n.thin = 1
)

##mcmc plot
mcmcplots::mcmcplot(fit)
