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
      left_join( .,cov_dat, by="site") %>%
      filter(year == 2)
    
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
    
    #  Data as a list
    jdat <- list(
      "species" = y_obs$sp,
      "nspecies" = length(unique(y_obs$sp)),
      "nsite" = length(unique(y_obs$site)),
      "ntimes" = y_obs %>% 
        group_by(site) %>% 
        summarise(ntimes = max(occ)) %>% 
        .$ntimes,
      "nsec" = 2,
      "nobs" = nrow(y_obs),
      "n_mu" = 10,
      "n_tau" = 1/100,
      "site" = new_site,
      "time" = y_obs$prim,
      "occs" = y_obs$sec,
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
    nn <- (matrix(N_init, nrow = 4) + 2) * 100
    
    #  Initial Values
    inits <- function(){
      list(
        Mean_n = rpois(1, mean(y_obs$cnt)),
        Mean_p = runif(1, 0.02, 0.5),
        N = nn
      )
    }
    
    #  Parameters to monitor
    parms <- c("Lambda", "Mean_n", "Mean_p", "P", "N", "Site_eff", "created_eff", "impacted_eff")
    
    fit <- jags(
      jdat,
      inits, 
      parms,
      model.file = "models/Nmix_sNtN_trapD.txt",
      n.chains = 3,
      n.iter = 500,
      n.burnin = 100, 
      n.thin = 1
    )
    
    
    
    
    
    
    
    
    
    
################################################################################
    #  Call model on grouped data, species by year
    #  Remove species grouping if using multi-species model
    #  This example call shows how to use a single covariate
    fit1 <- y_obs %>%
      filter(site != 22) %>%
      group_by(sp, year) %>%
      do(fit = 
        call_jags(
          x = .,
          covs = "n_trap",
          model.file = "models/Nmix_cN_trapD.txt",
          n.chains = 3,
          n.iter = 500,
          n.burnin = 100, 
          n.thin = 1
        )
      )

    #  Example call without covariates
    fit2 <- y_obs %>%
      group_by(sp, year) %>%
      do(fit = 
        try(call_jags(
          x = .,
          covs = NULL,
          model.file = "models/Nmix_cN_cD.txt",
          n.chains = 3,
          n.iter = 500,
          n.burnin = 100, 
          n.thin = 1
        ))
      )
      
    #  Fit 3=random effect of site, trap effect on detection
    
    fit3 <- y_obs %>%
      group_by(sp, year) %>%
      do(fit = 
        try(call_jags(
          x = .,
          covs = c("n_trap"),
          model.file = "models/Nmix_sN_trapD.txt",
          n.chains = 3,
          n.iter = 500,
          n.burnin = 100, 
          n.thin = 1
        ))
      )
    # Fit 4 = random effect of site, trap effect on detection, wetland type effect on abundance
    fit4 <- y_obs %>%
      group_by(sp, year) %>%
      do(fit = 
           try(call_jags(
             x = .,
             covs = c("n_trap","impacted","created"),
             model.file = "models/Nmix_sNtN_trapD.txt",
             n.chains = 3,
             n.iter = 500,
             n.burnin = 100, 
             n.thin = 1
           ))
      )
    # Fit 5 =  trap effect on detection, wetland type effect on abundance
    fit5 <- y_obs %>%
      group_by(sp, year) %>%
      do(fit = 
           try(call_jags(
             x = .,
             covs = c("n_trap","impacted","created"),
             model.file = "models/Nmix_tN_trapD.txt",
             n.chains = 3,
             n.iter = 500,
             n.burnin = 100, 
             n.thin = 1
           ))
      )
    #  
################################################################################
    #  TODO
    #   multi-species, optional
    #   random effect on site
    #   accommodate primary and secondary occassions
    #   add covariates
    ##look at zero inflation
    ##add occupancy part to account for zeros
    ##mcmc plots
