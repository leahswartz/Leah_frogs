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
    setwd("C:/Users/josh.nowak/Documents/GitHub/Leah_frogs")

    #  Source helper functions
    source("helpers/Nmix_utility_funs.R")

    #  Load data
    dat <- read_csv("C:/Users/josh.nowak/Documents/Leah_frogs/data.csv")
################################################################################
    #  Create lookup dictionary for sites
    uni_site <- unique(dat$Site)
    site_dic <- tibble(
      site_nm = uni_site,
      site_num = 1:length(uni_site)
    )
################################################################################
    #  Morph observation data
    y_obs <- dat %>%
      select(Site, TrappingOccasion, Year, contains("total")) %>%
      select(-TotalNumTrapsPredators, -TotalNumSweeps) %>%
      gather(
        Species, Count, -Site, -TrappingOccasion, -Year, -TotalNumTraps
      ) %>%
      transmute(
        sp = as.numeric(as.factor(Species)),
        site = site_dic$site_num[match(Site, site_dic$site_nm)],
        year = Year - min(Year) + 1,
        occ = TrappingOccasion,
        prim = id_primary(TrappingOccasion),
        cnt = Count,
        n_trap = TotalNumTraps 
      ) %>%
      group_by(sp, site, year, prim) %>%
      mutate(
        sec = 1:n()
      ) %>%
      select(sp:prim, sec, cnt, n_trap)
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
################################################################################
    #  TODO
    #   multi-species, optional
    #   random effect on site
    #   accommodate primary and secondary occassions
    #   add covariates
