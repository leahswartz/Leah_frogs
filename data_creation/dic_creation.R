    #  Dictionary/lookup creation for Leah frogs
    #  Josh Nowak
    #  03/2017
################################################################################
    require(readr)
    require(dplyr)
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
    
    save(site_dic, file = "data/site_dic.RData")
################################################################################
    #  End