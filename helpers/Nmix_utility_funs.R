    id_primary <- function(x){
    
      tmp <- cut(x, breaks = c(0, 2, 4, 6, 8), labels = 1:4)
      out <- as.integer(tmp)

    return(out)
    }
################################################################################
    init_N <- function(x){

      out <- x %>% 
        group_by(sp, site, year) %>% 
        summarise(
          N = max(cnt, na.rm = T)
        ) %>%
        group_by(sp, site) %>%
        mutate(
          N = replace(N, N == 0, round(mean(N, na.rm = T) + 2))
        ) %>%
        .$N

    return(out)
    }
################################################################################
    pretty_cov <- function(x, covar_nm){
      out <- x %>%
        ungroup() %>%
        select_("site", "occ", covar_nm) %>%
        spread_("occ", covar_nm) %>%
        select(-1)
    
    return(out)
    }
################################################################################
    call_jags <- function(x, covs = list(), ...){

      #  Initial Values
      inits <- function(){
        list(
          Mean_n = rpois(1, mean(x$cnt)),
          Mean_p = runif(1, 0.02, 0.5),
          N = init_N(x)
        )
      }

      #  Parameters to monitor
      parms <- c("Lambda", "Mean_n", "Mean_p", "P", "N", paste0(covs, "_eff"))
      
      #  Rename sites to catch missed visits
      new_site <- as.numeric(as.factor(x$site))
      
      #  Gather covariates in a named list
      if(is.null(covs)){
        cov_vals <- list()
      }else{
        cov_vals <- lapply(covs, function(y){
          pretty_cov(x, y)
        })
        names(cov_vals) <- covs
      }
      
      #  Data as a list
      jdat <- c(cov_vals, list(
        "nsite" = length(unique(x$site)),
        "ntimes" = x %>% 
          group_by(site) %>% 
          summarise(ntimes = max(occ)) %>% 
          .$ntimes,
        "nobs" = nrow(x),
        "n_mu" = 10,
        "n_tau" = 1/100,
        "site" = new_site,
        "time" = x$occ,
        "y" = x$cnt
      ))
      
      out <- jags(
        jdat, 
        inits,
        parms,
        ...
      )
      
      out$site_dic <- tibble(
        original = x$site,
        new = new_site
      )
      
    return(out)
    }
################################################################################
    #  End