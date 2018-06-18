#' Run linear mixed models to identify predictors of species' eigencentralities
#' across sites
#'
#'@importFrom magrittr %>%
#'@importFrom parallel makePSOCKcluster setDefaultCluster clusterExport stopCluster clusterEvalQ detectCores parLapply
#'
#'@param mods_list A \code{list} of models returned from \code{hurdleModel}
#'@param n_bootstraps \code{integer} representing the number of bootstrap replicates
#'to perform
#'@param n_cores Positive integer stating the number of processing cores to split the job across.
#'Default is \code{parallel::detect_cores() - 1}
#'
#'@details Species' network eigencentralities are used as the outcome variable in a linear mixed model
#'to identify environmental and biotic predictors. Only species with mean eigencentrality scores
#'in the top 25th percentile are considered, as these 'keystone' species are the ones that have the
#'strongest influences on community structure. The model pipeline is as follows:\cr
#'1: shuffle the dataset, impute any missing species' trait values and run a linear mixed model with
#'species and region as random effects (variable intercepts).
#'2: run a full model with all possible covariates as predictors
#'3: use stepwise variable selection (using AIC as the loss function) to identify key predictors
#'4: repeat \code{n_bootstraps} times to account for uncertainty in coefficients\cr
#'Models are then run to identify predictors of site-level interaction \emp{B}'os using a similar
#'pipeline
#'
#'@export
#'
networkModel = function(mods_list, n_bootstraps, n_cores){

  if(missing(n_cores)){
      n_cores <- parallel::detectCores() - 1
  }

  #### Function to impute missing centrality variables and randomly remove 25% of observations ####
  shuffle_rows <- function(empty){
    cent_mod_data %>%
      dplyr::group_by(species) %>%

      ## Impute missing nesting, habitat and clutch size information
      dplyr::mutate(Ground.nest = ifelse(is.na(Ground.nest),
                                         sample(c(0, 1), 1), Ground.nest),
                    Disturbed.hab = ifelse(is.na(Disturbed.hab),
                                           sample(c(0, 1), 1), Disturbed.hab),
                    Clutch.size = ifelse(is.na(Clutch.size),
                                         sample(seq(-1, 1, length.out = 100), 1),
                                         Clutch.size)) %>%
      dplyr::ungroup() %>%
      dplyr::sample_n(., nrow(cent_mod_data) * .75, FALSE)
  }

  #### Function to shuffle community-level beta diversity data and randomly remove 25% of observations ####
  shuffle_beta_rows <- function(empty){
    betadiv_mod_data %>%
      dplyr::sample_n(., nrow(betadiv_mod_data) * .75, FALSE)
  }

  #### Extract names of regions within the flyway and load species' traits ####
  region_names <- names(mods_list)
  data("Bird.traits")

  #### Determine which environmental covariates to use as predictors of
  # centrality and network beta diversity ####
  mod_covs <- lapply(seq_along(mods_list), function(x){
    cov_names <- names(mods_list[[x]]$network_metrics)
    cov_names <- cov_names[-which(cov_names %in% Bird.traits$species)]

  })
  all_covs <- unique(unlist(mod_covs))

  ## Only keep predictors that were retained in all regional models
  to_keep <- lapply(seq_along(all_covs), function(i){
    any(grepl(all_covs[i], mod_covs, fixed = TRUE) == FALSE)
  })

  covs_keep <- paste(c(all_covs[!unlist(to_keep)],
                       'species', 'centrality$'),
                        collapse = "$|^")
  covs_keep <- paste('^' , covs_keep, sep = '')

  #### Convert eigencentrality scores to long format, keeping only scores for
  # species with above-average centrality (top 25th percentile) ####
  cents <- lapply(seq_along(mods_list), function(x){

    lassoAbund <- mods_list[[x]]$abundance_mod
    sp_names <- rownames(lassoAbund$graph)

    ## Find species with high overall centrality (keystone species)
    abundance_cent <- mods_list[[x]]$network_metrics
    high_cent <- as.numeric(quantile(colMeans(abundance_cent[ , sp_names]),
                                     probs = 0.75))
    keystone_sp <- names(abundance_cent[, sp_names][which(colMeans(abundance_cent[ , sp_names])
                                                          >= high_cent)])

    ## For keystone species, convert centrality scores to long format
    cent_mod_data = abundance_cent %>%
      tidyr::gather('species' = sp_names, species, centrality) %>%
      dplyr::filter(centrality > 0) %>%
      dplyr::filter(species %in% keystone_sp) %>%

      ## logit-transform centrality, as it represents skewed proportional data
      dplyr::mutate(centrality = log(centrality / (1 - centrality))) %>%
      dplyr::filter(is.finite(centrality)) %>%
      dplyr::select(matches(covs_keep)) %>%
      dplyr::mutate(region = region_names[x]) %>%
      dplyr::left_join(Bird.traits) %>%

      ## Floating nest in this dataset is rank-deficient (all are NA's or zero)
      dplyr::select(-Floating.nest) %>%

      ## Many habitat and nesting traits are correlated, keep only
      # data on ground nesting and use of urbanised habitats as predictors
      dplyr::select(-Lower.canopy.nest, -Upper.canopy.nest,
                    -Forest.hab, -Grassland.hab)
  })

  ## Bind the regional centrality datasets together to create the species-level datast
  cent_mod_data = dplyr::bind_rows(cents) %>%
    dplyr::select(-beta_os)

  ## Now create the community-level beta diversity dataset
  betadiv_mod_data = dplyr::bind_rows(cents) %>%
    dplyr::select(-species,-centrality,-Clutch.size,-Ground.nest,-Disturbed.hab,
                  -Migrate.status,-Rarity,-Diet.diversity,-Hab.diversity) %>%
    dplyr::mutate(beta_os = as.vector(scale(beta_os)))

  #### Use the covariate names to build a linear mixed model formula for the species data ####
  fixed_effects <- paste(colnames(cent_mod_data[, !names(cent_mod_data) %in%
                                                  c('species','centrality',
                                                    'region')]),
                         collapse = '+')

  ## Model formulate should contain all fixed effects, as well as variable intercepts
  # for 'species' and 'region'
  full_formula <- paste(paste('centrality', '~', fixed_effects),
                        '(1 | species)', '(1 | region)', sep = " + ")

  #### Create a list of shuffled and imputed datasets to run models across ####
  booted_list <- vector('list', n_bootstraps)
  booted_datas <- lapply(booted_list, shuffle_rows)

  #### Create the formula and shuffled datasets for the community-level beta diversity data ####
  beta_fixed_effects <- paste(colnames(betadiv_mod_data[, !names(betadiv_mod_data) %in% c('beta_os',
                                                                                       'region')]),
                         collapse = '+')

  beta_full_formula <- paste(paste('beta_os', '~', beta_fixed_effects),
                             '(1 | region)', sep = " + ")

  beta_booted_list <- vector('list', n_bootstraps)
  beta_booted_datas <- lapply(beta_booted_list, shuffle_beta_rows)
  rm(beta_booted_list, booted_list)

  #### If n_cores > 1, check parallel library loading ####
  if(n_cores > 1){
    #Initiate the n_cores parallel clusters
    cl <- makePSOCKcluster(n_cores)
    setDefaultCluster(cl)

    #### Check for errors when directly loading a necessary library on each cluster ####
    test_load1 <- try(clusterEvalQ(cl, library(lmerTest)), silent = TRUE)

    #If errors produced, iterate through other options for library loading
    if(class(test_load1) == "try-error") {

      #Try finding unique library paths using system.file()
      pkgLibs <- unique(c(sub("/lmerTest$", "", system.file(package = "lmerTest"))))
      clusterExport(NULL, c('pkgLibs'), envir = environment())
      clusterEvalQ(cl, .libPaths(pkgLibs))

      #Check again for errors loading libraries
      test_load2 <- try(clusterEvalQ(cl, library(lmerTest)), silent = TRUE)

      if(class(test_load2) == "try-error"){

        #Try loading the user's .libPath() directly
        clusterEvalQ(cl,.libPaths(as.character(.libPaths())))
        test_load3 <- try(clusterEvalQ(cl, library(lmerTest)), silent = TRUE)

        if(class(test_load3) == "try-error"){

          #Give up and use lapply instead!
          parallel_compliant <- FALSE
          stopCluster(cl)

        } else {
          parallel_compliant <- TRUE
        }

      } else {
        parallel_compliant <- TRUE
      }

    } else {
      parallel_compliant <- TRUE
    }
  } else {
    #If n_cores = 1, set parallel_compliant to FALSE
    parallel_compliant <- FALSE
    warning('Parallel loading failed, calculations may crash!')
  }

  #### Estimate coefficient quantiles and rsquared ####
  if(parallel_compliant){

    #Export necessary data to each cluster
    clusterExport(NULL, c('n_bootstraps', 'booted_datas',
                          'beta_booted_datas',
                          'full_formula',
                          'beta_full_formula'),
                  envir = environment())

  all_cvs <- pbapply::pblapply(seq_len(n_bootstraps), function(x){

    ## Load lmerTest so that summary functions such as 'coef' and 'predict' function
    # properly
    require(lmerTest)

    ## Run the linear mixed model on the shuffled data
    train_mod <- suppressMessages(lmerTest::lmer(as.formula(full_formula),
                                                 data = booted_datas[[x]], REML = FALSE))

    ## Use step-wise selection (based on AIC scores) to choose the best predictors
    step_test <- lmerTest::step(train_mod, reduce.random = FALSE, steps = 200)
    rm(train_mod)

    best_mod <- lmerTest::get_model(step_test)
    rm(step_test)

    ## Calculate rsquared by comparing observed and predicted centralities
    booted_datas[[x]]$predict <- predict(best_mod)
    rsquared <- cor.test(booted_datas[[x]]$centrality, booted_datas[[x]]$predict)[[4]]

    ## Extract fixed and random effect coefficients
    sp_intercepts <- data.frame(coef(best_mod)$species)
    sp_intercepts <- data.frame(Species = rownames(sp_intercepts),
                                Intercept = sp_intercepts[,1])

    reg_intercepts <- data.frame(coef(best_mod)$region)
    reg_intercepts <- data.frame(Region = rownames(reg_intercepts),
                                 Intercept = reg_intercepts[,1])

    mod_coefs <- data.frame(Parameter = rownames(data.frame(summary(best_mod)$coefficients[-1, ])),
                            Coefficient = summary(best_mod)$coefficients[-1, 1])

    ## Run a linear mixed model using the community beta diversity output
    mod <- lme4::lmer(as.formula(beta_full_formula),
                      data = beta_booted_datas[[x]], REML = FALSE)

    beta_booted_datas[[x]]$predict <- predict(mod)
    beta_rsquared <- cor.test(beta_booted_datas[[x]]$beta_os, beta_booted_datas[[x]]$predict)[[4]]

    beta_reg_intercepts <- data.frame(coef(mod)$region)
    beta_reg_intercepts <- data.frame(Region = rownames(beta_reg_intercepts),
                                      Intercept = beta_reg_intercepts[,1])

    beta_mod_coefs <- data.frame(Parameter = rownames(data.frame(summary(mod)$coefficients[-1, ])),
                                 Coefficient = summary(mod)$coefficients[-1, 1])

    list(coefficients = mod_coefs,
         sp_intercepts = sp_intercepts,
         reg_intercepts = reg_intercepts,
         rsquared = rsquared,
         beta_coefficients = beta_mod_coefs,
         beta_reg_intercepts = beta_reg_intercepts,
         beta_rsquared = beta_rsquared)


  }, cl = cl)
  stopCluster(cl)

  } else {

    ## Use lapply if parallel loading fails
    all_cvs <- pbapply::pblapply(seq_len(n_bootstraps), function(x){

      ## Load lmerTest so that summary functions such as 'coef' and 'predict' function
      # properly
      require(lmerTest)

      ## Run the linear mixed model on the shuffled data
      train_mod <- suppressMessages(lmerTest::lmer(as.formula(full_formula),
                                                   data = booted_datas[[x]], REML = FALSE))

      ## Use step-wise selection (based on AIC scores) to choose the best predictors
      step_test <- lmerTest::step(train_mod, reduce.random = FALSE, steps = 200)
      rm(train_mod)

      best_mod <- lmerTest::get_model(step_test)
      rm(step_test)

      ## Calculate rsquared by comparing observed and predicted centralities
      booted_datas[[x]]$predict <- predict(best_mod)
      rsquared <- cor.test(booted_datas[[x]]$centrality, booted_datas[[x]]$predict)[[4]]

      ## Extract fixed and random effect coefficients
      sp_intercepts <- data.frame(coef(best_mod)$species)
      sp_intercepts <- data.frame(Species = rownames(sp_intercepts),
                                  Intercept = sp_intercepts[,1])

      reg_intercepts <- data.frame(coef(best_mod)$region)
      reg_intercepts <- data.frame(Region = rownames(reg_intercepts),
                                   Intercept = reg_intercepts[,1])

      mod_coefs <- data.frame(Parameter = rownames(data.frame(summary(best_mod)$coefficients[-1, ])),
                              Coefficient = summary(best_mod)$coefficients[-1, 1])

      ## Run a linear mixed model using the community beta diversity output
      mod <- lme4::lmer(as.formula(beta_full_formula),
                        data = beta_booted_datas[[x]], REML = FALSE)

      beta_booted_datas[[x]]$predict <- predict(mod)
      beta_rsquared <- cor.test(beta_booted_datas[[x]]$beta_os, beta_booted_datas[[x]]$predict)[[4]]

      beta_reg_intercepts <- data.frame(coef(mod)$region)
      beta_reg_intercepts <- data.frame(Region = rownames(beta_reg_intercepts),
                                        Intercept = beta_reg_intercepts[,1])

      beta_mod_coefs <- data.frame(Parameter = rownames(data.frame(summary(mod)$coefficients[-1, ])),
                                   Coefficient = summary(mod)$coefficients[-1, 1])

      list(coefficients = mod_coefs,
           sp_intercepts = sp_intercepts,
           reg_intercepts = reg_intercepts,
           rsquared = rsquared,
           beta_coefficients = beta_mod_coefs,
           beta_reg_intercepts = beta_reg_intercepts,
           beta_rsquared = beta_rsquared)

    })
  }

  #### Extract summary statistics of model coefficients ####
  coef.summary = do.call(rbind, purrr::map(all_cvs, 'coefficients')) %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(Lower.coef = round(quantile(Coefficient, probs = 0.025), 4),
                     Median.coef = round(quantile(Coefficient, probs = 0.5), 4),
                     Upper.coef = round(quantile(Coefficient, probs = 0.975), 4)) %>%
    dplyr::arrange(-abs(Median.coef))

  beta.coef.summary = do.call(rbind, purrr::map(all_cvs, 'beta_coefficients')) %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(Lower.coef = round(quantile(Coefficient, probs = 0.025), 4),
                     Median.coef = round(quantile(Coefficient, probs = 0.5), 4),
                     Upper.coef = round(quantile(Coefficient, probs = 0.975), 4)) %>%
    dplyr::arrange(-abs(Median.coef))

  sp.intercept.summary = do.call(rbind, purrr::map(all_cvs, 'sp_intercepts')) %>%
    dplyr::group_by(Species) %>%
    dplyr::mutate(Intercept = boot::inv.logit(Intercept)) %>%
    dplyr::summarise(Lower.intercept = round(quantile(Intercept, probs = 0.025), 4),
                     Median.intercept = round(quantile(Intercept, probs = 0.5), 4),
                     Upper.intercept = round(quantile(Intercept, probs = 0.975), 4)) %>%
    dplyr::arrange(-abs(Median.intercept))

  reg.intercept.summary = do.call(rbind, purrr::map(all_cvs, 'reg_intercepts')) %>%
    dplyr::group_by(Region) %>%
    dplyr::mutate(Intercept = boot::inv.logit(Intercept)) %>%
    dplyr::summarise(Lower.intercept = round(quantile(Intercept, probs = 0.025), 4),
                     Median.intercept = round(quantile(Intercept, probs = 0.5), 4),
                     Upper.intercept = round(quantile(Intercept, probs = 0.975), 4)) %>%
    dplyr::arrange(-abs(Median.intercept))

  beta.reg.intercept.summary = do.call(rbind, purrr::map(all_cvs, 'beta_reg_intercepts')) %>%
    dplyr::group_by(Region) %>%
    dplyr::summarise(Lower.intercept = round(quantile(Intercept, probs = 0.025), 4),
                     Median.intercept = round(quantile(Intercept, probs = 0.5), 4),
                     Upper.intercept = round(quantile(Intercept, probs = 0.975), 4)) %>%
    dplyr::arrange(-abs(Median.intercept))

  rsquared.summary = data.frame(r = do.call(rbind, purrr::map(all_cvs, 'rsquared')),
                                mod = 1) %>%
    dplyr::group_by(mod) %>%
    dplyr::summarise(Lower.rsquared = round(quantile(cor, probs = 0.025), 4),
                     Median.rsquared = round(quantile(cor, probs = 0.5), 4),
                     Upper.rsquared = round(quantile(cor, probs = 0.975), 4)) %>%
    dplyr::select(-mod)

  beta.rsquared.summary = data.frame(r = do.call(rbind, purrr::map(all_cvs, 'beta_rsquared')),
                                mod = 1) %>%
    dplyr::group_by(mod) %>%
    dplyr::summarise(Lower.rsquared = round(quantile(cor, probs = 0.025), 4),
                     Median.rsquared = round(quantile(cor, probs = 0.5), 4),
                     Upper.rsquared = round(quantile(cor, probs = 0.975), 4)) %>%
    dplyr::select(-mod)

  #### Return the long-format datasets and the model summaries ####
  return(list(sp_cent_mod_data = cent_mod_data,
              betadiv_mod_data = betadiv_mod_data,
              sp_cent_coefs = coef.summary,
              sp_cent_intercepts = sp.intercept.summary,
              sp_cent_rsquared = rsquared.summary,
              betadiv_coefs = beta.coef.summary,
              betadiv_intercepts = beta.reg.intercept.summary,
              betadiv_rsquared = beta.rsquared.summary))

}

