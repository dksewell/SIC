#' Perform inference for the SIC model 
#' 
#' Obtain posterior samples based on data fit to
#' the SIC model.  This will return output for the incidence and 
#' the clearance models separately.
#' 
#' @param formula Formula of the form <y> ~ <x1> + <x2> + ... + (<time variable> | <subject ID>)
#' @param data Data frame or tibble
#' @param seed a single value, interpreted as an integer.  Used for replication purposes.
#' @param min_count_to_estimate positive integer giving the minimum cell count in a 2x2 
#' contingency table required to regress one pathogen on another
#' @param CI_level numeric between 0 and 1 for the credible interval level.
#' @param verbose logical.  Whether to print messages as to where the model fitting algorithm is.
#' @param ... Further arguments to pass to stan_glm().
#' @returns  Object of class 'sic' which has the following elements:
#' * results tibble with columns telling the response variable, the model (incidence or 
#' clearance), the covariate, the regression coefficient estimate, the lower and upper CI 
#' bounds, and the probability of direction (i.e., max(Prob(beta > 0 | data),Prob(beta < 0 | data))).
#' *  bic the BIC of the fitted model (lower is better)
#' * estimability list with the valid response variables and the valid covariates 
#' for each response variable
#' * networks list of the PxP weighted adjacency matrix giving the regression coefficients 
#' as the edge weights.  E.g., the (i,j)-th element of incidence matrix is the effect of 
#' the presence of pathogen i on the log incidence rate of pathogen j.
#' * stanreg fits The individual fits of the incidence and clearance models
#' * CI_level
#' * min_count_to_estimate
#' 
#' @examples 
#' sic_data = sic_simulator(seed = 2023)
#' test = sic(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject),
#'                  sic_data$data[[1]] %>% select(-tau))
#' 
#' @export
#' 


if(FALSE){
  library(magrittr)
  library(rstanarm)
  library(dplyr)
  library(plot.matrix)
  library(viridis)
  import::from(ggplot2,ggtitle)
  library(adaptMCMC)
  source("~/SIC/R/sic_simulator.R")
  source("~/SIC/R/sic_frequentist.R")
  source("~/SIC/R/sic.R")
  source("~/SIC/R/plot.sic.R")
  
  sic_data = 
    sic_simulator(seed = 2023)
  
  # test = sic_lowrank(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject),
  #                    sic_data$data[[1]] %>% select(-tau),
  #                    verbose = TRUE,
  #                    seed = 2023)
  
  formula = cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject)
  data = sic_data$data[[1]] %>% select(-tau)
  seed = 2023
  
  min_count_to_estimate = 5
  CI_level = 0.95
  verbose = TRUE
  lrank = 2
  prior_beta_I = normal(location = 0, scale = 2.5, autoscale = TRUE)
  prior_beta_C = normal(location = 0, scale = 2.5, autoscale = TRUE)
  prior_U = normal(location = 0, scale = sqrt(2.5) / lrank^(1/4))
  prior_V_I = normal(location = 0, scale = sqrt(2.5) / lrank^(1/4))
  prior_V_C = normal(location = 0, scale = sqrt(2.5) / lrank^(1/4))
  n_draws = 5e3
}

sic_lowrank = function(formula,
                       data, 
                       seed = NULL,
                       min_count_to_estimate = 5,
                       CI_level = 0.95,
                       lrank = 2,
                       n_draws = 5e3,
                       prior_I = normal(location = 0, scale = 2.5, autoscale = TRUE),
                       prior_C = normal(location = 0, scale = 2.5, autoscale = TRUE),
                       prior_U = normal(location = 0, scale = sqrt(2.5) / lrank^(1/4)),
                       prior_V_I = normal(location = 0, scale = sqrt(2.5) / lrank^(1/4)),
                       prior_V_C = normal(location = 0, scale = sqrt(2.5) / lrank^(1/4)),
                       verbose = TRUE){
  
  if(!is.null(seed)){
    set.seed(seed)
  }else{
    warning("Make sure to set your seed!")
  }
  
  if(verbose) cat("\nCurating data and determining which responses/covariates are estimable\n")
  # Drop NAs
  data %<>% 
    select(all.vars(formula)) %>% 
    na.omit()
  # tidyr::drop_na(all.vars(formula))
  
  # get names of ...
  form_char = as.character(formula)[-1]
  form_char[1] = gsub("cbind","",form_char[1])
  form_char[1] = gsub("[[:punct:]]","",form_char[1])
  ## pathogen columns
  pathogen_vars = 
    strsplit(form_char[1]," ")[[1]]
  ## covariates
  X_vars = 
    gsub("\ \\+","",
         substr(form_char[2],
                start = 1,
                stop = gregexpr("\\(",form_char[2])[[1]][1] - 3)
    ) %>% 
    strsplit(" ")
  X_vars = X_vars[[1]]
  ## subject variable
  subject_var = 
    gsub('[[:punct:]]',"",
         substr(form_char[2],
                start = gregexpr("\\|",form_char[2])[[1]][1] + 2,
                stop = 1e4)
    )
  ## time variable
  time_var = 
    gsub('[[:punct:]]',"",
         substr(form_char[2],
                start = gregexpr("\\(",form_char[2])[[1]][1] + 1,
                stop = gregexpr("\\|",form_char[2])[[1]][1] - 2)
    )
  
  
  # Store basic quantities
  NT = nrow(data)
  N = 
    data %>% 
    select(all_of(subject_var)) %>% 
    unlist() %>% 
    unique() %>% 
    length()
  P = length(pathogen_vars)
  Q = length(X_vars) + 1 # + 1 for the intercept
  
  # arrange data by subject and by time
  data %<>%
    arrange(across(starts_with(subject_var)),
            across(starts_with(time_var)))
  
  # Compute time between observations (tau)
  data %<>%
    mutate(tau = 
             ifelse(
               c(TRUE,data[[subject_var]][-1] == 
                   data[[subject_var]][-NT]),
               diff(c(data[[time_var]][1],data[[time_var]])),
               0.0)
    )
  
  # Create lagged variables
  for(p in 1:P){
    data[[paste0(pathogen_vars[p],"_lagged")]] = 
      ifelse(
        c(FALSE,
          data[[subject_var]][-1] == data[[subject_var]][-NT]),
        c(-999,data[[pathogen_vars[p]]][-NT]),
        NA)
  }
  data %<>%
    na.omit()
  
  
  # Distinguish between incidence and clearance 
  incidence_index = 
    clearance_index = 
    list()
  for(p in 1:P){
    incidence_index[[p]] =
      which(data[[paste0(pathogen_vars[p],"_lagged")]] == 0)
    clearance_index[[p]] =
      which(data[[paste0(pathogen_vars[p],"_lagged")]] == 1)
  }
  
  # all_incidence = 
  #   cbind( unlist(incidence_index),
  #          rep(1:P,sapply(incidence_index,length)) 
  #   )
  # all_clearance = 
  #   cbind( unlist(clearance_index),
  #          rep(1:P,sapply(clearance_index,length)) 
  #   )
  
  
  
  # Compute tables to see if we can make inference
  valid_responses = list()
  
  for(p in 1:P){
    
    if( (length(unique(data[[pathogen_vars[p]]][ incidence_index[[p]] ] )) > 1) &
        (min(table(data[[pathogen_vars[p]]][ incidence_index[[p]]])) > min_count_to_estimate ) ){
      valid_responses$incidence =
        c(valid_responses$incidence,p)
    }
    
    if( (length(unique(data[[pathogen_vars[p]]][ clearance_index[[p]] ] )) > 1) &
        (min(table(data[[pathogen_vars[p]]][ clearance_index[[p]]])) > min_count_to_estimate ) ){
      valid_responses$clearance =
        c(valid_responses$clearance,p)
      
    }
    
  }
  
  
  
  # Reverse code clearance pathogen variables
  for(p in 1:P){
    data[[pathogen_vars[p]]][clearance_index[[p]]] =
      1 - data[[pathogen_vars[p]]][clearance_index[[p]]]
  }
  
  # Get X explicitly
  X_matrix = 
    model.matrix(as.formula(paste0("~ ",
                                   paste(X_vars,collapse = "+"))),
                 data = data)
  
  # Get y_(t-1) explicitly
  y_lagged = 
    data %>% 
    select(all_of(paste(pathogen_vars,"lagged",sep="_"))) %>% 
    as.matrix()
  
  # Get y^* (vectorized) and tau
  y_star = 
    data %>% 
    select(all_of(pathogen_vars)) %>% 
    unlist()
  # y_star_mat = 
  #   data %>% 
  #   select(all_of(pathogen_vars)) %>% 
  #   as.matrix()
  tau_rep = rep(data$tau,P)
  
  if(verbose) cat("\nInitializing parameters\n")
  
  init_fit = NULL
  warn = options()$warn
  options(warn = -1)
  try({
    init_fit = 
      suppressMessages(sic_frequentist(formula,data))
  },silent = TRUE)
  if(is.null(init_fit)){
    try({
      init_fit = 
        suppressMessages(sic_frequentist(formula,data, min_count_to_estimate = 10))
    },silent = TRUE)
  }
  if(is.null(init_fit)){
    init_fit = 
      suppressMessages(sic_frequentist(formula,data, min_count_to_estimate = 2 * P))
  }
  options(warn = warn); rm(warn)
  
  init_values = 
    list(Beta_I = matrix(0.0,Q,P),
         Beta_C = matrix(0.0,Q,P),
         U = matrix(0.0,P,lrank),
         V_I = matrix(0.0,P,lrank),
         V_C = matrix(0.0,P,lrank))
  
  temp = tibble(Covariate = c("(Intercept)",X_vars))
  for(j in 1:P){
    temp =
      temp %>% 
      left_join(init_fit$results %>% 
                  filter(Response == pathogen_vars[j],
                         Model == "Incidence") %>% 
                  select(Covariate,Estimate),
                by = "Covariate")
  }
  init_values$Beta_I = matrix(unlist(temp[,-1]),Q,P,
                              dimnames = list(Covariate = c("Intercept",X_vars),
                                              Response = pathogen_vars))
  temp = tibble(Covariate = c("(Intercept)",X_vars))
  for(j in 1:P){
    temp =
      temp %>% 
      left_join(init_fit$results %>% 
                  filter(Response == pathogen_vars[j],
                         Model == "Clearance") %>% 
                  select(Covariate,Estimate),
                by = "Covariate")
  }
  init_values$Beta_C = matrix(unlist(temp[,-1]),Q,P,
                              dimnames = list(Covariate = c("Intercept",X_vars),
                                              Response = pathogen_vars))
  
  temp = tibble(Covariate = paste(pathogen_vars,"lagged",sep="_"))
  for(j in 1:P){
    temp =
      temp %>% 
      left_join(init_fit$results %>% 
                  filter(Response == pathogen_vars[j],
                         Model == "Incidence") %>% 
                  select(Covariate,Estimate),
                by = "Covariate")
  }
  eta_I = 
    matrix(unlist(temp[,-1]),P,P,
           dimnames = list(Covariate = pathogen_vars,
                           Response = pathogen_vars))
  eta_I[which(is.na(eta_I))] = 0.0
  temp = tibble(Covariate = paste(pathogen_vars,"lagged",sep="_"))
  for(j in 1:P){
    temp =
      temp %>% 
      left_join(init_fit$results %>% 
                  filter(Response == pathogen_vars[j],
                         Model == "Clearance") %>% 
                  select(Covariate,Estimate),
                by = "Covariate")
  }
  eta_C = 
    matrix(unlist(temp[,-1]),P,P,
           dimnames = list(Covariate = pathogen_vars,
                           Response = pathogen_vars))
  eta_C[which(is.na(eta_C))] = 0.0
  eta_svd = svd(cbind(eta_I,eta_C))
  init_values$U = 
    eta_svd$u[,1:lrank] %*% diag(sqrt(eta_svd$d[1:lrank]))
  init_values$V_I = 
    eta_svd$v[1:P,1:lrank] %*% diag(sqrt(eta_svd$d[1:lrank]))
  init_values$V_C = 
    eta_svd$v[P + 1:P,1:lrank] %*% diag(sqrt(eta_svd$d[1:lrank]))
  rm(temp,eta_I,eta_C,eta_svd,form_char)
  
  
  # Set prior hyperparameter values
  prior_list = list(beta_I = list(mean = prior_beta_I$location),
                    beta_C = list(mean = prior_beta_I$location),
                    U = list(mean = prior_U$location,
                             sd = prior_U$scale),
                    V_I = list(mean = prior_V_I$location,
                               sd = prior_V_I$scale),
                    V_C = list(mean = prior_V_C$location,
                               sd = prior_V_C$scale))
  sd_x = apply(X_matrix[,-1],2,sd)
  if(prior_beta_I$autoscale){
    prior_list$beta_I$sd = 
      prior_beta_I$scale * matrix(abs(init_values$Beta_I[1,]),nr=1)
    for(j in 2:Q) prior_list$beta_I$sd = rbind(prior_list$beta_I$sd,
                                               prior_beta_I$scale * sd_x[j - 1])
  }else{
    prior_list$beta_I$sd = prior_beta_I$scale
  }
  if(prior_beta_C$autoscale){
    prior_list$beta_C$sd = 
      prior_beta_C$scale * matrix(abs(init_values$Beta_C[1,]),nr=1)
    for(j in 2:Q) prior_list$beta_C$sd = rbind(prior_list$beta_C$sd,
                                               prior_beta_C$scale * sd_x[j - 1])
  }else{
    prior_list$beta_C$sd = prior_beta_C$scale
  }
  
  
  # Prior function
  prior_logpdf = function(beta_I,beta_C,U,V_I,V_C){
    sum(dnorm(c(beta_I),
              mean = c(prior_list$beta_I$mean),
              sd = c(prior_list$beta_I$sd),
              log = TRUE)) + 
      sum(dnorm(c(beta_C),
                mean = c(prior_list$beta_C$mean),
                sd = c(prior_list$beta_C$sd),
                log = TRUE)) +
      sum(dnorm(c(U),
                mean = c(prior_list$U$mean),
                sd = c(prior_list$U$sd),
                log = TRUE)) +
      sum(dnorm(c(V_I),
                mean = c(prior_list$V_I$mean),
                sd = c(prior_list$V_I$sd),
                log = TRUE)) +
      sum(dnorm(c(V_C),
                mean = c(prior_list$V_C$mean),
                sd = c(prior_list$V_C$sd),
                log = TRUE))
  }
  
  # log posterior function
  
  lp = function(x){
    beta_I = matrix(x[1:(P*Q)], Q, P)
    beta_C = matrix(x[P*Q + 1:(P*Q)], Q, P)
    U = matrix(x[2*P*Q + 1:(P*lrank)], P, lrank)
    V_I = matrix(x[2*P*Q + P*lrank + 1:(P*lrank)], P, lrank)
    V_C = matrix(x[2*P*Q + 2*P*lrank + 1:(P*lrank)], P, lrank)
    
    eta_I = tcrossprod(U,V_I)
    eta_C = tcrossprod(U,V_C)
    
    Lambda = 
      exp(
        (1.0 - y_lagged) * ( X_matrix %*% beta_I + y_lagged %*% eta_I )  +
          (y_lagged) * ( X_matrix %*% beta_C + y_lagged %*% eta_C )
      )
    
    y_probs = 1 - exp(-c(Lambda) * tau_rep)
    
    
    return(
      sum(dbinom(y_star,1,y_probs,log = TRUE)) +
        prior_logpdf(beta_I,beta_C,U,V_I,V_C)
      )
  }
  
  # ----------------------
  
  if(verbose) cat("\nUsing Adaptive MCMC to fit SIC model\n")
  
  mcmc_draws = 
    MCMC(p = lp,
         n = n_draws,
         init = 
           c(c(init_values$Beta_I),
             c(init_values$Beta_C),
             c(init_values$U),
             c(init_values$V_I),
             c(init_values$V_C)
           ),
         acc.rate = 0.234,
         showProgressBar = verbose)
  
  
  heidel = 
    heidel.diag(mcmc(mcmc_draws$samples[,1:(2*P*Q)]))
  if(any(heidel[,"stest"] == 0)){
    # Second try
    cat("\nChain not converged- producing additional samples")
    mcmc_draws = 
      MCMC.add.samples(mcmc_draws,
                       n_draws)
    heidel = 
      heidel.diag(mcmc(mcmc_draws$samples[,1:(2*P*Q)]))
    
    if(any(heidel[,"stest"] == 0)){
      # Third try
      cat("\nChain still not converged- producing additional samples")
      
      mcmc_draws = 
        MCMC.add.samples(mcmc_draws,
                         n_draws)
      heidel = 
        heidel.diag(mcmc(mcmc_draws$samples[,1:(2*P*Q)]))
      
      if(any(heidel[,"stest"] == 0)){
        # Fourth and final try
        cat("\nChain still not converged- last attempt.")
        
        mcmc_draws = 
          MCMC.add.samples(mcmc_draws,
                           2 * n_draws)
        heidel = 
          heidel.diag(mcmc(mcmc_draws$samples[-c(1:(3 * n_draws)),1:(2*P*Q)]))
        
        
      }
      
    }
    
  }
  
  
  
  geweke =
    geweke.diag(mcmc(mcmc_draws$samples[-c(1:burnin),1:(2*P*Q)]))
  qqnorm(geweke$z)
  qqline(geweke$z)
  shapiro.test(geweke$z)
  
  beta_I = beta_R = 
    array(0.0,c(n_draws,Q,P))
  eta_I = eta_C = 
    array(0.0, c(n_draws,P,P))
  
  #################
  ### Auto-Tune ###
  #################
  accRate = 0
  autoTune_count = 0
  stepSize = 0.01/(2 * (Q*P + P*lrank) + P*lrank)^(1/6) #Roberts and Rosenthal (1998), setting ell = 0.01 (Journal of the Royal Statistical Society. Series B (Statistical Methodology), Vol. 60, No. 1. (1998), pp. 255-268.)
  stepSize = c(stepSize,NA,NA)
  turnNumber = 1
  gr = (1+sqrt(5))*0.5
  new_draws = 
    list(q = c(c(init_values$Beta_I),
               c(init_values$Beta_C),
               c(init_values$U),
               c(init_values$V_I),
               c(init_values$V_C)
               ),
         accept = 0)
  
  test = nlp_grad(new_draws$q)
  test1 = numDeriv::grad(nlp,new_draws$q)
  plot(test~test1,
       col = c(1:2)[1 + ((abs(test - test1) > 1) & test < 25)])
  test_index = which((abs(test - test1) > 1) & test < 25)
  
  
  while( ((accRate < autoTune_bounds[1]) | (accRate > autoTune_bounds[2])) & (autoTune_count < max_autoTune_count)){
    
    ### First turn of auto-tuning
    if(turnNumber == 1){
      accRate = 0
      for(it in 1:autoTune_n){
        new_draws = HMC(nlp,nlp_grad,stepSize[1],n_lf_steps,new_draws$q)
        accRate = accRate + new_draws$accept/autoTune_n
      }
      cat(paste0("--- Acceptance Rate = ",round(accRate,4)," ---\n"))
      autoTune_count = 1
    }
    
  }
  
  
  
  
  # Collect results
  estimates = 
    list(incidence = lapply(fits$incidence[which(!sapply(fits$incidence,is.null))],coef),
         clearance = lapply(fits$clearance[which(!sapply(fits$clearance,is.null))],coef)
    )
  
  CIs = 
    list(incidence = lapply(fits$incidence[which(!sapply(fits$incidence,is.null))],posterior_interval,prob = CI_level),
         clearance = lapply(fits$clearance[which(!sapply(fits$clearance,is.null))],posterior_interval,prob = CI_level)
    )
  
  helper_fun = function(x){
    x_matrix = as.matrix(x)
    results = apply(x_matrix,2,function(y) mean(y > 0))
    results = sapply(results, function(y) max(y,1-y))
    return(results)
  }
  
  prob_dir = 
    list(incidence = lapply(fits$incidence[which(!sapply(fits$incidence,is.null))],helper_fun),
         clearance = lapply(fits$clearance[which(!sapply(fits$clearance,is.null))],helper_fun)
    )
  
  
  
  names(estimates$incidence) = 
    names(CIs$incidence) = 
    names(prob_dir$incidence) = 
    pathogen_vars[valid_responses$incidence]
  names(estimates$clearance) = 
    names(CIs$clearance) = 
    names(prob_dir$clearance) = 
    pathogen_vars[valid_responses$clearance]
  
  
  results = 
    tibble(Response = character(),
           Model = character(),
           Covariate = character(),
           Estimate = numeric(),
           Lower = numeric(),
           Upper = numeric(),
           `Prob of Direction` = numeric())
  
  for(i in 1:P){
    if(i %in% valid_responses$incidence){
      results %<>%
        bind_rows(
          tibble(Response = pathogen_vars[i],
                 Model = "Incidence",
                 Covariate = rownames(CIs$incidence[[pathogen_vars[i]]]),
                 Estimate = estimates$incidence[[pathogen_vars[i]]],
                 Lower = CIs$incidence[[pathogen_vars[i]]][,1],
                 Upper = CIs$incidence[[pathogen_vars[i]]][,2],
                 `Prob of Direction` = prob_dir$incidence[[pathogen_vars[i]]])
        )
    }else{
      results %<>%
        bind_rows(
          tibble(Response = pathogen_vars[i],
                 Model = "Incidence",
                 Covariate = NA,
                 Estimate = NA,
                 Lower = NA,
                 Upper = NA,
                 `Prob of Direction` = NA)
        )
    }
    
    if(i %in% valid_responses$clearance){
      results %<>%
        bind_rows(
          tibble(Response = pathogen_vars[i],
                 Model = "Clearance",
                 Covariate = rownames(CIs$clearance[[pathogen_vars[i]]]),
                 Estimate = estimates$clearance[[pathogen_vars[i]]],
                 Lower = CIs$clearance[[pathogen_vars[i]]][,1],
                 Upper = CIs$clearance[[pathogen_vars[i]]][,2],
                 `Prob of Direction` = prob_dir$clearance[[pathogen_vars[i]]])
        )
    }else{
      results %<>%
        bind_rows(
          tibble(Response = pathogen_vars[i],
                 Model = "Clearance",
                 Covariate = NA,
                 Estimate = NA,
                 Lower = NA,
                 Upper = NA,
                 `Prob of Direction` = NA)
        )
    }
  }
  
  # Create pathogen-pathogen network
  incidence_network = 
    clearance_network = 
    matrix(NA, P, P,
           dimnames = list(source = pathogen_vars,
                           target = pathogen_vars))
  
  temp1 = 
    results %>% 
    filter(Model == "Incidence") %>% 
    select("Response","Covariate","Estimate") %>% 
    mutate(Covariate = gsub("_lagged","",Covariate)) %>% 
    filter(Covariate %in% pathogen_vars)
  incidence_network[cbind(temp1$Covariate,
                          temp1$Response)] = temp1$Estimate
  temp1 = 
    results %>% 
    filter(Model == "Clearance") %>% 
    select("Response","Covariate","Estimate") %>% 
    mutate(Covariate = gsub("_lagged","",Covariate)) %>% 
    filter(Covariate %in% pathogen_vars)
  clearance_network[cbind(temp1$Covariate,
                          temp1$Response)] = temp1$Estimate
  
  if(verbose) cat("\nComputing the BIC\n")
  # Compute BIC
  bic_helper = function(x){
    ll = max(rowSums(log_lik(x)))
    num_parms = length(coef(x))
    
    -2 * ll + num_parms * log(nrow(x$data))
  }
  bic_values = 
    list(incidence = lapply(fits$incidence[which(!sapply(fits$incidence,is.null))],bic_helper),
         clearance = lapply(fits$clearance[which(!sapply(fits$clearance,is.null))],bic_helper)
    )
  
  bic = sum(unlist(bic_values))
  
  
  valid_responses$incidence = 
    pathogen_vars[valid_responses$incidence]
  valid_responses$clearance = 
    pathogen_vars[valid_responses$clearance]
  for(i in 1:P){
    if(!is.null(valid_covariates$incidence[[i]])) valid_covariates$incidence[[i]] = pathogen_vars[valid_covariates$incidence[[i]]]
    if(!is.null(valid_covariates$clearance[[i]])) valid_covariates$clearance[[i]] = pathogen_vars[valid_covariates$clearance[[i]]]
  }
  names(valid_covariates$incidence) = pathogen_vars
  names(valid_covariates$clearance) = pathogen_vars
  
  object = 
    list(results = results,
         bic = bic,
         estimability = list(response = valid_responses,
                             covariates = valid_covariates),
         networks = list(incidence = incidence_network,
                         clearance = clearance_network),
         stanreg_fits = fits,
         CI_level = CI_level,
         min_count_to_estimate = min_count_to_estimate)
  
  class(object) = "sic"
  
  return(object)
}














