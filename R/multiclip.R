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
  library(dplyr)
  library(plot.matrix)
  library(viridis)
  import::from(ggplot2,ggtitle)
  # library(mvtnorm)
  source("~/SIC/R/sic_simulator.R")
  # source("~/SIC/R/sic_frequentist.R")
  # source("~/SIC/R/sic.R")
  # source("~/SIC/R/plot.sic.R")

  sic_data =
    sic_simulator(seed = 2023)

  # test = sic_lowrank(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject),
  #                    sic_data$data[[1]] %>% select(-tau),
  #                    verbose = TRUE,
  #                    seed = 2023)

  formula_clearance =  cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 + (time | subject)
  formula_incidence =  ~ x1 + x2 + (time | subject)
  formula_prevalence =  ~ x1 + (time | subject)
  data = sic_data$data[[1]] %>% select(-tau)
  seed = 2023

  min_count_to_estimate = 5
  CI_level = 0.95
  verbose = TRUE
  prior_beta_I = list(location = 0, scale = 2.5, autoscale = TRUE)
  prior_beta_C = list(location = 0, scale = 2.5, autoscale = TRUE)
  prior_eta_I = list(location = 0, scale = 2.5, autoscale = FALSE)
  prior_eta_C = list(location = 0, scale = 2.5, autoscale = FALSE)
  prior_prev = list(location = 0, scale = 2.5, autoscale = TRUE)
  prior_sensitivity = list(a = 18.5, b = 3.9)
  prior_specificity = list(a = 50, b = 1)
  n_draws = 5e4
  MH_scalars = list(clearance = 1/5, incidence = 1/5, prevalence = 1/5, z_itp = 0.5)


  temp_fun = function(x){
    (qbeta(0.025,x[1], x[2]) - 0.65)^2 +
       (qbeta(0.975,x[1], x[2]) - 0.95)^2
  }
  optim(c(200,10),temp_fun)
  curve(dbeta(x,18.5, 3.9))
  qbeta(c(0.025,0.975),18.5, 3.9)
  qbeta(c(0.025,0.975),50, 1)
  rm(temp_fun)
}

multiclip = function(formula_clearance =  cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 + (time | subject),
                     formula_incidence =  ~ x1 + x2 + (time | subject),
                     formula_prevalence =  ~ x1 + (time | subject),
                     data,
                     seed = NULL,
                     min_count_to_estimate = 5,
                     CI_level = 0.95,
                     n_draws = 5e4,
                     prior_beta_I = list(location = 0, scale = 2.5, autoscale = TRUE),
                     prior_beta_C = list(location = 0, scale = 2.5, autoscale = TRUE),
                     prior_eta_I = list(location = 0, scale = 2.5, autoscale = FALSE),
                     prior_eta_C = list(location = 0, scale = 2.5, autoscale = FALSE),
                     prior_prev = list(location = 0, scale = 2.5, autoscale = TRUE),
                     prior_sensitivity = list(a = 18.5, b = 3.9), # Puts 95% prior probability between 0.65 and 0.95
                     prior_specificity = list(a = 50, b = 1), # Puts 95% prior probability between 0.929 and 0.999
                     MH_scalars = list(clearance = 1/5, incidence = 1/5, prevalence = 1/5, z_itp = 0.5),
                     verbose = TRUE){

  if(!is.null(seed)){
    set.seed(seed)
  }else{
    warning("Make sure to set your seed!")
  }


  # Curating data -----------------------------------------------------------

  if(verbose) cat("\nCurating data and determining which responses/covariates are estimable\n")
  # Drop NAs
  data_aug =
    data %>%
    select(c(all.vars(formula_prevalence),
             all.vars(formula_incidence),
             all.vars(formula_clearance))) %>%
    na.omit()

 # get names of ...
  form_char = as.character(formula_clearance)[-1]
  form_char[1] = gsub("cbind","",form_char[1])
  form_char[1] = gsub("[[:punct:]]","",form_char[1])
  ## pathogen columns
  pathogen_vars =
    strsplit(form_char[1]," ")[[1]]
  ## covariates
  X_vars = list()
  X_vars$clearance =
    gsub("\ \\+","",
         substr(form_char[2],
                start = 1,
                stop = gregexpr("\\(",form_char[2])[[1]][1] - 3)
    ) %>%
    strsplit(" ")
  X_vars$clearance = X_vars$clearance[[1]]

  form_char = as.character(formula_incidence)[-1]
  X_vars$incidence =
    gsub("\ \\+","",
         substr(form_char,
                start = 1,
                stop = gregexpr("\\(",form_char)[[1]][1] - 3)
    ) %>%
    strsplit(" ")
  X_vars$incidence = X_vars$incidence[[1]]

  form_char = as.character(formula_prevalence)[-1]
  X_vars$prevalence =
    gsub("\ \\+","",
         substr(form_char,
                start = 1,
                stop = gregexpr("\\(",form_char)[[1]][1] - 3)
    ) %>%
    strsplit(" ")
  X_vars$prevalence = X_vars$prevalence[[1]]


  ## subject variable
  subject_var =
    gsub('[[:punct:]]',"",
         substr(form_char,
                start = gregexpr("\\|",form_char)[[1]][1] + 2,
                stop = 1e4)
    )
  ## time variable
  time_var =
    gsub('[[:punct:]]',"",
         substr(form_char,
                start = gregexpr("\\(",form_char)[[1]][1] + 1,
                stop = gregexpr("\\|",form_char)[[1]][1] - 2)
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
  Q = lapply(X_vars, function(x) length(x) + 1) # + 1 for the intercept

  # arrange data by subject and by time
  data_aug %<>%
    arrange(across(starts_with(subject_var)),
            across(starts_with(time_var)))

  # Compute time between observations (tau)
  data_aug %<>%
    mutate(time0 = c(TRUE,data_aug[[subject_var]][-1] !=
                       data_aug[[subject_var]][-NT])) %>%
    mutate(tau =
             ifelse(!time0,
                    diff(c(data_aug[[time_var]][1],data_aug[[time_var]])),
                    0.0)
    )

  # Create lagged variables
  for(p in 1:P){
    data_aug[[paste0(pathogen_vars[p],"_lagged")]] =
      ifelse(!data_aug$time0,
             c(-999,data_aug[[pathogen_vars[p]]][-NT]),
             NA)
  }

  # Create place to store estimated true values
  for(p in 1:P){
    data_aug[[paste0(pathogen_vars[p],"_est")]] =
      data_aug[[pathogen_vars[p]]]
  }

  # data_aug_old = data_aug

  # Filter data by time = 0 or time >= 1
  data_gr0 =
    data_aug %>%
    na.omit()
  data_0 =
    data_aug %>%
    filter(time0)



  # Determine estimability --------------------------------------------------

  incidence_index =
    clearance_index =
    list()
  for(p in 1:P){
    incidence_index[[p]] =
      which(data_gr0[[paste0(pathogen_vars[p],"_lagged")]] == 0)
    clearance_index[[p]] =
      which(data_gr0[[paste0(pathogen_vars[p],"_lagged")]] == 1)
  }

  valid_responses = list()
  valid_covariates = list()
  valid_covariates$incidence =
    valid_covariates$clearance = list()

  tables = list()
  tables$incidence =
    tables$clearance =
    list()

  for(p in 1:P){
    tables$incidence[[pathogen_vars[p]]] =
      tables$clearance[[pathogen_vars[p]]] =
      list()

    if( (length(unique(data_gr0[[pathogen_vars[p]]][ incidence_index[[p]] ] )) > 1) &
        (min(table(data_gr0[[pathogen_vars[p]]][ incidence_index[[p]]])) > min_count_to_estimate ) ){
      valid_responses$incidence =
        c(valid_responses$incidence,p)

      no_valid_covars = TRUE
      for(p2 in c(1:P)[-p]){
        tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]] =
          table(data_gr0[[pathogen_vars[p]]][ incidence_index[[p]] ] ,
                data_gr0[[paste0(pathogen_vars[p2],"_lagged")]][ incidence_index[[p]] ] )

        names(dimnames(tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]] )) =
          pathogen_vars[c(p,p2)]

        if( (ncol(tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]]) > 1) &
            (min(tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]]) >= min_count_to_estimate) ){
          if(no_valid_covars){
            valid_covariates$incidence[[ p ]] = p2
            no_valid_covars = FALSE
          }else{
            valid_covariates$incidence[[ p ]] =
              c(valid_covariates$incidence[[ p ]],
                p2)
          }
        }

      }

    }

    if( (length(unique(data_gr0[[pathogen_vars[p]]][ clearance_index[[p]] ] )) > 1) &
        (min(table(data_gr0[[pathogen_vars[p]]][ clearance_index[[p]]])) >= min_count_to_estimate ) ){
      valid_responses$clearance =
        c(valid_responses$clearance,p)

      no_valid_covars = TRUE
      for(p2 in c(1:P)[-p]){
        tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]] =
          table(data_gr0[[pathogen_vars[p]]][ clearance_index[[p]] ] ,
                data_gr0[[paste0(pathogen_vars[p2],"_lagged")]][ clearance_index[[p]] ] )

        names(dimnames(tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]] )) =
          pathogen_vars[c(p,p2)]

        if( (ncol(tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]]) > 1) &
            (min(tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]]) >= min_count_to_estimate) ){
          if(no_valid_covars){
            valid_covariates$clearance[[ p ]] = p2
            no_valid_covars = FALSE
          }else{
            valid_covariates$clearance[[ p ]] =
              c(valid_covariates$clearance[[ p ]],
                p2)
          }
        }

      }

    }
  }

  tables$prevalence =
    data_0 %>%
    dplyr::select(all_of(pathogen_vars)) %>%
    as.matrix() %>%
    colSums()

  valid_responses$prevalence =
    which(tables$prevalence >= min_count_to_estimate)

  valid_responses$combined =
    intersect(valid_responses$incidence,
              valid_responses$clearance)


  # Get more basic quantities -----------------------------------------------

  # Get tau as replicated columns in a matrix
  tau_matrix = matrix(data_gr0$tau,nrow(data_gr0),P)

  # Get design matrices
  X = list()
  X$clearance =
    model.matrix(as.formula(paste0("~ ",
                                   paste(X_vars$clearance,collapse = "+"))),
                 data = data_gr0)
  X$incidence =
    model.matrix(as.formula(paste0("~ ",
                                   paste(X_vars$incidence,collapse = "+"))),
                 data = data_gr0)
  X$prevalence =
    model.matrix(as.formula(paste0("~ ",
                                   paste(X_vars$prevalence,collapse = "+"))),
                 data = data_0)

  # Get y_{itp} matrix
  y_i0p = z_i0p =
    data_0 %>%
    select(all_of(pathogen_vars)) %>%
    as.matrix()
  y_itp = z_itp =
    data_gr0 %>%
    select(all_of(pathogen_vars)) %>%
    as.matrix()

  # Get lagged
  get_lagged_z = function(dataset,z_i0p,z_itp){

    dataset[dataset$time0,grep("_est",colnames(dataset))] =
      z_i0p
    dataset[!dataset$time0,grep("_est",colnames(dataset))] =
      z_itp

    for(p in 1:P){
      dataset[[paste0(pathogen_vars[p],"_lagged")]] =
        ifelse(!dataset$time0,
               c(-999,dataset[[paste0(pathogen_vars[p],"_est")]][-NT]),
               NA)
    }

    dataset %>%
      na.omit() %>%
      select(all_of(paste(pathogen_vars,"lagged",sep="_"))) %>%
      as.matrix()
  }
  z_lagged = get_lagged_z(data_aug,z_i0p,z_itp)


  # Initialize --------------------------------------------------------------

  if(verbose) cat("\nInitializing parameters\n")
  fits = list()

  ## Clearance
  fits$clearance = list()
  data_aug_reverse_coded =
    data_aug
  for(p in 1:P){
    data_aug_reverse_coded[[pathogen_vars[p]]][clearance_index[[p]]] =
      1 - data_aug_reverse_coded[[pathogen_vars[p]]][clearance_index[[p]]]
  }
  for(p in valid_responses$combined){
    formula_p = NULL

    try({
      formula_p =
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(c(X_vars$clearance,
                  paste(pathogen_vars[valid_covariates$clearance[[p]]],
                        "_lagged",
                        sep = "")),
                collapse = " + ")
        ))
    },silent = TRUE)

    if(is.null(formula_p)){
      formula_p =
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(X_vars$clearance,
                collapse = " + ")
        ))
    }

    try({
      fits$clearance[[p]] =
        glm(formula_p,
            data =
              data_aug_reverse_coded[clearance_index[[p]],],
            family = binomial("cloglog"))
    },silent = TRUE)
  }

  ## Fit incidence models
  fits$incidence = list()
  for(p in valid_responses$combined){
    formula_p = NULL

    try({
      formula_p =
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(c(X_vars$incidence,
                  paste(pathogen_vars[valid_covariates$incidence[[p]]],
                        "_lagged",
                        sep = "")),
                collapse = " + ")
        ))
    },silent = TRUE)

    if(is.null(formula_p)){
      formula_p =
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(X_vars$incidence,
                collapse = " + ")
        ))
    }

    try({
      fits$incidence[[p]] =
        glm(formula_p,
            data = data_aug[incidence_index[[p]],],
            family = binomial("cloglog"))
    },silent = TRUE)
  }

  ## Fit prevalence models
  fits$prevalence = list()
  for(p in valid_responses$prevalence){
    formula_p =
      as.formula(paste0(
        pathogen_vars[p],
        " ~ ",
        paste(X_vars$prevalence,
              collapse = " + ")
      ))
    try({
      fits$prevalence[[p]] =
        glm(formula_p,
            data = data_0,
            family = binomial())
    },silent = TRUE)
  }

  ## Fill in values from glm's
  beta_C = array(0.0,
                 c(Q$clearance,P,n_draws),
                 dimnames = list(c("(Intercept)",X_vars$clearance),
                                 pathogen_vars,
                                 NULL))
  beta_I = array(0.0,
                 c(Q$incidence,P,n_draws),
                 dimnames = list(c("(Intercept)",X_vars$incidence),
                                 pathogen_vars,
                                 NULL))
  beta_pr = array(0.0,
                  c(Q$prevalence,P,n_draws),
                  dimnames = list(c("(Intercept)",X_vars$prevalence),
                                  pathogen_vars,
                                  NULL))
  eta_I = eta_C = array(0.0,
                        c(P,P,n_draws),
                        dimnames = list(paste(pathogen_vars,"lagged",sep="_"),
                                        pathogen_vars,
                                        NULL))
  for(p in 1:P){
    try({
      if(p %in% valid_responses$combined){
        beta_C[names(coef(fits$clearance[[p]]))[1:Q$clearance],p,1] =
          coef(fits$clearance[[p]])[1:Q$clearance]
        try({
          eta_C[names(coef(fits$clearance[[p]]))[-c(1:Q$clearance)],p,1] =
            coef(fits$clearance[[p]])[-c(1:Q$clearance)]
        },silent=TRUE)
      }
    },silent=TRUE)

    try({
      if(p %in% valid_responses$combined){
        beta_I[names(coef(fits$incidence[[p]]))[1:Q$incidence],p,1] =
          coef(fits$incidence[[p]])[1:Q$incidence]
        try({
          eta_I[names(coef(fits$incidence[[p]]))[-c(1:Q$incidence)],p,1] =
            coef(fits$incidence[[p]])[-c(1:Q$incidence)]
        },silent=TRUE)
      }
    },silent=TRUE)

    try({
      if(p %in% valid_responses$prevalence){
        beta_pr[,p,1] =
          coef(fits$prevalence[[p]])
      }
    },silent=TRUE)
  }


  ## Finally, do sensitivity and specificity
  Se = Sp = numeric(n_draws)

  Se[1] = prior_sensitivity$a / sum(unlist(prior_sensitivity))
  Sp[1] = prior_specificity$a / sum(unlist(prior_specificity))



  # Set prior hyperparameter values -----------------------------------------

  sd_x = lapply(X,function(x) apply(as.matrix(x[,-1]),2,sd))
  names(sd_x) = names(X)
  sd_zlagged = apply(z_lagged,2,sd)
  prior_list = list(beta_I = list(mean = prior_beta_I$location,
                                  sd =
                                    prior_beta_I$scale *
                                    c(1,sd_x$incidence)^prior_beta_I$autoscale),
                    beta_C = list(mean = prior_beta_C$location,
                                  sd =
                                    prior_beta_C$scale *
                                    c(1,sd_x$clearance)^prior_beta_C$autoscale),
                    beta_pr = list(mean = prior_prev$location,
                                   sd =
                                     prior_prev$scale *
                                     c(1,sd_x$prevalence)^prior_prev$autoscale),
                    eta_C = list(mean = prior_eta_C$location,
                                 sd =
                                   prior_eta_C$scale *
                                   sd_zlagged^prior_eta_C$autoscale),
                    eta_I = list(mean = prior_eta_I$location,
                                 sd =
                                   prior_eta_I$scale *
                                   sd_zlagged^prior_eta_I$autoscale))
  rm(sd_x,sd_zlagged)


  # Create functions for MH-within-Gibbs sampler ----------------------------

  get_tau_lambda = function(z_lagged,beta_I,beta_C,eta_I,eta_C){
    log_lambda =
      (1 - z_lagged) * (X$incidence %*% beta_I + z_lagged %*% eta_I) +
      (z_lagged) * (X$clearance %*% beta_C + z_lagged %*% eta_C)

    tau_matrix * exp(log_lambda)
  }
  tau_lambda = get_tau_lambda(z_lagged = z_lagged,
                              beta_I = beta_I[,,1],
                              beta_C = beta_C[,,1],
                              eta_I = eta_I[,,1],
                              eta_C = eta_C[,,1])


  logpost_prevalence = function(z_i0p,rho,Se,Sp){

    sum(z_i0p[,valid_responses$prevalence] *
          log(rho[,valid_responses$prevalence])) +
      sum((1 - z_i0p[,valid_responses$prevalence]) *
            log(1-rho[,valid_responses$prevalence])) +
      sum(z_i0p[,valid_responses$prevalence] *
            y_i0p[,valid_responses$prevalence]) * log(Se) +
      sum(z_i0p[,valid_responses$prevalence] *
            (1 - y_i0p[,valid_responses$prevalence])) * log(1 - Se) +
      sum((1 - z_i0p[,valid_responses$prevalence]) *
            (1 - y_i0p[,valid_responses$prevalence])) * log(Sp) +
      sum((1 - z_i0p[,valid_responses$prevalence]) *
            y_i0p[,valid_responses$prevalence]) * log(1 - Sp)
  }

  logpost_clear_incid = function(tau_lambda,z_itp,z_lagged,Se,Sp){
    z_star = abs(z_itp - z_lagged)

    sum(z_itp[,valid_responses$combined] *
          y_itp[,valid_responses$combined]) * log(Se) +
      sum(z_itp[,valid_responses$combined] *
            (1 - y_itp[,valid_responses$combined])) * log(1 - Se) +
      sum((1 - z_itp[,valid_responses$combined]) *
            (1 - y_itp[,valid_responses$combined])) * log(Sp) +
      sum((1 - z_itp[,valid_responses$combined]) *
            y_itp[,valid_responses$combined]) * log(1 - Sp) +
      sum(z_star[,valid_responses$combined] *
            log(1 - exp(-tau_lambda[,valid_responses$combined]))) -
      sum((1 - z_star[,valid_responses$combined]) *
            tau_lambda[,valid_responses$combined])
  }



  prior_logpdf = function(beta_I,beta_C,beta_pr,eta_C,eta_I){
    sum(dnorm(c(beta_I),
              mean = c(prior_list$beta_I$mean),
              sd = c(prior_list$beta_I$sd),
              log = TRUE)) +
      sum(dnorm(c(beta_C),
                mean = c(prior_list$beta_C$mean),
                sd = c(prior_list$beta_C$sd),
                log = TRUE)) +
      sum(dnorm(c(beta_pr),
                mean = c(prior_list$beta_pr$mean),
                sd = c(prior_list$beta_pr$sd),
                log = TRUE)) +
      sum(dnorm(c(eta_C),
                mean = c(prior_list$eta_C$mean),
                sd = c(prior_list$eta_C$sd),
                log = TRUE)) +
      sum(dnorm(c(eta_I),
                mean = c(prior_list$eta_I$mean),
                sd = c(prior_list$eta_I$sd),
                log = TRUE))
  }


  draw_beta_pr = expression({
    new_beta_pr =
      beta_pr[,,iter]
    for(p in valid_responses$prevalence){
      new_beta_pr[,p] =
        new_beta_pr[,p] +
        rnorm(Q$prevalence) %*% MH$beta_pr[[p]] * MH_scalars$prevalence
    }

    rho_new = rho

    for(p in valid_responses$prevalence){
      rho_new[,p] =
        1 / (1 + exp(-X$prevalence %*% beta_pr[,p,1]))

      acc_probability =
        logpost_prevalence(z_i0p,rho_new,Se[iter],Sp[iter]) -
        logpost_prevalence(z_i0p,rho,Se[iter],Sp[iter]) +
        prior_logpdf(beta_I[,,iter],beta_C[,,iter],new_beta_pr,eta_C[,,iter],eta_I[,,iter]) -
        prior_logpdf(beta_I[,,iter],beta_C[,,iter],beta_pr[,,iter],eta_C[,,iter],eta_I[,,iter])

      if(runif(1) < exp(acc_probability)){
        acc_rates$prevalence[p] =
          acc_rates$prevalence[p] + 1L

        rho[,p] = rho_new[,p]
        beta_pr[,p,iter] =
          new_beta_pr[,p]
      }else{
        rho_new[,p] = rho[,p]
      }

    }
  })


  draw_beta_eta_I = expression({
    tau_lambda = get_tau_lambda(z_lagged = z_lagged,
                                beta_I = beta_I[,,iter],
                                beta_C = beta_C[,,iter],
                                eta_I = eta_I[,,iter],
                                eta_C = eta_C[,,iter])
    new_beta_eta =
      rbind(beta_I[,,iter],
            eta_I[,,iter])
    for(p in valid_responses$combined){
      new_beta_eta[get_rows(p,"incidence"),p] =
        new_beta_eta[get_rows(p,"incidence"),p] +
        rnorm(Q$incidence + length(valid_covariates$incidence[[p]])) %*%
        MH$beta_eta_I[[p]] * MH_scalars$incidence

      new_tau_lambda = get_tau_lambda(z_lagged = z_lagged,
                                      beta_I = new_beta_eta[1:Q$incidence,],
                                      beta_C = beta_C[,,iter],
                                      eta_I = new_beta_eta[-c(1:Q$incidence),],
                                      eta_C = eta_C[,,iter])

      acc_probability =
        logpost_clear_incid(new_tau_lambda,z_itp,z_lagged,Se[iter],Sp[iter]) -
        logpost_clear_incid(tau_lambda,z_itp,z_lagged,Se[iter],Sp[iter]) +
        prior_logpdf(new_beta_eta[1:Q$incidence,],beta_C[,,iter],beta_pr[,,iter],eta_C[,,iter],new_beta_eta[-c(1:Q$incidence),]) -
        prior_logpdf(beta_I[,,iter],beta_C[,,iter],beta_pr[,,iter],eta_C[,,iter],eta_I[,,iter])


      if(runif(1) < exp(acc_probability)){
        acc_rates$incidence[p] =
          acc_rates$incidence[p] + 1L

        tau_lambda = new_tau_lambda
        beta_I[,p,iter] =
          new_beta_eta[1:Q$incidence,p]
        eta_I[,p,iter] =
          new_beta_eta[-c(1:Q$incidence),p]
      }else{
        new_beta_eta[,p] =
          c(beta_I[,p,iter],
            eta_I[,p,iter])
      }
    }
  })


  draw_beta_eta_C = expression({
    tau_lambda = get_tau_lambda(z_lagged = z_lagged,
                                beta_I = beta_I[,,iter],
                                beta_C = beta_C[,,iter],
                                eta_I = eta_I[,,iter],
                                eta_C = eta_C[,,iter])
    new_beta_eta =
      rbind(beta_C[,,iter],
            eta_C[,,iter])
    for(p in valid_responses$combined){
      new_beta_eta[get_rows(p,"clearance"),p] =
        new_beta_eta[get_rows(p,"clearance"),p] +
        rnorm(Q$clearance + length(valid_covariates$c[[p]])) %*%
        MH$beta_eta_C[[p]] * MH_scalars$clearance

      new_tau_lambda = get_tau_lambda(z_lagged = z_lagged,
                                      beta_I = beta_I[,,iter],
                                      beta_C = new_beta_eta[1:Q$clearance,],
                                      eta_I = eta_I[,,iter],
                                      eta_C = new_beta_eta[-c(1:Q$clearance),])

      acc_probability =
        logpost_clear_incid(new_tau_lambda,z_itp,z_lagged,Se[iter],Sp[iter]) -
        logpost_clear_incid(tau_lambda,z_itp,z_lagged,Se[iter],Sp[iter]) +
        prior_logpdf(beta_I[,,iter],new_beta_eta[1:Q$clearance,],beta_pr[,,iter],new_beta_eta[-c(1:Q$clearance),],eta_I[,,iter]) -
        prior_logpdf(beta_I[,,iter],beta_C[,,iter],beta_pr[,,iter],eta_C[,,iter],eta_I[,,iter])


      if(runif(1) < exp(acc_probability)){
        acc_rates$c[p] =
          acc_rates$clearance[p] + 1L

        tau_lambda = new_tau_lambda
        beta_C[,p,iter] =
          new_beta_eta[1:Q$clearance,p]
        eta_C[,p,iter] =
          new_beta_eta[-c(1:Q$clearance),p]
      }else{
        new_beta_eta[,p] =
          c(beta_C[,p,iter],
            eta_C[,p,iter])
      }
    }
  })


  draw_z0 = expression({
    prob_matrix =
      rho * (Se[iter]^y_i0p * (1 - Se[iter])^(1 - y_i0p))
    prob_matrix =
      prob_matrix /
      (prob_matrix +
         (1 - rho) * (Sp[iter]^(1 - y_i0p) * (1 - Sp[iter])^y_i0p))
    for(p in valid_responses$prevalence){
      z_i0p[,p] =
        rbinom(nrow(z_i0p),1,prob_matrix[,p])
    }

    z_i0_estimate =
      z_i0_estimate + z_i0p
  })

  draw_zt = expression({
    e_neg_tau_lambda = exp(-tau_lambda)
    prob_matrix_new =
      MH_scalars$z_itp * Se[iter]^y_itp * (1 - Se[iter])^(1 - y_itp) *
      (1 - e_neg_tau_lambda)^(1 - z_lagged) * e_neg_tau_lambda^z_lagged
    prob_matrix_new =
      prob_matrix_new /
      (prob_matrix_new +
         Sp[iter]^(1 - y_itp) * (1 - Sp[iter])^y_itp *
         (1 - e_neg_tau_lambda)^z_lagged * e_neg_tau_lambda^(1 - z_lagged) )

    z_itp_new = z_itp
    for(p in valid_responses$combined){
      z_itp_new[,p] =
        rbinom(nrow(z_itp),1,prob_matrix_new[,p])
    }

    z_lagged_new =
      get_lagged_z(data_aug,z_i0p,z_itp_new)
    tau_lambda_new = get_tau_lambda(z_lagged_new,beta_I[,,iter],beta_C[,,iter],eta_I[,,iter],eta_C[,,iter])

    e_neg_tau_lambda_new = exp(-tau_lambda_new)
    prob_matrix_new =
      MH_scalars$z_itp * Se[iter]^y_itp * (1 - Se[iter])^(1 - y_itp) *
      (1 - e_neg_tau_lambda_new)^(1 - z_lagged_new) * e_neg_tau_lambda_new^z_lagged_new
    prob_matrix_new =
      prob_matrix_new /
      (prob_matrix_new +
         Sp[iter]^(1 - y_itp) * (1 - Sp[iter])^y_itp *
         (1 - e_neg_tau_lambda_new)^z_lagged_new * e_neg_tau_lambda_new^(1 - z_lagged_new) )

    acc_probability =
      logpost_clear_incid(tau_lambda_new,z_itp_new,z_lagged_new,Se[iter],Sp[iter]) -
      logpost_clear_incid(tau_lambda,z_itp,z_lagged,Se[iter],Sp[iter]) +
      sum(dbinom(c(z_itp[,valid_responses$combined]),
                 1,
                 c(prob_matrix_new[,valid_responses$combined]),
                 log = TRUE)) -
      sum(dbinom(c(z_itp_new[,valid_responses$combined]),
                 1,
                 c(prob_matrix[,valid_responses$combined]),
                 log = TRUE))

    if(runif(1) < exp(acc_probability)){
      acc_rates$Z_t = acc_rates$Z_t + 1L
      z_itp = z_itp_new
      z_lagged = z_lagged_new
      tau_lambda = tau_lambda_new
    }

    z_it_estimate =
      z_it_estimate + z_itp
  })


  draw_sens_spec = expression({
    one_minus_y =
      list(1 - y_i0p,
           1 - y_itp)
    one_minus_z =
      list(1 - z_i0p,
           1 - z_itp)
    a_spec =
      prior_specificity$a +
      sum(one_minus_y[[1]] * one_minus_z[[1]]) +
      sum(one_minus_y[[2]] * one_minus_z[[2]])
    b_spec =
      prior_specificity$b +
      sum(y_i0p * one_minus_z[[1]]) +
      sum(y_itp * one_minus_z[[2]])

    a_sens =
      prior_sensitivity$a +
      sum(y_i0p * z_i0p) +
      sum(y_itp * z_itp)
    b_sens =
      prior_sensitivity$b +
      sum(one_minus_y[[1]] * z_i0p) +
      sum(one_minus_y[[2]] * z_itp)

    Sp[iter] = rbeta(1,a_spec,b_spec)
    Se[iter] = rbeta(1,a_sens,b_sens)
  })



  # Create objects for MCMC  ------------------------------------------------


  MH = list()
  MH$beta_eta_C =
    lapply(fits$clearance,
           function(x) try({chol(summary(x)$cov.scaled)},silent=T))
  MH$beta_eta_I =
    lapply(fits$incidence,
           function(x) try({chol(summary(x)$cov.scaled)},silent=T))
  MH$beta_pr =
    lapply(fits$prevalence,
           function(x) try({chol(summary(x)$cov.scaled)},silent=T))

  # Check for nonsense proposals from bad initial unregularized fits
  prior_determinants = proposal_determinants = list()
  proposal_determinants$beta_eta_C =
    sapply(MH$beta_eta_C, function(x) try({det(crossprod(x))},silent=T))
  proposal_determinants$beta_eta_I =
    sapply(MH$beta_eta_I, function(x) try({det(crossprod(x))},silent=T))
  proposal_determinants$beta_pr =
    sapply(MH$beta_pr, function(x) try({det(crossprod(x))},silent=T))
  prior_determinants$beta_eta_C =
    sapply(1:P,function(p) try({prod(c(prior_list$beta_C$sd^2,
                                       prior_list$eta_C$sd[ valid_covariates$clearance[[p]] ]^2))},
                               silent=T))
  prior_determinants$beta_eta_I =
    sapply(1:P,function(p) try({prod(c(prior_list$beta_I$sd^2,
                                       prior_list$eta_I$sd[ valid_covariates$incidence[[p]] ]^2))},
                               silent=T))
  prior_determinants$beta_pr =
    prod(prior_list$beta_pr$sd^2)
  # Left off here
  for(p in valid_responses$combined){
    if(proposal_determinants$beta_eta_C[[p]] >= prior_determinants$beta_eta_C){
      MH$beta_eta_C[[p]] =
        sapply()c(prior_list$beta_C$sd,prior_list$eta_C$sd)

    }
  }



  acc_rates = list()
  acc_rates$clearance =
    acc_rates$incidence =
    acc_rates$prevalence =
    integer(P)
  acc_rates$Z_t =
    0L

  rho = 1 / (1 + exp(-X$prevalence %*% beta_pr[,,1]))

  get_rows = function(p,type){
    if(is.null(valid_covariates[[type]][[p]])){
      c("(Intercept)",X_vars[[type]])
    }else{
      c("(Intercept)",X_vars[[type]],
        paste(pathogen_vars[valid_covariates[[type]][[p]]],"lagged",sep="_"))
    }
  }

  z_i0_estimate = matrix(0.0,nrow(z_i0p),P)
  z_it_estimate = matrix(0.0,nrow(z_itp),P)



  # Perform MH-within-Gibbs sampler -----------------------------------------

  if(verbose) pb = txtProgressBar(0,n_draws,style=3)
  for(iter in 2:n_draws){
    beta_C[,,iter] = beta_C[,,iter - 1]
    beta_I[,,iter] = beta_I[,,iter - 1]
    eta_C[,,iter] = eta_C[,,iter - 1]
    eta_I[,,iter] = eta_I[,,iter - 1]
    beta_pr[,,iter] = beta_pr[,,iter - 1]
    Se[iter] = Se[iter - 1]
    Sp[iter] = Sp[iter - 1]

    eval(draw_beta_pr)
    eval(draw_beta_eta_I)
    eval(draw_beta_eta_C)
    eval(draw_z0)
    eval(draw_zt)
    eval(draw_sens_spec)

    setTxtProgressBar(pb,iter)
  }


}
