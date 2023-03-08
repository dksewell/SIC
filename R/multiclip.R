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
  # library(adaptMCMC)
  source("~/SIC/R/sic_simulator.R")
  source("~/SIC/R/sic_frequentist.R")
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
  prior_eta_I = list(location = 0, scale = sqrt(2.5))
  prior_eta_C = list(location = 0, scale = sqrt(2.5))
  prior_prev = list(location = 0, scale = sqrt(2.5), autoscale = TRUE)
  n_draws = 5e4

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
                     prior_eta_I = list(location = 0, scale = sqrt(2.5)),
                     prior_eta_C = list(location = 0, scale = sqrt(2.5)),
                     prior_prev = list(location = 0, scale = sqrt(2.5), autoscale = TRUE),
                     verbose = TRUE){

  if(!is.null(seed)){
    set.seed(seed)
  }else{
    warning("Make sure to set your seed!")
  }

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

  # Filter data by time = 0 or time >= 1
  data_gr0 =
    data_aug %>%
    na.omit()
  data_0 =
    data_aug %>%
    filter(time0)


  # Determine what is estimable
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
    tables$clearance = list()

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
        (min(table(data_gr0[[pathogen_vars[p]]][ clearance_index[[p]]])) > min_count_to_estimate ) ){
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
  get_lagged_z = function(){

    data_aug[data_aug$time0,grep("_est",colnames(data_aug))] =
      z_i0p
    data_aug[!data_aug$time0,grep("_est",colnames(data_aug))] =
      z_itp

    for(p in 1:P){
      data_aug[[paste0(pathogen_vars[p],"_lagged")]] =
        ifelse(!data_aug$time0,
               c(-999,data_aug[[paste0(pathogen_vars[p],"_est")]][-NT]),
               NA)
    }

    data_aug %>%
      na.omit() %>%
      select(all_of(paste(pathogen_vars,"lagged",sep="_"))) %>%
      as.matrix()
  }
  z_lagged = get_lagged_z()

  # Reverse coding
  z_star = abs(z_itp - z_lagged)



  #asdf
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


}
