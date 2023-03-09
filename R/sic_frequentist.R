#' Perform frequentist inference for the SIC model 
#' 
#' Compute likelihood based frequentist inference on data fit to
#' the SIC model.  This will return output for the incidence and 
#' the clearance models separately.
#' 
#' @param formula Formula of the form <y> ~ <x1> + <x2> + ... + (<time variable> | <subject ID>)
#' @param data Data frame or tibble
#' @param min_count_to_estimate positive integer giving the minimum cell count in a 2x2 
#' contingency table required to regress one pathogen on another
#' @param CI_level numeric between 0 and 1 for the confidence interval level.
#' @returns  Object of class 'sic_freq' which has the following elements:
#' * results tibble with columns telling the response variable, the model (incidence or 
#' clearance), the covariate, the regression coefficient estimate, the lower and upper CI 
#' bounds, and (begrudgingly) the 2-sided p-value.
#' *  logLik the log likelihood of the fitted model.  Note that this only includes the 
#' model fits for those response variables and corresponding covariates that met the 
#' min_count_to_estimate criterion.
#' * estimability list with the valid response variables and the valid covariates 
#' for each response variable
#' * networks list of the PxP weighted adjacency matrix giving the regression coefficients 
#' as the edge weights.  E.g., the (i,j)-th element of incidence matrix is the effect of 
#' the presence of pathogen i on the log incidence rate of pathogen j.
#' * CI_level 
#' * min_count_to_estimate
#' 
#' @examples 
#' sic_data = sic_simulator(seed = 2023)
#' test = sic_frequentist(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject),
#'                        sic_data$data[[1]] %>% dplyr::select(-tau))
#' 
#' @export
#' 

if(FALSE){
  library(magrittr)
  library(dplyr)
  source("~/SIC/R/sic_simulator.R")
  sic_data = 
    sic_simulator(seed = 1)
  # data = 
  #   sic_data$data[[1]] %>% 
  #   dplyr::select(-tau)
  # 
  # formula = cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject)
  # 
  # min_count_to_estimate = 5
  # CI_level = 0.95
  
  test = sic_frequentist(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject),
                         sic_data$data[[1]] %>% 
                           dplyr::select(-tau))
                         
}

sic_frequentist = function(formula,
                           data, 
                           min_count_to_estimate = 5,
                           CI_level = 0.95){
  
  # Drop NAs
  data %<>% 
    dplyr::select(all.vars(formula)) %>% 
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
    dplyr::select(all_of(subject_var)) %>% 
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
  
  
  # Compute tables to see if we can make inference
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
    
    if( (length(unique(data[[pathogen_vars[p]]][ incidence_index[[p]] ] )) > 1) &
        (min(table(data[[pathogen_vars[p]]][ incidence_index[[p]]])) > min_count_to_estimate ) ){
      valid_responses$incidence = 
        c(valid_responses$incidence,p)
      
      no_valid_covars = TRUE
      for(p2 in c(1:P)[-p]){
        tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]] = 
          table(data[[pathogen_vars[p]]][ incidence_index[[p]] ] ,
                data[[paste0(pathogen_vars[p2],"_lagged")]][ incidence_index[[p]] ] )
        
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
    
    if( (length(unique(data[[pathogen_vars[p]]][ clearance_index[[p]] ] )) > 1) &
        (min(table(data[[pathogen_vars[p]]][ clearance_index[[p]]])) > min_count_to_estimate ) ){
      valid_responses$clearance = 
        c(valid_responses$clearance,p)
      
      no_valid_covars = TRUE
      for(p2 in c(1:P)[-p]){
        tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]] = 
          table(data[[pathogen_vars[p]]][ clearance_index[[p]] ] ,
                data[[paste0(pathogen_vars[p2],"_lagged")]][ clearance_index[[p]] ] )
        
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
  
  
  # Reverse code clearance pathogen variables
  for(p in 1:P){
    data[[pathogen_vars[p]]][clearance_index[[p]]] = 
      1 - data[[pathogen_vars[p]]][clearance_index[[p]]]
  }
  
  
  # Get model fits
  fits = list()
  fits$incidence = 
    fits$clearance = list()
  ## Fit incidence models
  for(p in valid_responses$incidence){
    formula_p = NULL
    
    try({
      formula_p = 
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(c(X_vars,
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
          paste(X_vars,
                collapse = " + ")
        ))
    }
    
    try({ 
      fits$incid[[p]] =
        glm(formula_p,
            data = data[incidence_index[[p]],],
            family = binomial("cloglog"))
    },silent = TRUE)
  }
  ## Fit clearance models
  for(p in valid_responses$clearance){
    formula_p = NULL
    
    try({
      formula_p = 
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(c(X_vars,
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
          paste(X_vars,
                collapse = " + ")
        ))
    }
     
    try({ 
      fits$clearance[[p]] =
        glm(formula_p,
            data = data[clearance_index[[p]],],
            family = binomial("cloglog"))
    },silent = TRUE)
  }
  
  
  # Collect results
  helper_fun = function(x){
    temp = summary(x)$coef
    cbind(Estimate = temp[,"Estimate"],p = temp[,"Pr(>|z|)"])
  }
  estimates_pvals = 
    list(incidence = lapply(fits$incidence[which(!sapply(fits$incidence,is.null))],helper_fun),
         clearance = lapply(fits$clearance[which(!sapply(fits$clearance,is.null))],helper_fun)
    )
  
  CIs = 
    list(incidence = lapply(fits$incidence[which(!sapply(fits$incidence,is.null))],confint,level = CI_level),
         clearance = lapply(fits$clearance[which(!sapply(fits$clearance,is.null))],confint,level = CI_level)
         )
  
  names(estimates_pvals$incidence) = 
    names(CIs$incidence) = 
    pathogen_vars[valid_responses$incidence]
  names(estimates_pvals$clearance) = 
    names(CIs$clearance) = 
    pathogen_vars[valid_responses$clearance]
  
  
  results = 
    tibble(Response = character(),
           Model = character(),
           Covariate = character(),
           Estimate = numeric(),
           Lower = numeric(),
           Upper = numeric(),
           `p-value` = numeric())
  
  for(i in 1:P){
    if(i %in% valid_responses$incidence){
      results %<>%
        bind_rows(
          tibble(Response = pathogen_vars[i],
                 Model = "Incidence",
                 Covariate = rownames(CIs$incidence[[pathogen_vars[i]]]),
                 Estimate = estimates_pvals$incidence[[pathogen_vars[i]]][,1],
                 Lower = CIs$incidence[[pathogen_vars[i]]][,1],
                 Upper = CIs$incidence[[pathogen_vars[i]]][,2],
                 `p-value` = estimates_pvals$incidence[[pathogen_vars[i]]][,2])
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
                 `p-value` = NA)
        )
    }
    
    if(i %in% valid_responses$clearance){
      results %<>%
        bind_rows(
          tibble(Response = pathogen_vars[i],
                 Model = "Clearance",
                 Covariate = rownames(CIs$clearance[[pathogen_vars[i]]]),
                 Estimate = estimates_pvals$clearance[[pathogen_vars[i]]][,1],
                 Lower = CIs$clearance[[pathogen_vars[i]]][,1],
                 Upper = CIs$clearance[[pathogen_vars[i]]][,2],
                 `p-value` = estimates_pvals$clearance[[pathogen_vars[i]]][,2])
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
                 `p-value` = NA)
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
    dplyr::select("Response","Covariate","Estimate") %>% 
    mutate(Covariate = gsub("_lagged","",Covariate)) %>% 
    filter(Covariate %in% pathogen_vars)
  incidence_network[cbind(temp1$Covariate,
                          temp1$Response)] = temp1$Estimate
  temp1 = 
    results %>% 
    filter(Model == "Clearance") %>% 
    dplyr::select("Response","Covariate","Estimate") %>% 
    mutate(Covariate = gsub("_lagged","",Covariate)) %>% 
    filter(Covariate %in% pathogen_vars)
  clearance_network[cbind(temp1$Covariate,
                          temp1$Response)] = temp1$Estimate
  
  # Store log likelihood
  ll = 
    sum(c(sapply(fits$incidence[which(!sapply(fits$incidence,is.null))],logLik),
          sapply(fits$clearance[which(!sapply(fits$clearance,is.null))],logLik))
    )
  
  
  valid_responses$incidence = 
    pathogen_vars[valid_responses$incidence]
  valid_responses$clearance = 
    pathogen_vars[valid_responses$clearance]
  for(i in 1:P){
    try({ 
      valid_covariates$incidence[[i]] = pathogen_vars[valid_covariates$incidence[[i]]]
      names(valid_covariates$incidence)[i] = pathogen_vars[i]
      }, silent = T)
    try({
      valid_covariates$clearance[[i]] = pathogen_vars[valid_covariates$clearance[[i]]]
      names(valid_covariates$clearance)[i] = pathogen_vars[i]
      }, silent = T)
  }
  
  
  object = 
    list(results = results,
         logLik = ll,
         estimability = list(response = valid_responses,
                             covariates = valid_covariates),
         networks = list(incidence = incidence_network,
                         clearance = clearance_network),
         CI_level = CI_level,
         min_count_to_estimate = min_count_to_estimate)
  
  class(object) = "sic_freq"
  
  return(object)
}













