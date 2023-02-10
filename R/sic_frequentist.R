#' Perform frequentist inference for the SIC model 
#' 
#' Compute likelihood based frequentist inference on data fit to
#' the SIC model.  This will return output for the incidence and 
#' the clearance models separately.
#' 

if(FALSE){
  library(magrittr)
  library(dplyr)
  source("~/SIC/R/sic_simulator.R")
  sic_data = 
    sic_simulator(seed = 2023)
  data = 
    sic_data$data[[1]] %>% 
    select(-tau)
  
  X = model.matrix(formula,data)
  X = X[,-ncol(X)]
  
  formula = cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject)
  
}

sic_frequentist = function(formula,data){
  
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
    
    if(length(unique(data[[pathogen_vars[p]]][ incidence_index[[p]] ] )) > 1){
      valid_responses$incidence = 
        c(valid_responses$incidence,p)
      
      for(p2 in c(1:P)[-p]){
        tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]] = 
          table(data[[pathogen_vars[p]]][ incidence_index[[p]] ] ,
                data[[paste0(pathogen_vars[p2],"_lagged")]][ incidence_index[[p]] ] )
        
        names(dimnames(tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]] )) = 
          pathogen_vars[c(p,p2)]
        
        if(ncol(tables$incidence[[pathogen_vars[p]]][[pathogen_vars[p2]]]) > 1){
          if(p2 > min(c(1:P)[-p])){
            valid_covariates$incidence[[ p ]] =
              c(valid_covariates$incidence[[ p ]],
                p2)
          }else{
            valid_covariates$incidence[[ p ]] = p2
          }
        }
        
      }
        
    }
    
    if(length(unique(data[[pathogen_vars[p]]][ clearance_index[[p]] ] )) > 1){
      valid_responses$clearance = 
        c(valid_responses$clearance,p)
      
      for(p2 in c(1:P)[-p]){
        tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]] = 
          table(data[[pathogen_vars[p]]][ clearance_index[[p]] ] ,
                data[[paste0(pathogen_vars[p2],"_lagged")]][ clearance_index[[p]] ] )
        
        names(dimnames(tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]] )) = 
          pathogen_vars[c(p,p2)]
        
        if(ncol(tables$clearance[[pathogen_vars[p]]][[pathogen_vars[p2]]]) > 1){
          if(p2 > min(c(1:P)[-p])){
            valid_covariates$clearance[[ p ]] =
              c(valid_covariates$clearance[[ p ]],
                p2)
          }else{
            valid_covariates$clearance[[ p ]] = p2
          }
        }
        
      }
      
    }
    
  }
  
  stop("I stopped here.  Need to move forward taking into account which variables can be modeled with which other covariates")
  
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
    try({
      formula_p = 
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(c(X_vars,
                  paste(pathogen_vars[valid_covariates$incidence[[ ]]],
                        "_lagged",
                        sep = "")),
                collapse = " + ")
          ))
      fits$incidence[[p]] =
        glm(formula_p,
            data = data[incidence_index[[p]],],
            family = binomial("cloglog"))
    },silent = TRUE)
  }
  ## Fit clearance models
  for(p in 1:P){
    try({
      formula_p = 
        as.formula(paste0(
          pathogen_vars[p],
          " ~ ",
          paste(c(X_vars,
                  paste(pathogen_vars,"_lagged",sep = "")),
                collapse = " + ")
        ))
      fits$clearance[[p]] =
        glm(formula_p,
            data = data[clearance_index[[p]],],
            family = binomial("cloglog"))
    },silent = TRUE)
  }
  
  
  # Collect results
  ## \beta_I
  beta_I = sapply(fits$incidence,coef)
  
  
  # Store log likelihood
  ll = 
    sum(c(sapply(fits$incidence,logLik),
          sapply(fits$clearance,logLik))
    )
  
}













