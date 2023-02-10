#' Simulate data from the SIC model
#' 
#' @param nsim integer. The number of simulated data sets to create.
#' @param seed a single value, interpreted as an integer.  Used for replication purposes.
#' @param P positive integer.  Number of pathogen types
#' @param N_subj positive integer. Number of subjects in data sets
#' @param N_time Either positive integer or vector of integers of length N_subj. 
#' Number of time points per subject.
#' @param bl_prevalences positive numeric vector of length P, giving the 
#' baseline prevalences for each pathogen type.
#' @param bl_clearance_rates positive numeric vector of length P, giving the 
#' overall rate of clearance for each pathogen type
#' 
#' @export

sic_simulator = function(nsim = 1,
                         seed = NULL,
                         P = 10,
                         N_subj = 100,
                         N_time = 5,
                         bl_prevalences = rexp(P,1/0.05),
                         bl_clearance_rates = runif(P,3,10)){
  
  # Set seed for reproducibility
  if(!is.null(seed)){
    set.seed(seed)
  }else{
    warning("Make sure to set your seed!")
  }
  
  # Set the number of time points for each subject, if not done already
  if(length(N_time) == 1) N_time = rep(N_time,N_subj)
  
  # Compute incidence rates
  bl_incidence_rates = bl_prevalences / (1 - bl_prevalences) / bl_clearance_rates
  
  
  # Set up objects to be returned
  parameters = sic_data = list()
  
  # Run nsim simulations
  for(it in 1:nsim){
    
    # Set up tibble with subject and time information
    sic_data[[it]] = 
      tibble(subject = rep(1:N_subj,N_time),
             time = rgamma(sum(N_time),
                           shape = 2,
                           rate = 2 / 7)) %>% 
      arrange(subject,
              time)
    
    # Compute time between observations
    sic_data[[it]] %<>%
      mutate(tau = 
               ifelse(
                 c(TRUE,sic_data[[it]]$subject[-1] == sic_data[[it]]$subject[-nrow(sic_data[[it]])]),
                 diff(c(sic_data[[it]]$time[1],sic_data[[it]]$time)),
                 0.0)
             )
    
    # Set up columns for pathogen presence/absence at each time point for each subject
    for(p in 1:P){
      eval(parse(text = paste0(
        "sic_data[[it]] %<>%
          mutate(p",
        p,
        " = 0L)"
        )))
    }
    
    # Generate one binary and one continuous covariate
    sic_data[[it]] %<>%
      mutate(x1 = rbinom(nrow(sic_data[[it]]),
                         1,
                         0.25),
             x2 = rnorm(nrow(sic_data[[it]])))
    
    # Initialize y using baseline prevalences
    ind = which(sic_data[[it]]$tau == 0)
    for(p in 1:P){
      sic_data[[it]][ind,3 + p] = 
        rbinom(length(ind),
               1,
               bl_prevalences[p])
    }
    
    # Set up parameter values
    beta_I = 
      beta_C = 
      matrix(0.0,3,P)
    
    beta_I[1,] = log(bl_incidence_rates)
    beta_I[-1,] = rnorm(P*2,sd = 2.5)
    beta_C[1,] = log(bl_clearance_rates)
    beta_C[-1,] = rnorm(P*2,sd = 2.5)
    
    eta_I = 
      eta_C =
      matrix(rnorm(P^2,sd=2.5),P,P)
    
    parameters[[it]] = 
      list(beta_I = beta_I,
           beta_C = beta_C,
           eta_I = eta_I,
           eta_C = eta_C)
    
    # Create design matrix
    X =
      sic_data[[it]][,3 + P + 1:2] %>% 
      as.matrix()
    X = cbind(1,X)
    
    # Create helper vector
    prob1 = numeric(P)
    
    # Simulate data for each subject
    for(i in 1:N_subj){
      
      # Pull out appropriate rows of sic_data[[it]]
      i_index = which(sic_data[[it]]$subject == i)
      
      # Find out what their initial infection statuses are
      y_tm1 = as.matrix(sic_data[[it]][i_index[1],3 + 1:P])
      
      # Simulate over each time point
      for(tt in 2:length(i_index)){
        
        # Create rate (start on log scale and then exponentiate)
        lambda_it = 
          (1 - y_tm1) * 
          (
            X[i_index[tt-1],] %*% beta_I
            + y_tm1 %*% eta_I
          ) + 
          y_tm1 * 
          (
            X[i_index[tt-1],] %*% beta_C
            + y_tm1 %*% eta_C
          )
        lambda_it = exp(lambda_it)
        
        # Compute the probability of detecting each pathogen.  Will
        # vary based on the pathogen presence previously
        prob1[!y_tm1] =
          1 - exp(-sic_data[[it]]$tau[i_index[tt]] * 
                    lambda_it[!y_tm1])
        prob1[as.logical(y_tm1)] =
          exp(-sic_data[[it]]$tau[i_index[tt]] * 
                lambda_it[y_tm1])
        
        # Simulate new data at time tt
        y_tm1 = # Yes, this is actually y_t, but it will serve next as y_{t-1}
          rbinom(P,1,prob1)
        sic_data[[it]][i_index[tt],3 + 1:P] = as.list(y_tm1)
      }
      
    }
    
  }
  
  # Return data and stochastically generated parameters
  return(list(data = sic_data,
              parameters = parameters))
}
