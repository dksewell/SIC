#'



if(FALSE){
  library(magrittr)
  library(rstanarm)
  library(dplyr)
  library(plot.matrix)
  library(viridis)
  source("~/SIC/R/sic_simulator.R")
  sic_data = 
    sic_simulator(seed = 2023)
  
  test = sic(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) ~ x1 + x2 +(time | subject),
             sic_data$data[[1]] %>% select(-tau),
             prior = student_t(df = 5, scale = 2.5),
             verbose = TRUE,
             seed = 2023)
}

plot.sic = function(x,
                    type = c("trace","matrix")[1]){
  
  #----
  if(type == "trace"){
    cat("\nIncidence models:\n")
    
    for(i in 1:length(x$estimability$response$incidence)){
      cat(paste0("---",x$estimability$response$incidence[i],"\n"))
      
      print(
        plot(x$stanreg_fits$incidence[[i]],
             plotfun = "trace") +
          ggplot2::ggtitle(paste0(x$estimability$response$incidence[i],
                                  " (Incidence)"))
      )
      
      cat("Hit ENTER for next plot\n")
      readline()
    }
    
    for(i in 1:length(x$estimability$response$clearance)){
      cat(paste0("---",x$estimability$response$clearance[i],"\n"))
      
      print(
        plot(x$stanreg_fits$clearance[[i]],
             plotfun = "trace") +
          ggplot2::ggtitle(paste0(x$estimability$response$clearance[i],
                                  " (Clearance)"))
      )
      
      cat("Hit ENTER for next plot\n")
      readline()
    }
    
  }
  
  
  #----
  if(type == "matrix"){
    plot(x$networks$incidence,
         main = "Incidence",
         col = viridis::viridis)
         # col = colorRampPalette(c("darkslategray3","darkslategray4","brown4","brown1")))
    plot(x$networks$clearance,
         main = "Clearance",
         col = function(x) viridis(n = x,begin=1,end=0))
         # col = colorRampPalette(rev(c("darkslategray3","darkslategray4","brown4","brown1"))))
  }
  
}