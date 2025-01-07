# Plotting Multiple Simulations for Plots
library(dplyr)

#create dataframe
multisim <- data.frame(genint = numeric(),
                       group = character())

multiforwardsims <- function(iterations, vectortype, N, beta, gamma){
  
  for (i in 1:iterations){ 
      newsim = forwardsim(0,N,beta,gamma) # run a simulation from stochastic regime
      datapoints = newsim$genint
      print(datapoints)
      output = c(datapoints, rep(paste0("sim_", i), times = length(datapoints)))
      multisim <- rbind(multisim,output)
      print(multisim)
  }
}
print(multisim)