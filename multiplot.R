# Plotting Multiple Simulations for Plots
library(dplyr)

#create dataframe
multisim <- data.frame(genint = numeric(),
                       group = character())
output=c()

multiforwardsims <- function(iterations, vectortype, N, beta, gamma){
  
  for (i in 1:iterations){ 
      newsim = forwardsim(0,N,beta,gamma) # run a simulation from stochastic regime
      output <- cbind(newsim$genint,paste0("sim_", i))
      print(output)
      #multisim <- rbind(multisim,output)
      #print(multisim)
  }
}
#print(multisim)