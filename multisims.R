# Plotting for Multiple Forward Generation Interval Simulations

# Represents base case -> Exponential Individual Infectiousness Profile

# GOAL: Create epidemiologically motivated figures to explore 
#       forward generation interval simulations

# Loading packages
library(ggplot2)
library(tidyverse)
library(dplyr)
library(hrbrthemes)

multisim <- data.frame(generation = numeric(),
                       simnumb = character()
                       )
# Note: simulation isolated in sims.R
plotsim <- function(iterations,vectortype,cohortsize,transmissionrate,recoveryrate){
  for (k in 1:iterations){
  
  # Generate a number of infected individuals and the times at which they were infected  
  # Refer to this as the initial cohort
  # Note: transmissionrate : beta :: recoveryrate : (gamma from SIR framework)
  if (vectortype == 0){
    vector <- rep(0,cohortsize) # all initial cohort infections occur at time=0
  } else if (vectortype == 1){
    vector <- runif(cohortsize, min=0, max=1) # all initial cohort infections occur between 0 and 1 uniformally
  } else {
    print("ERROR in vectortype")
  }
  vector <- sort(vector) # put vector of initial cohort into chronological order 
  
  # Create data frame this goes into
  contactdata <- data.frame(personnumber = numeric(),
                            initcohorttime = numeric(),
                            secondinfecttime = numeric(),
                            genint = numeric())
  
  # Execute the contact regime for one generation of secondary infections to be added to data frame
  for (i in 1:length(vector)){
    contacts = rpois(1,transmissionrate/recoveryrate) # generate number of secondary contacts
    ##print(paste0("Number of contacts:", contacts,"for individual:", i," from initial cohort."))
    # If no contacts produced...
    if (contacts == 0){
      # Include individual in data frame but denote no contacts
      ##print(paste0("There were no contacts produced by individual:", i," from initial cohort."))
      output = c(i,vector[i],NA,NA)
      contactdata <- rbind(contactdata, output)
      next
    } else {
      # If nonzero contacts produced...
      for (j in 1:contacts){
        # Compute whether transmission occurs for each contact
        transmission <- rbinom(1,1,transmissionrate)
        if (transmission == 0){
          # If the contact does not produce an infection...
          ##print(paste0("Contact:", j,"for individual:", i,"did NOT yield an infection. Transmission prob =", transmission,"."))
          # Does not add to dataframe to avoid clutter - would we need to track contacts that do not become infections?
          next
        } else {
          # If the contact does produce an infection...
          infecttime = rexp(1,recoveryrate) #timing of secondary infection
          ##print(paste0("Timing of secondary infection:", infecttime,"for contact:", j," from initial cohort member", i, "."))
          output = c(i,vector[i],vector[i] + infecttime,infecttime) 
          contactdata <- rbind(contactdata, output)
        }
      }
    }
    
  }
  # Data Wrangling Stuff...
  names(contactdata) <- c('personnumber','initcohorttime', 'secondinfecttime','genint') # insert names in contactdata
  gidata = sort(contactdata$genint) # extract the generation interval data and sort it
  namedintervals = cbind(gidata,paste0("sim_", k)) # append simulation group
  multisim <- rbind(multisim,namedintervals)
  }
  names(multisim) <- c('generation', 'simnumb') # insert names in multisim
  multisim <- multisim[!is.na(multisim$generation),] #remove NA from no contact or other
  multisim <- multisim %>%
    gather(key="simnumb", value="generation") %>%
    mutate(simnumb = as.character(simnumb)) %>%
    mutate(generation = as.numeric(generation))
  
  # Add Cumulative Sum by Simulation
  multisim$count <- ave(multisim$generation,multisim$simnumb, FUN=seq_along) # each time recieves number of infection it was by sim
  # Normalize the Cumulative Sum
  multisim <- multisim %>%
    group_by(simnumb) %>%
    mutate(normcsum = count/max(count, na.rm=TRUE)) # normalize the values of indidivuals infected at every time per sim

  # # Specify the Theoretical Distribution Characteristics
  # maxgen = max(multisim$generation)
  # extramaxgen = maxgen + 3 # can factor as fit
  # x = seq(0,extramaxgen,length=length(multisim$generation))
  # x = multisim$generation
  # # PDF
  # exppdf = recoveryrate*exp(-(recoveryrate)*x)
  # theoretical <- cbind(x,exppdf)
  # # CDF
  # expcdf = 1 - exp(-(recoveryrate)*x)
  # theoretical <- cbind(theoretical,expcdf)
  
  
  # Create CDF and PDF using the timing of secondary infections which in this case is defined by the generation interval
  multisim <- multisim %>%
    group_by(simnumb) %>%
    mutate(expncdf = 1 - exp(-(recoveryrate)*generation))
  
  multisim <- multisim %>%
    group_by(simnumb) %>%
    mutate(expnpdf = recoveryrate*exp(-(recoveryrate)*generation))
  
  print(multisim)
  
  # Plot the histogram of each simulation with density and theoretical distribution curves
  p1 <- multisim %>%
    mutate(simnumb = fct_reorder(simnumb, generation)) %>%
    ggplot(aes(x=generation)) +
    geom_histogram(aes(y=after_stat(density), color = simnumb),alpha=0.1, binwidth = 0.25) + 
    geom_density() + #adjust value? ex. adjust = 1.5
    geom_line(aes(y=expnpdf), linetype = "dashed") +
    theme_ipsum() +
     theme(
       legend.position="none",
       panel.spacing = unit(0.1, "lines"),
       strip.text.x = element_text(size = 8)
     ) +
    xlab("Generation Interval") +
    ylab("Number of Secondary Infections") +
    ggtitle(paste0("ExpPDF FGIs with Beta = ", transmissionrate, " Gamma = ", recoveryrate)) +
    facet_wrap(~simnumb)
  
  #Plot the CDF of each simulation and theoretical
  p2 <- multisim %>%
    ggplot(aes(x = generation)) +
    geom_line(aes(y = normcsum, color = "red")) +
    geom_line(aes(y=expncdf), linetype = "dashed") +
    #stat_function(fun = pnorm) +
    theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Generation Interval") +
    ylab("Cumulative Secondary Infections") +
    ggtitle(paste0("ExpCDF FGIs with Beta = ", transmissionrate, " Gamma = ", recoveryrate)) +
    facet_wrap(~simnumb)
  
  print(p1) # eventually make output type as part of function call
}


