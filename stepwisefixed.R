# Forward Generation Interval Simulations with Stepwise Individual Infectiousness Profile
# Stepwise with Beta fixed implying variable R0

# Loading packages
library(ggplot2)
library(tidyverse)
library(dplyr)
library(hrbrthemes)

fswsim <- data.frame(generation = numeric(),
                       simnumb = character()
)
# Note: simulation isolated in sims.R
fixedstepwise <- function(iterations,vectortype,cohortsize,transmissionrate,recoveryrate){
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
      tau = rexp(1,recoveryrate) #exp distributed secondary infection (generation interval)
      reproductive = (transmissionrate/recoveryrate)*(transmissionrate*tau) # auc of profile
      contacts = rpois(1,reproductive) # generate number of secondary contacts
      ##print(paste0("Number of contacts: ", contacts," for individual: ", i," from initial cohort infected at time: ", vector[i]," ."))
      # If no contacts produced...
      if (contacts == 0){
        # Include individual in data frame but denote no contacts
        ##print(paste0("There were no contacts produced by individual: ", i," from initial cohort."))
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
            ##print(paste0("Contact: ", j," for individual: ", i," did NOT yield an infection. Transmission prob = ", transmission,"."))
            output = c(i,vector[i],NA,NA)
            contactdata <- rbind(contactdata, output)
            next
          } else {
            # If the contact does produce an infection...
            infection = runif(1,min = 0, max = tau)
            infecttime = vector[i] + infection #timing of secondary infection based on init cohort member timing
            ##print(paste0("Generation interval: ", infection," and timing of secondary infection: ", infecttime," for contact: ", j," from initial cohort member ", i, "."))
            output = c(i,vector[i],infecttime,infection)
            contactdata <- rbind(contactdata, output)
          }
        }
      }
    }
    # Data Wrangling Stuff...
    names(contactdata) <- c('personnumber','initcohorttime', 'secondinfecttime','genint') # insert names in contactdata
    contactdata <- contactdata[order(contactdata$genint, decreasing = FALSE),] # order based on timing of generation interval length 
    gidata = contactdata$genint # extract the generation interval data
    if (length(gidata) == 0){
      print(paste0("NOTE: Simulation number ", k," did not produce any secondary infections!"))
    } else {
      namedintervals = cbind(gidata,paste0("sim_", k)) # append simulation group
      fswsim <- rbind(fswsim,namedintervals)}
  }
  names(fswsim) <- c('generation', 'simnumb') # insert names in multisim
  fswsim <- fswsim[!is.na(fswsim$generation),] #remove NA from no contact or other
  fswsim <- fswsim %>%
    gather(key="simnumb", value="generation") %>%
    mutate(simnumb = as.character(simnumb)) %>%
    mutate(generation = as.numeric(generation))
  
  # Add Cumulative Sum by Simulation
  fswsim$count <- ave(fswsim$generation,fswsim$simnumb, FUN=seq_along) # each time receives number of infection it was by sim
  # Normalize the Cumulative Sum
  fswsim <- fswsim %>%
    group_by(simnumb) %>%
    mutate(normcsum = count/max(count, na.rm=TRUE)) # normalize the values of individuals infected at every time per sim
  
  # Specify the Theoretical Distribution Characteristics
  maxgen = max(fswsim$generation)
  x = seq(0,maxgen,length.out = 1000) #length(multisim$generation) instead?
  # PDF
  exppdf = recoveryrate*exp(-(recoveryrate)*x)
  theoretical <- cbind(x,exppdf)
  # CDF
  expcdf = 1 - exp(-(recoveryrate)*x)
  theoretical <- cbind(theoretical,expcdf)
  
  # Plot the histogram of each simulation with density and theoretical distribution curves
  p1 <- fswsim %>%
    mutate(simnumb = fct_reorder(simnumb, generation)) %>%
    ggplot(aes(x=generation)) +
    geom_histogram(aes(y=after_stat(density), color = simnumb),alpha=0.1, binwidth = 0.25) + 
    geom_density() + #adjust value? ex. adjust = 1.5
    geom_line(data =theoretical, aes(x = x, y=exppdf), linetype = "dashed") +
    theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Generation Interval") +
    ylab("Number of Secondary Infections") +
    ggtitle(paste0("FswPDF FGIs with Beta = ", transmissionrate, " Gamma = ", recoveryrate)) +
    facet_wrap(~simnumb)
  
  #Plot the CDF of each simulation and theoretical
  p2 <- fswsim %>%
    ggplot(aes(x = generation)) +
    geom_line(aes(y = normcsum, color = "red")) +
    geom_line(data = theoretical, aes(x = x, y=expcdf), linetype = "dashed") +
    theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    ) +
    xlab("Generation Interval") +
    ylab("Cumulative Secondary Infections") +
    ggtitle(paste0("FswCDF FGIs with Beta = ", transmissionrate, " Gamma = ", recoveryrate)) +
    facet_wrap(~simnumb)
  
  print(p2) # eventually make output type as part of function call
}


