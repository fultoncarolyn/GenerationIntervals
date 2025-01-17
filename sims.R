# Simulation of Forward Generation Interval Measurement

# GOAL: Create df with original timing of infection, timing of secondary infection, 
#       and the different between for one generation/layer

# Loading packages
library(ggplot2)

forwardsim <- function(vectortype,cohortsize,transmissionrate,recoveryrate){
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
  vector <- sort(vector) # put vector of initial cohort into chronological order ##necessary?##
  
# Create data frame this goes into
contactdata <- data.frame(personnumber = numeric(),
                          initcohorttime = numeric(),
                          secondinfecttime = numeric(),
                          genint = numeric())

# Execute the contact regime for one generation of secondary infections to be added to data frame
for (i in 1:length(vector)){
  contacts = rpois(1,transmissionrate/recoveryrate) # generate number of secondary contacts
  print(paste0("Number of contacts:", contacts,"for individual:", i," from initial cohort."))
  # If no contacts produced...
  if (contacts == 0){
    # Include individual in data frame but denote no contacts
    print(paste0("There were no contacts produced by individual:", i," from initial cohort."))
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
        print(paste0("Contact:", j,"for individual:", i,"did NOT yield an infection. Transmission prob =", transmission,"."))
        ## Does not add to dataframe to avoid clutter - would we need to track contacts that do not become infections?
        next
      } else{
        # If the contact does produce an infection...
        infecttime = rexp(1,recoveryrate) #timing of secondary infection
        print(paste0("Timing of secondary infection:", infecttime,"for contact:", j," from initial cohort member", i, "."))
        output = c(i,vector[i],vector[i] + infecttime,infecttime) 
        contactdata <- rbind(contactdata, output)
      }
    }
  }
}
names(contactdata) <- c('personnumber','initcohorttime', 'secondinfecttime','genint') # insert names in contactdata
print(contactdata)

# Generate plots

#Forward-Scheme figure (only really useful for scenario 1)
#ggplot(contactdata, aes(initcohorttime, genint)) + geom_point() + labs(x = 'Time of Initial Cohort Infection', y = 'Generation Interval', title = 'Forward Model Schematic')

par(mfrow=c(1,2))

print(contactdata)

#Histogram Plot
plot <- contactdata %>%
  ggplot( aes(x = genint)) + 
  geom_histogram(bins = 100) + 
  geom_density(color="blue") +
  geom_density(adjust = 0.5, color="green", linetype="dotted")
  labs(x = 'Generation Interval Length', y = 'Number of Secondary Infections', title = 'Forward Model Histogram')
print(plot)

#Density Plot
#ggplot(contactdata, aes(x = genint)) + geom_density() + labs(x = 'Generation Interval Length', y = 'Frequency of Secondary Infections', title = 'Forward Model Density Distribution')

}



