# Simulation of Forward Generation Interval Measurement

# GOAL: Create df with original timing of infection, timing of secondary infection, 
#       and the different between for one generation/layer

# Loading packages
library(ggplot2)

# Parameters (for now)
beta = 1.75
gamma = 0.35

###### Step One: Generate a number of infected individuals and the times at which they were infected ###### 

# Refer to this as the initial cohort

initcohort = c()
initcohort <- sort(c(0.2,0.1,0.3,0.4,0.5,0.6,0.7,0.8)) # for now define strict first individuals and their timing (diff vs same)

####### Step Two: Order the initial cohort chronologically (if necessary) ###### 

####### Step Three: Loop over initial cohort to generate secondary infections ######
#######             AND map to data frame for manipulation ###### 

#Assumptions: All infections stick -> contact implies infected
#             The timing of infection for secondary is NOT identical


funloop <- function(vector){

# Create data frame this goes into
contactdata <- data.frame(personnumber = numeric(),
                          initcohorttime = numeric(),
                          secondinfecttime = numeric(),
                          genint = numeric())

# Execute the contact regime to be added to data frame
for (i in 1:length(vector)){
  contacts = rpois(1,beta/gamma) # generate number of secondary contacts
  print(paste0("Number of contacts:", contacts,"for individual:", i," from initial cohort."))
  if (contacts == 0){
    #do something here for df or is next fine?
    next
  } else {
    for (j in 1:contacts){
      infecttime = rexp(1,gamma) #timing of secondary infection
      print(paste0("Timing of secondary infection:", infecttime,"for contact:", j," from initial cohort member", i, "."))
      output = c(i,vector[i],vector[i] + infecttime,infecttime) 
      contactdata <- rbind(contactdata, output)
    }
  }
}
names(contactdata) <- c('personnumber','initcohorttime', 'secondinfecttime','genint') # names aren't working but not worried about right now...

print(contactdata)

###### Step Four: Generate plots #####

#Forward Scheme
ggplot(contactdata, aes(initcohorttime, genint)) + geom_point() + labs(x = 'Time of Initial Cohort Infection', y = 'Generation Interval', title = 'Forward Model Schematic')

#Density/Distribution
ggplot(contactdata, aes(x = genint)) + geom_histogram(bins = 25) + labs(x = 'Generation Interval', y = 'Number of Secondary', title = 'Forward Model Histogram')


}



