# Simulation of Forward Generation Interval Measurement

#GOAL: Create df with original timing of infection, timing of secondary infection, 
#       and the different between for one generation/layer

#Step One: Generate a number of infected individuals and the times at which they were infected
#           Refer to this as the initial cohort

# Parameters (for now)
beta = 1.75
gamma = 0.35

initcohort = c()
initcohort <- sort(c(0.2,0.1)) #define strict first individuals and their timing

#Step Two: Order the initial cohort chronologically

###initcohort = sort(init)

#Step Three: Loop over initial cohort to generate secondary infections

#Assumptions: All infections stick -> contact implies infected
#             The timing of infection for secondary is identical


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
    #do something here for df
    
  } else {
    for (j in 1:contacts){
      timeofinfect = rexp(1,gamma) #timing of secondary infection
      print(paste0("Timing of secondary infection:", timeofinfect,"for contact:", j," from initial cohort member", i, "."))
      output = c(i,vector[i],timeofinfect,timeofinfect - vector[i]) #gen int timing vs timeofinfect
      contactdata <- rbind(contactdata, output)
    }
  }
}
names(contactdata) <- c('personnumber','initcohorttime', 'secondinfecttime','genint') # names aren't working but not worried about right now...

print(contactdata)
}

#Step Four: Map these to a data frame for easier manipulation

