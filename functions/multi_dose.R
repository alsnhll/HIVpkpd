multi_dose <- function(dose,ka,ke,vd,doses,dosinginterval,drugname,mw,IC50,m,R0){

## One Compartment Model
one_comp <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), { dDose<- -ka*Dose
  dComp <- (ka*Dose)/Vd - ke*Comp
  list(c(dDose, dComp))
  })
}

drug_conc <- NA
pk_times <-NA
names <- NA
DrugName<-NA
pd_prof <-NA


concs <- 0
parameters <- c(ka = ka, ke = ke, Vd = vd)
state      <- c(Dose = dose, Comp = 0)
times      <- seq(0, dosinginterval, by = 0.1)
names <- paste0(drugname)

# Run the ODE solver for selected drugs / doses
for(j in 1:doses){
  out <- ode(y = state, times = times, func = one_comp, parms = parameters) # Solve ODE
  
  output<- as.data.frame(out) #Change output type from deSolve to data frame
  concs <- c(concs,output$Comp[2:length(output$Comp)]) #Record drug concentrations in a vector
  state <- c(Dose = dose+last(output$Dose), Comp = last(output$Comp)) #Set the starting state to the end of the last state + newly administered dose
}
concs <- concs*1000/mw
pharmdyn <- R0/(1 + (concs/IC50)^m)

drug_conc <- c(drug_conc, concs)
names <- rep(names,times = length(concs))
pk_times <- c(pk_times, seq(0,doses*dosinginterval, by = 0.1))
DrugName <- c(DrugName, names)
pd_prof <- c(pd_prof, pharmdyn)


# Remove NaNs
pk_times <- pk_times[!is.na(pk_times)]
drug_conc <- drug_conc[!is.na(drug_conc)]
DrugName <- DrugName[!is.na(DrugName)]
pd_prof <- pd_prof[!is.na(pd_prof)]

dataframe <- data.frame(pk_times,drug_conc,pd_prof,DrugName)

}