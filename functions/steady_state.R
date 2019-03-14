steady_state<- function(dose,ka,ke,vd,doses,dosinginterval,drugname,mw,IC50,m,R0){
  
SteadyState <- function(t,parameters){
  with(as.list(parameters),
       { pk <- ((ka*dose)/(vd*(ka-ke))) * ((exp(-ke*t%%dosinginterval))/(1-exp(-ke*dosinginterval)) -
                                             (exp(-ka*t%%dosinginterval))/(1-exp(-ka*dosinginterval)) )}
  )}

drug_conc <- NA
pk_times <-NA
names <- NA
DrugName<-NA
pd_prof <-NA


names <- paste0(drugname)
parameters <- c(ka = ka, dose = dose,vd = vd,ke = ke,dosinginterval = dosinginterval)
times      <- seq(0, dosinginterval*doses, by = 0.1)

concs <- SteadyState(t= times,parameters = parameters)
concs <- concs*1000/mw
pharmdyn <- R0/(1 + (concs/IC50)^m)

names <- rep(names,times = length(concs))
pk_times <- c(pk_times, times)
DrugName <- c(DrugName, names)
drug_conc <- c(drug_conc, concs)
pd_prof <- c(pd_prof, pharmdyn)
dataframe <- data.frame(pk_times,drug_conc,pd_prof,DrugName)
}