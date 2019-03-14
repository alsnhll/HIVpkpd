library(plotly)
library(deSolve)
library(dplyr)
library(ggplot2)
library(rootSolve)
library(DT)
require(visNetwork)
# Load functions to save space:
source('functions/MissingMutationFinder.R')
source('functions/MutRateComputer.R')
source('functions/getMutMatrix.R')
source('functions/multi_dose.R')
source('functions/steady_state.R')
source('functions/single_mut_plot.R')

# Load data
drugdata <- read.table("data/drug_params.txt",header=TRUE,row.names=1)
drugdata$ke <- log(2)/drugdata$halflife # define ke constant from known relationship to terminal halflife
mparams = read.table("data/mut_params.txt",header=TRUE)


# Define server logic required 
shinyServer(function(input, output,session) {
  
  
  
  #get the drug values for selected drug
  drugValues <- reactive({
    
    drugdata[input$drugs,]
    
  })
  

   # Create adjusted table for known PK/PD values output
  formattedDrugValues <- reactive({
    
    d<-drugdata[input$drugs,]
    # choose variables to include in table
    d=d[c("drugname","dose","mw","drugclass","IC50", "m","cmax","tmax","halflife","dosinginterval")]
    
    # express cmax and IC50 in correct units
    d$cmax <- if(input$units=="ng/ml")d$cmax*d$mw else if (input$units=="µg/ml") d$cmax*d$mw/1000 else d$cmax
    d$IC50 <- if(input$units=="ng/ml")d$IC50*d$mw else if (input$units=="µg/ml") d$IC50*d$mw/1000 else d$IC50
    
    #rename
    names(d)[names(d)=="halflife"] <- "Half-life (h)"
    names(d)[names(d)=="dose"] <- "Dose (mg)"
    names(d)[names(d)=="mw"] <- "Molecular weight (g/mol)"
    names(d)[names(d)=="drugname"] <- "Drug"
    names(d)[names(d)=="dosinginterval"] <- "Dosing interval (h)"
    names(d)[names(d)=="m"] <- "Hill slope (m)"
    names(d)[names(d)=="IC50"] <- paste("IC50 (",input$units,")",sep="")
    names(d)[names(d)=="cmax"] <-paste("Cmax (",input$units,")",sep = "")
    names(d)[names(d)=="drugclass"] <- "Drug class"
    names(d)[names(d)=="tmax"] <- "Time to Cmax (h)"
    d
    
  }) 
  
  
  
  # Show the values using an HTML table
  output$drugvalues <-renderTable(
    formattedDrugValues(),hover = T,bordered = T,striped = F
  )
  
  

  
  # Create table of fitted params from simple PK model
  formattedFittedParams <- reactive({
    d <- drugdata[input$drugs,]
    
    d = d[c("drugname","ka","ke","vd")]
    names(d)[names(d)=="drugname"] <- "Drug"
    names(d)[names(d)=="ka"] <- "Absorption rate (1/h)"
    names(d)[names(d)=="ke"] <- "Elimination rate (1/h)"
    names(d)[names(d)=="vd"] <- "Volume of distribution (L)"
    d
  })
  
  # output adjusted table
  output$fittedparams <- renderTable(
    formattedFittedParams(),hover = T,bordered = T,striped = F
  )
  
  
  
  # Get the mutation values based on chosen drug (to be used for calculations)
  mutValues <- reactive({ 
    spec_mparams=mparams[mparams$DrugName==input$drugs,]
    rownames(spec_mparams)=spec_mparams$MutName 
    spec_mparams$DrugName=NULL
    spec_mparams
  })
  
  
  # define mutant parameters to be used in table
  formattedMutValues <- reactive({ 
    m=mutValues()
    
    # convert from s to cost 'c' (r0_mut = r0_wt/s --> r0_mut = r0_wt*c)
    m$s <- 1/m$s 

    names(m)[names(m)=="MutName"] <- "Mutant strain"
    names(m)[names(m)=="s"] <- "Relative infectivity"
    names(m)[names(m)=="rho"] <- "Fold-change in IC50"
    names(m)[names(m)=="sigma"] <- "Fractional change in slope"

    m
    
  })
  
  # Show the values using an HTML table
  output$mutvalues <- renderTable({
    formattedMutValues()
  },hover = T,bordered = T,striped = F) 
  
  
  
  #_______________________________________________GENERATE PK and PD DATA__________________________________________#
  
  
  # Generate PK/PD data to use for PK plot and PK/PD plot
  dataframe <- reactive({
    
    d<-drugdata[input$drugs,] # import drug data
    
    
    if(input$model == "No Drug"){  # if intial condition is no drug, run the multi_dose model
      dataframe <- multi_dose(d$dose,d$ka,d$ke,d$vd,input$doses,d$dosinginterval,d$drugname,d$mw,d$IC50,d$m,input$R0)
    }
    else{ # if steady state is assumed
      dataframe <- steady_state(d$dose,d$ka,d$ke,d$vd,input$doses,d$dosinginterval,d$drugname,d$mw,d$IC50,d$m,input$R0)
    }
    
      pk_times <- dataframe$pk_times[!is.na(dataframe$pk_times)]/24
      drug_conc <- dataframe$drug_conc[!is.na(dataframe$drug_conc)]
      DrugName <- dataframe$DrugName[!is.na(dataframe$DrugName)]
      pd_prof <- dataframe$pd_prof[!is.na(dataframe$pd_prof)]
      dataframe <- data.frame(pk_times,drug_conc,pd_prof,DrugName)
      
    
    dataframe
    
  })
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pharmacokinetic Profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
    output$plotPK <- renderPlotly({
    d <- drugdata[input$drugs,]
    pk_times <- dataframe()$pk_times
    
    #adjust concentrations based on units input
    drug_conc <- if(input$units=="µM") dataframe()$drug_conc else if(input$units=="µg/ml")dataframe()$drug_conc*d$mw/(10^3) else dataframe()$drug_conc*d$mw
    
    DrugName <- dataframe()$DrugName
    
    dataframe <- data.frame(pk_times,drug_conc,DrugName) #store relevant info in dataframe
    
    
    #generate plot output:
    p=plot_ly(data = dataframe,x=pk_times,y=drug_conc,type='scatter',mode='lines', line=list(width=4),width = 1100, height = 550,
              name=input$drugs)
    
    # set tick marks to 8 hours if dosing interval is less than 1 day
    if(d$dosinginterval < 24){
      tval <- seq(0,max(pk_times),by=0.25)
      ttxt <- rep("",length(tval))  
      ttxt[seq(1,length(ttxt),by=2)] <- as.character(tval)[seq(1,length(ttxt),2)] 
      
    }else{
    tval <- seq(0,max(pk_times),by=0.25)
    ttxt <- rep("",length(tval))  
    ttxt[seq(1,length(ttxt),by=4)] <- as.character(tval)[seq(1,length(ttxt),4)] }
    
    
    p=layout(p,xaxis=list(zeroline=F, showline=T, ticks = "outside",tickvals = tval,ticktext = ttxt,title='Time (days)'),yaxis=list(zeroline=F, showline=T,title=paste0("Drug Concentration (", input$units,")"),range = c(0,max(drug_conc)*1.15)),hovermode="closest",
             title=paste("Plasma pharmacokinetics with",input$drugs,"treatment"),barmode='stack',
             margin=list(l=100,t=82,b=100), autosize = F)
    
    
    
  })
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PK/PD Profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
    output$plotpd <- renderPlotly({
    
    
    d<-drugdata[input$drugs,]

    pd_times <-dataframe()$pk_times
    vals <- dataframe()$pd_prof
    drugnames <-as.character(dataframe()$DrugName)
    pk_vals <- dataframe()$drug_conc
    
    pddata <- data.frame(pd_times,vals,drugnames)
    
    times      <- seq(0, input$doses *d$dosinginterval, by = 0.1)/24
    
    # set tick marks (8hrs if dosing interval < 24 hours)
    if(d$dosinginterval < 24){
      tval <- seq(0,max(pd_times),by=0.25)
      ttxt <- rep("",length(tval))  
      ttxt[seq(1,length(ttxt),by=2)] <- as.character(tval)[seq(1,length(ttxt),2)] 
      
    }else{
      tval <- seq(0,max(pd_times),by=0.25)
      ttxt <- rep("",length(tval))  
      ttxt[seq(1,length(ttxt),by=4)] <- as.character(tval)[seq(1,length(ttxt),4)] }
    
    # adjust y-values depending on UI input
    if(input$fitness_disp=='r0'){
      wt_vals <- vals # using intantaneous R0
    }else{
      wt_vals <- log10(1 + (pk_vals/d$IC50)^(d$m)) #using IIP
    }
    
    
    # generate general plot of WT in presence of drug
    p=plot_ly(data = pddata,x=pd_times,y=wt_vals,type='scatter',mode='lines', line=list(width=4),width = 1100, height = 550,
              name=paste0("Wild Type"))
    if(input$fitness_disp == 'r0'){
    p=layout(p,xaxis=list(zeroline=F, showline=T, ticks = "outside",tickvals = tval,ticktext = ttxt,title='Time (days)'),yaxis=list(zeroline=F, showline=T,title=HTML(paste("Instantaneous viral fitness (R", tags$sub(0),')', sep = "")),type = "log"),hovermode="closest",
             title=paste("PK/PD profile with", input$drugs,"treatment"),barmode='stack',
             margin=list(l=120,t=82,b=90), autosize = F)
    }else{
      p=layout(p,xaxis=list(zeroline=F, showline=T, ticks = "outside",tickvals = tval,ticktext = ttxt,title='Time (days)'),yaxis=list(zeroline=F, 
                showline=T,title="Instantaneous inhibitory potential (IIP)"),hovermode="closest",
               title=paste("PK/PD profile with", input$drugs,"treatment"),barmode='stack',
               margin=list(l=120,t=82,b=90), autosize = F)
      
    }
    
    # ADD MUTANT STRAINS:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (input$show_muts=="All Known Mutants"){ #include all mutants
      
      m_all=mutValues()
      
      R0_mut=matrix(0,nrow(m_all),length(times))
      i=1
      for (j in 1:nrow(m_all)){
        m=m_all[j,]
        if(input$disp_nas==TRUE){ # display mutants with incomplete data
            # replace missing values
        m$s <- if(is.nan(m$s))1/input$s else m$s
        m$s <- if(m$s <= 1) 1/input$maxfit else m$s 
        m$rho <- if(is.nan(m$rho))input$rho else m$rho
        m$sigma <- if(is.nan(m$sigma)) input$sigma else m$sigma
        if(input$fitness_disp == 'r0'){ #display y values as R0
          R0_mut[i,]=input$R0/(m$s)/(1+(pk_vals/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
        }else{ #display y values as IIP
          R0_mut[i,] = log10(1 + (pk_vals/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
        }
        p=add_trace(p,x=times,y=R0_mut[i,],name=row.names(m))
        }else{
          #exclude mutants with missing data
          if(is.na(m$s)|is.na(m$rho)|is.na(m$sigma)) next
          
          m$s <- if(m$s <=1) 1/input$maxfit else m$s #fix fitness if unrealistic
          
          if(input$fitness_disp == 'r0'){ #display y values as R0
            R0_mut[i,]=input$R0/(m$s)/(1+(pk_vals/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
          }else{ #display y values as IIP
            R0_mut[i,] = log10(1 + (pk_vals/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
          }
          p=add_trace(p,x=times,y=R0_mut[i,],name=row.names(m))
        }
        i=i+1 # use i as index to avoid misnaming variables if skipped
      }
    }
    else if(input$show_muts=="Test Mutant"){ # test mutant selected
      # mutant params = inputs:
      rho <- input$rho
      s <- if(input$s >=1) 1/input$maxfit else 1/input$s # fix fitness if unrealistic
      sigma <- input$sigma
      if(input$fitness_disp=='r0'){ #display y values as R0
        R0_mut=input$R0/(s)/(1+(pk_vals/(d$IC50*rho))^(d$m*(1+sigma)))
      }else{ #display y values as IIP
        R0_mut = log10(1 + (pk_vals/(d$IC50*rho))^(d$m*(1+sigma)))
      }
      p=add_trace(p,x=times,y=R0_mut,name=paste0("Test Mutant"))
      
    }
    else if(input$show_muts=="Single Known Mutant"){ # choose a single mutant strain to analyze
      m_all=mutValues()
        index <- match(input$mut_select,m_all$MutName) #find mutant values
        s <- m_all$s[index]
        rho <-m_all$rho[index]
        sigma <-m_all$sigma[index]
        # fix values if missing or unrealistic:
        s <- if(is.na(s))1/input$s else s
        s <- if(s <= 1) 1/input$maxfit else s
        rho <- if(is.na(rho))input$rho else rho
        sigma <- if(is.na(sigma))input$sigma else sigma
        
        if(input$fitness_disp=='r0'){
          R0_mut=input$R0/(s)/(1+(pk_vals/(d$IC50*rho))^(d$m*(1+sigma)))
        }else{
          R0_mut = log10(1 + (pk_vals/(d$IC50*rho))^(d$m*(1+sigma)))
        }
        p=add_trace(p,x=times,y=R0_mut,name=input$mut_select)
    }
    if(input$fitness_disp=='r0'){
    p=add_trace(p,x=times,y=rep(1,length(times)),showlegend=FALSE,
                line=list(color=rgb(200, 200, 200,max=255),dash=7),text='R_0=1, efficacy threshold')
    }
    p
    
  })
  
  


  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pharmacodynamics Profile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$drPlot <- renderPlotly({
    
    d=drugValues()
    
    #define concentration range
    conc_range=10^seq(from=input$conc_range[1],to=input$conc_range[2],by=0.025)
    
    #adjust units - IC50s etc all in molar units (uM), but can adjust to grams instead  
    plot_conc_range <- if(input$units=="µM")conc_range else if(input$units=="µg/ml")conc_range*d$mw/(10^3)else conc_range*d$mw 
    plot_Cmax <-  if (input$units == "ng/ml") d$cmax*d$mw else if (input$units == "µg/ml") d$cmax*d$mw/1000 else d$cmax
    
    #get fitness values for WT
    if(input$fitness_disp == 'iip'){
      R0_range = log10(1 + (conc_range/d$IC50)^d$m ) # define fitness by IIP
    }else{  
      R0_range=input$R0/(1+(conc_range/d$IC50)^d$m) # define fitness as R0 
      }
    
    #put in data frame
    drdata=data.frame(plot_conc_range,R0_range)
    
    # plot WT first
    p=plot_ly(data=drdata,x=plot_conc_range,y=R0_range,type='scatter',mode='lines', line=list(width=4),width = 1100, height = 550,
              name="Wild Type")
    
    if(input$fitness_disp == 'r0'){
    p=layout(p,xaxis=list(type='log',title=paste0("Drug Concentration (", input$units,")")),
             yaxis=list(title=HTML(paste("Viral fitness (R", tags$sub(0),')', sep = ""))),hovermode="closest",
             title=paste("Dose-response curve with", input$drugs,"treatment"),margin=list(l=120,t=110,b=80), autosize = F)
    }else{
      p=layout(p,xaxis=list(type='log',title=paste0("Drug Concentration (", input$units,")")),
               yaxis=list(title="Instantaneous inhibitory potential (IIP)"),hovermode="closest",
               title=paste("Dose-response curve with", input$drugs,"treatment"),margin=list(l=120,t=110,b=80), autosize = F)
    }
    
    # ADD MUTANT STRAINS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (input$show_muts=="All Known Mutants"){ #all mutants added
      
      m_all=mutValues()
      mut_names=as.character(m_all$MutName)
      
      R0_range_mut=matrix(0,nrow(m_all),length(conc_range))
      i=1
      for (j in 1:nrow(m_all)){
        m=m_all[j,]
        if(input$disp_nas==TRUE){
          # adjust unrealistic or missing mutant values
        m$s <- if(is.na(m$s))1/input$s else m$s
        m$s <- if(m$s <= 1) 1/input$maxfit else m$s
        m$rho <- if(is.na(m$rho))input$rho else m$rho
        m$sigma <- if(is.na(m$sigma)) input$sigma else m$sigma
        
        if(input$fitness_disp == 'iip'){
          # fitness as IIP
          R0_range_mut[i,] = log10(1 + (conc_range/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
        }else{ #fitness as R0
          R0_range_mut[i,]=input$R0/(m$s)/(1+(conc_range/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
          }
        
        # add to plot
        p=add_trace(p,x=plot_conc_range,y=R0_range_mut[i,],name=row.names(m))
        
        }else{ #exclude missing data
          if(is.na(m$s)|is.na(m$rho)|is.na(m$sigma)) next
          
          m$s <- if(m$s <= 1)1/input$maxfit else m$s
          
          if(input$fitness_disp == 'iip'){
            R0_range_mut[i,] = log10(1 + (conc_range/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
          }else{R0_range_mut[i,]=input$R0/(m$s)/(1+(conc_range/(d$IC50*m$rho))^(d$m*(1+m$sigma)))}
          
          p=add_trace(p,x=plot_conc_range,y=R0_range_mut[i,],name=row.names(m))
          
        }
        i=i+1
        }
    }
    else if(input$show_muts=="Test Mutant"){ # use inputs for test mutant strain
      rho <- input$rho
      s <- if(input$s < 1 ) 1/input$s else 1/input$maxfit
      sigma <- input$sigma
      
      if(input$fitness_disp == 'iip'){
        R0_range_mut = log10(1 + (conc_range/(d$IC50*rho))^(d$m*(1+sigma)))
      }else{R0_range_mut=input$R0/(s)/(1+(conc_range/(d$IC50*rho))^(d$m*(1+sigma)))}
      
      p=add_trace(p,x=plot_conc_range,y=R0_range_mut,name=paste0("Test Mutant"))
      
    }else if(input$show_muts=="Single Known Mutant"){ # choose single known mutant
      m_all=mutValues()
      index <- match(input$mut_select,m_all$MutName)
      s <- m_all$s[index]
      rho <-m_all$rho[index]
      sigma <-m_all$sigma[index]
      
      #adjust unrealistic or missing parameters
      s <- if(is.na(s))1/input$s else s
      s <- if(s <= 1) 1/input$maxfit else s
      rho <- if(is.na(rho))input$rho else rho
      sigma <- if(is.na(sigma))input$sigma else sigma
      
      
      if(input$fitness_disp == 'iip'){
        R0_range_mut = log10(1 + (conc_range/(d$IC50*rho))^(d$m*(1+sigma)))
      }else{R0_range_mut=input$R0/(s)/(1+(conc_range/(d$IC50*rho))^(d$m*(1+sigma)))}
      
      p=add_trace(p,x=plot_conc_range,y=R0_range_mut,name=input$mut_select)
      
    }
    

    #add a trace at Cmax
    p=add_trace(p,x=rep(plot_Cmax,2),y=c(0,input$R0),showlegend=FALSE,
                line=list(color=rgb(150, 150, 150,max=255),dash=3),text='Cmax for daily dosing')
    
    if(input$fitness_disp == 'r0'){
    #add trace at R_0 = 1
    p=add_trace(p,x=plot_conc_range,y=rep(1,length(plot_conc_range)),showlegend=FALSE,
                line=list(color=rgb(200, 200, 200,max=255),dash=7),text='R_0=1, efficacy threshold')
    }
    
    p
    
  })
  
  
  
  
  
  
  
  #--------------------------------------------------------- Mutant selection Window Plot:-------------------------------------
  
  output$mswPlot <- renderPlotly({
    
    d=drugValues()
    m_all=mutValues()
    max_iter=10^5
    mut_names=as.character(m_all$MutName)
    #find drug level at which wild type is suppressed (R0=1)
    Cwt=d$IC50*((input$R0-1)^(1/d$m))
    twt=max(0,d$halflife*log(d$cmax/Cwt)/log(2))/24
    
    tmax=ceiling(twt*3/5)*5
    cmax= d$cmax
    
    if(input$show_muts == "All Known Mutants"){
      
      #for each mutant
      Cint_top = rep(0,nrow(m_all))
      Cint_bot = rep(0,nrow(m_all))
      C_intersect = rep(0,nrow(m_all))
      Cmut = rep(0,nrow(m_all))
      
      tbot_int = rep(0,nrow(m_all))
      ttop_int = rep(0,nrow(m_all))
      tmut = rep(0,nrow(m_all))
      plot_type = rep(0,nrow(m_all))
      skips = NA
      i=1 #use for indexing when skipping (to keep strain names correctly aligned)
      for (j in 1:nrow(m_all)){
        
        m=m_all[j,]
        
        if(input$disp_nas == TRUE){
          #adjust missing or unrealistic data
           m$s <- if(is.nan(m$s))1/input$s else m$s
           m$s <- if(m$s <=1) 1/input$maxfit else m$s
          m$rho <- if(is.nan(m$rho))input$rho else m$rho
          m$sigma <- if(is.nan(m$sigma)) input$sigma else m$sigma
        
          } else{
            # exclude mutants with missing data
            if(is.na(m$s)){skips = c(skips,as.character(mut_names[j]))}
            if(is.na(m$s))next # skip if NA
            if(is.na(m$rho)){skips = c(skips,as.character(mut_names[j]))}
            if(is.na(m$rho))next # skip if NA
            if(is.na(m$sigma)){skips = c(skips,as.character(mut_names[j]))}
            if(is.na(m$sigma))next # skip if NA
            m$s <- if(m$s <=1) 1/input$maxfit else m$s # Replace if mutant baseline fitness is greater than WT
      }
          
        
        #find drug level at which each strain is suppressed (R0_mut=1)
        
        Cmut[i]= d$IC50*m$rho*((input$R0/m$s-1))^(1/(d$m*(1+m$sigma)))
        Cmut[i] <- if(is.nan(Cmut[i])) 10^-10 else Cmut[i] # Fixes error if baseline mutant fitness < 1. Also store this number
        
        if (Cmut[i]>Cwt){ # case 1: MSW definitely exists
          
          plot_type[i] <- 1
          
          #find drug level at which each strain is equal to wild type
          fun = function(x) {
            input$R0/(1+(x/d$IC50)^d$m)-
              (input$R0/m$s)/(1+(x/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
          }
          
          root_struct=uniroot(fun,c(0,Cmut[i]),maxiter=max_iter)
          
          C_intersect[i]= if (root_struct$iter<max_iter) root_struct$root else 0
          
          ttop_int[i]= max(0,d$halflife*log(d$cmax/Cmut[i])/log(2))/24 # Time to reach R0 = 1 for mutant
          tbot_int[i]= min(max(0,d$halflife*log(d$cmax/C_intersect[i])/log(2))/24,tmax) # Time to reach intersect with WT
         i=i+1 
        } else { # case 2: Possibly no MSW - must test to confirm:
          # Test to see if strain is ever equal to wild type
          fun = function(x) {
            input$R0/(1+(x/d$IC50)^d$m)-
              (input$R0/m$s)/(1+(x/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
          }
          
          root_struct=uniroot.all(fun,c(0,Cwt),maxiter=max_iter)
          
          Cint_bot[i] = root_struct[1]
          Cint_top[i] = root_struct[2]
          
          if(is.na(Cint_bot[i])){ # Confirmed no MSW:
            plot_type[i] <- 2
            i=i+1

            }else{ # delayed MSW
            plot_type[i] <- 3
            
            tmut[i]= max(0,d$halflife*log(d$cmax/Cmut[i])/log(2))/24 # time to reach R0=1 for mutant
            tbot_int[i]= min(max(0,d$halflife*log(d$cmax/Cint_bot[i])/log(2))/24,tmax) # Time to reach lower intersection
            ttop_int[i]= min(max(0,d$halflife*log(d$cmax/Cint_top[i])/log(2))/24,tmax) # Time to reach upper intersection
            
            i=i+1
          }
        }

}
      
      
      skips = skips[!is.na(skips)]
      
      mut_names=mut_names[!mut_names %in% skips]
      
      

  if(input$window_units == "days"){ #display results in terms of time since last dose
    l <- list(
      font = list(
        family = "sans-serif",
        size = 12,
        color = "#000"),
      orientation = 'v',
      bgcolor = "white",
      bordercolor = "white",
      traceorder = "normal")

    tval <- seq(0,tmax,by=0.25)
    ttxt <- rep("",length(tval))
    ttxt[seq(1,length(ttxt),by=4)] <- as.character(tval)[seq(1,length(ttxt),4)]


    p = plot_ly(y="Wild Type",x = twt,type = 'bar',orientation = 'h',width = 1100, height = 550,
                  marker = list(color = rgb(162, 208, 121,max=255)),
                  name = paste('','Suppressive therapy','',sep = '<br>'))
    p=add_trace(p,y="Wild Type",x=0,type='bar',orientation='h',
                marker=list(color=rgb(127, 20, 22,max=255)),
                name=paste('','Only MUT can grow','',sep='<br>'))
    p = add_trace(p,y ="Wild Type",x = 0,type = 'bar',orientation = 'h',
                  marker = list(color = 'indianred'),
                  name = paste('','Only WT can grow', '',sep='<br>'))
    p = add_trace(p,y = "Wild Type",x = 0,
                  marker=list(color=rgb(237, 28, 36,max=255)),
                  name=paste('','MUT + WT growth', 'MUT selected',sep='<br>'))
    p = add_trace(p,y = "Wild Type",x = tmax-twt,
                  marker = list(color=rgb(202, 216, 224,max=255)),
                  name = paste('','MUT + WT growth', 'WT selected',sep='<br>'))
    p=layout(p,xaxis=list(zeroline=F,  ticks = "outside",tickvals = tval,ticktext = ttxt,title='Time since last dose (days)'),yaxis=list(title=''),hovermode="closest",
             title=paste("Selection regimes with", input$drugs,"treatment"),barmode='stack',
           margin=list(l=120,t=82),showlegend=TRUE, legend = l, autosize = F)


    for (i in 1:length(mut_names)){

      # Plot Outcomes:
      if(plot_type[i] == 1){ #case 1 - MSW
        p=add_trace(p,y=mut_names[i],x=ttop_int[i],type='bar',orientation='h',
                    marker = list(color = rgb(162, 208, 121,max=255)),showlegend=FALSE,
                    name='Suppressive therapy')
        p=add_trace(p,y=mut_names[i],x=twt-ttop_int[i],type='bar',orientation='h',
                    marker=list(color=rgb(127, 20, 22,max=255)),showlegend=FALSE,
                    name=paste('','Only MUT can grow','',sep='<br>'))
        p=add_trace(p,y=mut_names[i],x=tbot_int[i]-twt,type='bar',orientation='h',
                    marker=list(color=rgb(237, 28, 36,max=255)),showlegend=FALSE,
                    name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
        p=add_trace(p,y=mut_names[i],x=tmax-tbot_int[i],type='bar',orientation='h',
                    marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                    name=paste('','MUT + WT grow', 'WT selected',sep='<br>'))

        p}else if(plot_type[i] == 2){ # case 2 - no MSW

          p=add_trace(p,y=mut_names[i],x=twt,type='bar',orientation='h',
                      marker = list(color = rgb(162, 208, 121,max=255)),showlegend=FALSE,
                      name='Suppressive therapy')
          p=add_trace(p,y=mut_names[i],x=tmax-twt,type='bar',orientation='h',
                      marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                      name=paste('WT Growth Selected',sep='<br>'))
          p

        }else if(plot_type[i]==3){ # case 3 - delayed MSW

          p = add_trace(p,y=mut_names[i],x = twt,type = 'bar',orientation = 'h',
                        marker = list(color = rgb(162, 208, 121,max=255)),showlegend=FALSE,
                        name = 'Suppressive therapy')
          p = add_trace(p,y = mut_names[i],x = tmut[i]-twt,type = 'bar',orientation = 'h',
                        marker = list(color = 'indianred'),showlegend=FALSE,
                        name = paste('','Only WT can grow', '',sep='<br>'))
          p = add_trace(p,y = mut_names[i],x = ttop_int[i]-tmut[i],type = 'bar',orientation = 'h',
                        marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                        name = paste('','MUT + WT grow', 'WT selected',sep='<br>'))
          p = add_trace(p,y = mut_names[i],x = tbot_int[i] - ttop_int[i],type = 'bar',orientation = 'h',
                        marker=list(color=rgb(237, 28, 36,max=255)),showlegend=FALSE,
                        name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
          p = add_trace(p,y = mut_names[i],x = tmax-tbot_int[i],type = 'bar',orientation = 'h',
                        marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                        name = paste('','MUT + WT grow', 'WT selected',sep='<br>'))

          p


        }}

    }else{ # If Selection Window Input is Concentrations:

              C_intersect <- if(input$units == "µM")C_intersect else if(input$units == "µg/ml")C_intersect*d$mw/1000 else C_intersect*d$mw
              Cint_bot <- if(input$units == "µM")Cint_bot else if(input$units == "µg/ml")Cint_bot*d$mw/1000 else Cint_bot*d$mw
              Cint_top <- if(input$units == "µM")Cint_top else if(input$units == "µg/ml")Cint_top*d$mw/1000 else Cint_top*d$mw
              Cwt <- if(input$units == "µM")Cwt else if(input$units == "µg/ml")Cwt*d$mw/1000 else Cwt*d$mw
              Cmut <- if(input$units == "µM")Cmut else if(input$units == "µg/ml")Cmut*d$mw/1000 else Cmut*d$mw
              cmax <- if(input$units == "µM")d$cmax else if(input$units == "µg/ml")d$cmax*d$mw/1000 else d$cmax*d$mw

              l <- list(
                font = list(
                  family = "sans-serif",
                  size = 12,
                  color = "#000"),
                orientation = 'v',
                bgcolor = "white",
                bordercolor = "white",
                traceorder = "reversed"
                )

              # Plot Wild Type (to create plotly)

              p=plot_ly(y="Wild Type",x=Cwt,type='bar',width = 1100, height = 550,
                        marker = list(color=rgb(202, 216, 224,max=255)),
                        name=paste('','MUT + WT growth', 'WT selected',sep='<br>'))
              p = add_trace(p,y = "Wild Type",x = 0,
                            marker=list(color=rgb(237, 28, 36,max=255)),
                            name=paste('','MUT + WT growth', 'MUT selected',sep='<br>'))
              p = add_trace(p,y ="Wild Type",x = 0,type = 'bar',
                            marker = list(color = 'indianred'),
                            name = paste('','Only WT can grow', '',sep='<br>'))
              p=add_trace(p,y="Wild Type",x=0,type='bar',
                          marker=list(color=rgb(127, 20, 22,max=255)),
                          name=paste('','Only MUT can grow','',sep='<br>'))
              p=add_trace(p,y="Wild Type",x=cmax-Cwt,type='bar',
                          marker = list(color = rgb(162, 208, 121,max=255)),
                          name=paste('','Suppressive therapy','',sep='<br>'))
              p=layout(p,xaxis=list(range = c(0,cmax),title=paste0("Drug Concentration (", input$units,")")),yaxis=list(title=''),hovermode="closest",
                       title=paste("Selection regimes with", input$drugs,"treatment"),barmode='stack',
                       margin=list(l=120,t=82),showlegend = TRUE,legend = l, autosize = F)


              for (i in 1:length(mut_names)){

                if(plot_type[i] == 1){ # case 1 - MSW

                    p=add_trace(p,y=mut_names[i],x=C_intersect[i],
                                marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                                name=paste('','MUT + WT grow', 'WT selected',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=Cwt-C_intersect[i],
                                marker=list(color=rgb(237, 28, 36,max=255)),showlegend=FALSE,
                                name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=Cmut[i]-Cwt,
                                marker=list(color=rgb(127, 20, 22,max=255)),showlegend=FALSE,
                                name=paste('','Only MUT can grow','',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=cmax-Cmut[i],
                                marker = list(color = rgb(162, 208, 121,max=255)),showlegend=FALSE,
                                name='Suppressive therapy')

                    p


                  }else{if(plot_type[i] ==2){ # case 2 - no MSW


                    p=add_trace(p,y=mut_names[i],x=Cwt,
                                marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                                name=paste('','MUT + WT grow', 'WT selected',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=cmax-Cwt,
                                marker = list(color = rgb(162, 208, 121,max=255)),showlegend=FALSE,
                                name='Suppressive therapy')

                    p


                  }else{ # case 3 - delayed MSW

                    p=add_trace(p,y=mut_names[i],x=Cint_bot[i],
                                marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                                name=paste('','MUT + WT grow', 'WT selected',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=Cint_top[i]-Cint_bot[i],
                                marker=list(color=rgb(237, 28, 36,max=255)),showlegend=FALSE,
                                name=paste('','MUT + WT grow', 'MUT selected',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=Cmut[i]-Cint_top[i],
                                marker = list(color=rgb(202, 216, 224,max=255)),showlegend=FALSE,
                                name = paste('','MUT + WT grow', 'WT selected',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=Cwt-Cmut[i],
                                marker = list(color = 'indianred'),showlegend=FALSE,
                                name = paste('','Only WT can grow', '',sep='<br>'))
                    p=add_trace(p,y=mut_names[i],x=cmax-Cwt,
                                marker = list(color = rgb(162, 208, 121,max=255)),showlegend=FALSE,
                                name='Suppressive therapy')

                    p



                  }}}
              p
        }
      }
    
    # REPEAT DATA AND PLOT GENERATING FOR TEST MUTANT AND SINGLE KNOWN MUTANT
    
    
    #---------------------------------------------------------Test Mutant MSW Plot--------------------------------------------------------------
    
    
    
    else if(input$show_muts=="Test Mutant"){
      m_all <- mutValues()
      # If "Show Mutant Strains" isn't selected, show the Test Mutant Selection Window
     
      s <- 1/input$s
      s <- if(input$s>=1) 1/input$maxfit else s
      rho <- input$rho
      sigma <- input$sigma
      m <- data.frame(s,rho,sigma)
      
      
      # Find drug level at which WT is suppressed:
      Cwt=d$IC50*((input$R0-1)^(1/d$m))
      twt=max(0,d$halflife*log(d$cmax/Cwt)/log(2))/24 # Time (since reaching Cmax) to R0 = 1 for WT strain
      
      tmax=ceiling(twt*3/5)*5
      
      #find drug level at which test strain is suppressed (R0_mut=1)
      Cmut= d$IC50*m$rho*((input$R0/m$s-1))^(1/(d$m*(1+m$sigma)))
      Cmut <- if(is.nan(Cmut)) 10^-10 else Cmut # Fixes error if baseline mutant fitness < 1
      
      if (Cmut>Cwt){ #MSW definitely exists
        plot_type <- 1
        
        # Find drug level at which test strain is equal to wild type
        fun = function(x) {
          input$R0/(1+(x/d$IC50)^d$m)-
            (input$R0/m$s)/(1+(x/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
        }
        
        root_struct=uniroot(fun,c(0,Cmut),maxiter=max_iter)
        
        C_intersect= if (root_struct$iter<max_iter) root_struct$root else 0
        
        ttop_mut= max(0,d$halflife*log(d$cmax/Cmut)/log(2))/24 # Time to reach R0 = 1 for mutant
        tbot_mut= min(max(0,d$halflife*log(d$cmax/C_intersect)/log(2))/24,tmax) # Time to reach intersect of two strains
        
        
      } else { # Possibly no MSW - must test to confirm:
        
        # Test to see if strain is equal to wild type
        fun = function(x) {
          input$R0/(1+(x/d$IC50)^d$m)-
            (input$R0/m$s)/(1+(x/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
        }
        
        root_struct=uniroot.all(fun,c(0,Cmut),maxiter=max_iter)
        
        Cint_bot = root_struct[1]
        Cint_top = root_struct[2]
        
        if(is.na(Cint_bot)){ # Confirmed no MSW:
          plot_type <- 2
          
        }
        else{ # delayed MSW
          plot_type <- 3
          
          tmut= max(0,d$halflife*log(d$cmax/Cmut)/log(2))/24 # time to reach R0=1 for mutant
          tbot_int= min(max(0,d$halflife*log(d$cmax/Cint_bot)/log(2))/24,tmax) # Time to reach lower intersection
          ttop_int= min(max(0,d$halflife*log(d$cmax/Cint_top)/log(2))/24,tmax) # Time to reach upper intersection
          
          
        }
      }
      
      # since it's only one trace being added, define data frame based on MSW plot type
      
      if(input$window_units == "days"){ # plot as time since last dose
        if(plot_type == 1){  #type 1 - MSW
          df <- data.frame(twt,ttop_mut,tbot_mut,tmax)
        }else{
          if(plot_type ==2){ # type 2 - no MSW
            df <- data.frame(twt,tmax)
          }else{ # type 3 - delayed MSW
            df <- data.frame(twt,tmut,ttop_int,tbot_int,tmax)
          }
        }
      }else{ # plot as concentration ranges
        if(plot_type ==1){ #type 1
          df <- data.frame(Cwt,C_intersect,Cmut,cmax)
        }else{
          if(plot_type ==2){ # type 2
            df <- data.frame(Cwt,cmax)
          }else{ # type 3
            df <- data.frame(Cwt,Cint_bot,Cint_top,Cmut,cmax)
          }
        }
      }
      
      p <- single_mut_plot(plot_type,input$window_units,d,Cwt,twt,tmax,df,input$units,"Test Mutant") # call plot function for single mutant strain
      
      
      #------------------------------------------- Single Known Mutant MSW ----------------------------------------------------
      # same as test mutant, but with predefined parameters
      
      
    }else if(input$show_muts=="Single Known Mutant"){
      m_all <- mutValues()

        index <- match(input$mut_select,m_all$MutName)
        s <- m_all$s[index]
        rho <-m_all$rho[index]
        sigma <-m_all$sigma[index]
        
        s <- if(is.na(s))1/input$s else s
        s <- if(s<=1) 1/input$maxfit else s
        
        rho <- if(is.na(rho))input$rho else rho
        sigma <- if(is.na(sigma))input$sigma else sigma
        
      
        
        m <-data.frame(s,rho,sigma)

      # Find drug level at which WT is suppressed:
      Cwt=d$IC50*((input$R0-1)^(1/d$m))
      twt=max(0,d$halflife*log(d$cmax/Cwt)/log(2))/24 # Time (since reaching Cmax) to R0 = 1 for WT strain
      
      tmax=ceiling(twt*3/5)*5
      
      #find drug level at which test strain is suppressed (R0_mut=1)
      Cmut= d$IC50*m$rho*((input$R0/m$s-1))^(1/(d$m*(1+m$sigma)))
      Cmut <- if(is.na(Cmut)) 10^-10 else Cmut # Fixes error if baseline mutant fitness < 1
      
      
      if (Cmut>Cwt){ #MSW definitely exists
        plot_type <- 1
        
        # Find drug level at which test strain is equal to wild type
        fun = function(x) {
          input$R0/(1+(x/d$IC50)^d$m)-
            (input$R0/m$s)/(1+(x/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
        }
        
        root_struct=uniroot(fun,c(0,Cmut),maxiter=max_iter)
        
        C_intersect= if (root_struct$iter<max_iter) root_struct$root else 0
        
        ttop_mut= max(0,d$halflife*log(d$cmax/Cmut)/log(2))/24 # Time to reach R0 = 1 for mutant
        tbot_mut= min(max(0,d$halflife*log(d$cmax/C_intersect)/log(2))/24,tmax) # Time to reach intersect of two strains
        
        
      } else { # Possibly no MSW - must test to confirm:
        
        # Test to see if strain is equal to wild type
        fun = function(x) {
          input$R0/(1+(x/d$IC50)^d$m)-
            (input$R0/m$s)/(1+(x/(d$IC50*m$rho))^(d$m*(1+m$sigma)))
        }
        
        root_struct=uniroot.all(fun,c(0,Cmut),maxiter=max_iter)
        
        Cint_bot = root_struct[1]
        Cint_top = root_struct[2]
        
        if(is.na(Cint_bot)){ # Confirmed no MSW:
          plot_type <- 2
          
        }
        else{ # delayed MSW
          plot_type <- 3
          
          tmut= max(0,d$halflife*log(d$cmax/Cmut)/log(2))/24 # time to reach R0=1 for mutant
          tbot_int= min(max(0,d$halflife*log(d$cmax/Cint_bot)/log(2))/24,tmax) # Time to reach lower intersection
          ttop_int= min(max(0,d$halflife*log(d$cmax/Cint_top)/log(2))/24,tmax) # Time to reach upper intersection
          
          
        }
      }
      
      # define data set based on type of plot
      if(input$window_units == "days"){ # plot as time since last dose
        if(plot_type == 1){ # MSW exists
          df <- data.frame(twt,ttop_mut,tbot_mut,tmax)
        }else{
          if(plot_type ==2){ # no MSW
            df <- data.frame(twt,tmax)
          }else{ #delayed MSW
            df <- data.frame(twt,tmut,ttop_int,tbot_int,tmax)
          }
        }
      }else{ # plot as conc ranges
        if(plot_type ==1){ # MSW exists
          df <- data.frame(Cwt,C_intersect,Cmut,cmax)
        }else{
          if(plot_type ==2){ # no MSW
            df <- data.frame(Cwt,cmax)
          }else{ # delayed MSW
            df <- data.frame(Cwt,Cint_bot,Cint_top,Cmut,cmax)
          }
        }
      }      
        
      # call plotting function for single mutant strain
      p <- single_mut_plot(plot_type,input$window_units,d,Cwt,twt,tmax,df,input$units,input$mut_select)

}


p

})
  
  
  
  
  
    
    
    
    #-------------------------------------------  generate Q Matrix and network visualization  ----------------------------------------------------
    
    
    
    

    # Generate Q matrix for use in outputs
    mut_matrix <- reactive({
      
      #import all mutant names
      m_all=mutValues()
      mut_names=as.character(m_all$MutName)
      
      mut_names_all <- NA
      # reformat names to enter into getMutMatrix (replace '/' with '_')
      for(i in 1:length(mut_names)){
        split_output <- t(unlist(strsplit(mut_names[i],'/')))
        adj_mutname <- paste(split_output,collapse = '_')
        mut_names_all <- c(mut_names_all,adj_mutname)
      }
      mut_names_all <- mut_names_all[!is.na(mut_names_all)]
      
      # define params for Q matrix generation
      back_mutation_on <-0
      direct_multi_hit <-0
      direct_multi_multi <-0
      if(input$backmut == TRUE){
        back_mutation_on <-1
      }
      if(input$directmulti ==TRUE){
        direct_multi_hit <-1
      }
      if(input$directmultmult==TRUE){
        direct_multi_multi <-1
      }
      
      #run getMutMatrix function with reactive inputs
      Q <- getMutMatrix(mut_names_all,back_mutation_on,direct_multi_hit,direct_multi_multi)
      
      #reformat Q col/row names to replace '_' with '/'
      mutnames <- colnames(Q)
      fixed_mut_names <- NA
      for(i in 1:length(mutnames)){
        split_output <- t(unlist(strsplit(mutnames[i],'_')))
        adj_mutname <- paste(split_output,collapse = '/')
        fixed_mut_names <- c(fixed_mut_names,adj_mutname)
      }
      fixed_mut_names <- fixed_mut_names[!is.na(fixed_mut_names)]
      colnames(Q) <- fixed_mut_names
      rownames(Q) <- fixed_mut_names
      Q
    })
    
    
    
    # generate warning if someone chooses direct multi multi without back mutations
    output$matwarning <- renderText({ 
      error <- NA
      if(input$directmultmult == TRUE){
        if(input$backmut == FALSE){
          error <- paste("Must have back mutations enabled!") 
        }
      }
      error
    })
    
    
    
  
    

    # Create Mutation Matrix Table Output

  output$Qmatrix <- renderTable({
    Q <- mut_matrix()

    Qdata <-Q[,1] 
    for(i in 2:ncol(Q)){
        Qdata <- cbind(Qdata,Q[,i])
    }
    rownames(Qdata) <- rownames(Q)
    colnames(Qdata) <- colnames(Q)
    
    # eliminate 0's for easier viewing (otherwise they come up as exponentials and expand the table)
    for(i in 1:ncol(Q)){
      for(j in 1:nrow(Q)){
        if(Qdata[i,j] ==0){
          Qdata[i,j] = NA
        }
      }
      if(input$selfmut == FALSE){
        Qdata[i,i] = NA
      }
    }
    Qdata <- as.data.frame(Qdata)
    Qdata
  },rownames = T,striped = F,bordered = T,hover = T,digits = -3,na = 0)
      

    
    
    
    
    
  #-------------------------------------------  Mutation Network Visualization  ----------------------------------------------------
  
    output$network <- renderVisNetwork({
      
      #import relevant data
      m_all=mutValues()
      d<-drugdata[input$drugs,]
      
      mut_names=as.character(m_all$MutName)
      Q <- mut_matrix()
      
      
      
      # create colors for nodes based on infectivity at Cmax
      mut_MIC <- matrix(0,nrow = nrow(m_all)+1,ncol = 1)
      k = 2 # index over k to avoid WT (first index)
      for(i in 1:length(mut_names)){
        m <- m_all[i,]
        #adjust missing or unrealistic data
        m$s <- if(is.nan(m$s))1/input$s else m$s
        m$s <- if(m$s <=1) 1/input$maxfit else m$s
        m$rho <- if(is.nan(m$rho))input$rho else m$rho
        m$sigma <- if(is.nan(m$sigma)) input$sigma else m$sigma
        
        mut_MIC[k,] = d$IC50*m$rho*((input$R0/m$s-1))^(1/(d$m*(1+m$sigma))) # calculate concentration at which R0_mut = 1
        
        # mut_MIC[k,] = (input$R0/m$s)/(1 + (d$cmax/(d$IC50*m$rho))^(d$m*(1+m$sigma))) # calculate R0 at cmax
        k = k+1
      }
      mut_MIC[1,] <- d$IC50*((input$R0-1))^(1/(d$m)) # add WT
      # mut_MIC[1,] = (input$R0)/(1 + (d$cmax/(d$IC50))^(d$m)) # calculate R0 at cmax for WT
      
      rownames(mut_MIC) = c('WT',mut_names) #name rows appropriately
      
      #Create vector of infectivity including newly added mutants (using assumed parameters):
      mut_MIC_all = matrix(0,nrow = nrow(Q),ncol = 1)
      rownames(mut_MIC_all) = rownames(Q)
      for(i in 1:length(mut_MIC)){
        index_muts <- which(rownames(mut_MIC_all) %in% rownames(mut_MIC)[i]) # find index of known mutant 
        mut_MIC_all[index_muts,] <- mut_MIC[i] # assign its fitness
      }
      
      # add assumed fitness for intermediate mutants based on input values
      s <- if(1/input$s <=1) 1/input$maxfit else 1/input$s
      rho <- input$rho
      sigma <- input$sigma
      mut_MIC_all[-which(rownames(mut_MIC_all)%in%rownames(mut_MIC))] <- d$IC50*rho*((input$R0/s-1))^(1/(d$m*(1+sigma)))
      # mut_MIC_all[-which(rownames(mut_MIC_all)%in%rownames(mut_MIC))] <- (input$R0/s)/(1 + (d$cmax/(d$IC50*rho))^(d$m*(1+sigma))) # calculate R0 at cmax
      
      mut_MIC_all <- if(input$units=="µM") mut_MIC_all else if(input$units=='µg/ml') mut_MIC_all*d$mw/(10^3) else mut_MIC_all*d$mw
      
      #Create a function to generate a continuous color palette from white to red
      bwPal <- colorRampPalette(c('blue','white','red')) # bwPal(x) creates a list of x colors that range from blue to red (separated by white to avoid purple middle ground)
      
      #Use function to define colors of known mutant strains
      color <- as.matrix(bwPal(nrow(Q))[as.numeric(factor(mut_MIC_all))]) # assigns colors for mutants based on fitness 
      rownames(color) = rownames(mut_MIC_all)
      
      
      # create NODES data frame
      id <- c(1:nrow(Q))
      label <- rownames(Q)
      
      nodes <- data.frame(id,label,color = color,title = paste0("Minimum inhibitory <br>","concentration (",input$units,"): <br>",round(mut_MIC_all,3)), stringsAsFactors = FALSE)
      
      
      
      # create EDGES  data frame
      from <- NA
      to <- NA
      label <- NA
      for(i in 1:length(rownames(Q))){
        for(j in 1:length(colnames(Q))){
          if(Q[i,j]!=0){
            if(input$selfmut == TRUE){
              from <- c(from,i)
              to <- c(to,j)
              label <- c(label, Q[i,j])
            }else
              if(i != j){
                from <- c(from,i)
                to <- c(to,j)
                label <- c(label, Q[i,j])
              }
          }
        }
      }
      from <- from[!is.na(from)]
      to <- to[!is.na(to)]
      label <- label[!is.na(label)]
      
      
      # edges <- data.frame(from,to,label) #includes edge lables
      edges <- data.frame(from,to)
      
      
      
      visNetwork(nodes, edges, main = paste("Resistant mutant formation in response to ",input$drugs,sep = ''), margin=list(l=100,t=80,b=100)) %>%
        
        visEdges(shadow = TRUE,
                 color = '#000000',
                 arrows =list(to = list(enabled = TRUE))) %>%
        
        
        visOptions(highlightNearest = TRUE) %>%
        visLayout(randomSeed = 123) %>%
        visInteraction(dragNodes = T,
                       dragView = F,
                       zoomView = F,
                       tooltipDelay = 0)
      
      # visLegend(width = 0.1, position = "right", main = HTML(paste("Viral fitness (R", tags$sub(0),') ','at C',tags$sub('max'), sep = "")), addNodes = col_legend,useGroups = FALSE)
      
      
    })
    

    # Text explaining colors
  output$visText <- renderText({
    text <- paste("Color scale used:","Blue: lower minimum inhibitory concentration","Red: higher minimum inhibitory concentration","* Color scale is relative to other strains. Hover over strains with cursor to view MIC *",sep = '\n')
  })
  
  
  
  # Warning for assuming parameters
  output$param_warning <- renderText({
    #import relevant data
    m_all=mutValues()
    mut_names=as.character(m_all$MutName)
    Q <- mut_matrix()
    
    
    new_muts <-rownames(Q)[-c(1,which(rownames(Q)%in% mut_names))] # store list of newly introduced mutants
    text <- if(length(new_muts) !=0)paste("Assuming parameters for the following mutants: ",paste(new_muts,collapse = ', '),'.',sep = '')
    text

  })
  
  
  
    
    # Mutation Matrix warnings (need two copies because they must be printed separately to panel and vis network)
    output$matwarning <- renderText({ 
      
      if(input$directmultmult == TRUE){
        if(input$backmut == FALSE){
          error <- paste("Error: cannot have direct mutations without back mutations enabled.") 
        }
      }
      
    })
    
   output$viswarning <- renderText({
     
     if(input$directmultmult == TRUE){
       if(input$backmut == FALSE){
         error <- paste("Error: cannot have direct mutations without back mutations enabled.") 
       }
     }
   })
    
    
    
    
    
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   #-------------------------------------------  reactive outputs for UI  (errors / warning messages) ----------------------------------------------------
   
   
   # Record instances when assuming new fitness cost:
   skips = reactive({
     
     m_all=mutValues()
     
     mut_names=as.character(m_all$MutName)
     skips <- NA
     if(input$show_muts == "All Known Mutants"){
       for (i in 1:nrow(m_all)){
         m=m_all[i,]
         if(input$disp_nas==TRUE){
           m$s <- if(is.na(m$s)) 1/input$s  else m$s
           if(m$s<=1){skips <- c(skips, mut_names[i])}
         }else
           if(is.na(m$s))next
         
         if(m$s<=1){skips <- c(skips, mut_names[i])}
       }
     }else
       if(input$show_muts== "Single Known Mutant"){
         index <- match(input$mut_select,m_all$MutName)
         s <- m_all$s[index]
         
         s <- if(is.na(s)) 1/input$s else s
         
         if(s<=1){skips <- c(skips, mut_names[index])}
         
       }
     
     skipped <- skips[!is.na(skips)]
     skipped = unique(skipped)
   })
   
   
   # Record if baseline viral fitness is too high
   fitskips = reactive({
     
     d=drugValues()
     fitskips <- NA
     
     Cwt=d$IC50*((input$R0-1)^(1/d$m)) 
     fitskips <- if(Cwt>d$cmax)1 else 0
     
     
   })
   
   
   
   
   
   # Allow user to select individual known mutants for a drug if "Single Known Mutant" is selected
   outVar = reactive({
     drugname = input$drugs
     m_all=mutValues()
     skips = NA
     if(input$disp_nas==TRUE){
       names = as.character(m_all$MutName)}
     else{
       if(input$show_muts == "Single Known Mutant"){
         
         m_all = m_all[complete.cases(m_all),]
       }
       names = as.character(m_all$MutName)
       
       
     }
     names
   })
   observe({
     updateSelectInput(session,"mut_select",
                       choices = outVar()
     )})
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Create Error Messages if replacing Values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # PD Plot
   output$warning <- renderText({ 
     m_all=mutValues()
     mut_names=as.character(m_all$MutName)
     s_errors = NA
     rho_errors = NA
     sigma_errors = NA
     missing_data = NA
     
     if(input$show_muts=="All Known Mutants"){
       if(input$disp_nas==TRUE){
         for (i in 1:nrow(m_all)){
           m=m_all[i,]
           s_errors = if(is.nan(m$s))c(s_errors,as.character(mut_names[i]))else s_errors
           rho_errors = if(is.nan(m$rho))c(rho_errors,as.character(mut_names[i]))else rho_errors
           sigma_errors = if(is.nan(m$sigma))c(sigma_errors,as.character(mut_names[i]))else sigma_errors}}
       else
         for (i in 1:nrow(m_all)){
           m=m_all[i,]
           if(is.na(m$s)|is.na(m$rho)|is.na(m$sigma)){missing_data = c(missing_data, as.character(mut_names[i]))} 
         }
       missing_data = missing_data[!is.na(missing_data)]
       miss = missing_data
       miss = unique(miss)
       miss = miss[!is.na(miss)]
       len_skipped <- length(miss)
       
       skipped = paste(miss,collapse=", ")
     }else
       if(input$show_muts=="Single Known Mutant"){
         len_skipped = 0
         index <- match(input$mut_select,m_all$MutName)
         s <- m_all$s[index]
         rho <-m_all$rho[index]
         sigma <-m_all$sigma[index]
         
         s_errors = if(is.nan(s))c(s_errors,as.character(mut_names[index])) else s_errors
         rho_errors = if(is.nan(rho))c(rho_errors,as.character(mut_names[index])) else rho_errors
         sigma_errors = if(is.nan(sigma))c(sigma_errors,as.character(mut_names[index])) else sigma_errors
         
       }else{          len_skipped = 0
       }
     
     s_errors <- s_errors[!is.na(s_errors)]
     sigma_errors <- sigma_errors[!is.na(sigma_errors)]
     rho_errors <- rho_errors[!is.na(rho_errors)]
     s<- if(length(s_errors)>0) paste("Assuming relative fitness for the following mutant(s): ",paste(s_errors,collapse = ", "),sep='\n')else NULL
     sigma<-if(length(sigma_errors)>0)paste("Assuming fractional slope change for the following mutant(s): ",paste(sigma_errors,collapse = ", "),sep='\n')else NULL
     rho<-if(length(rho_errors)>0) paste("Assuming fold-change of IC50 for the following mutant(s): ",paste(rho_errors,collapse = ", "),sep='\n')else NULL
     
     skips <-if(len_skipped>0)paste("The following mutants were removed due to incomplete data: ",skipped,sep = '\n')else NULL
     
     fit_fix <- paste(skips(),collapse = ", ")
     len_fitfix <- length(skips())
     
     too_fit <- if(len_fitfix>0)paste0("Limiting relative fitness of the following mutants: ", fit_fix)else NULL
     # Combine warnings and omissions:
     if(input$disp_nas==TRUE){paste(sigma,s,rho,too_fit,sep = '\n')}
     else paste(skips,too_fit,sep = '\n')
     
     
     
   })
   
   
   #PK/PD Plot
   output$warning2 <- renderText({
     
     
     
     m_all=mutValues()
     mut_names=as.character(m_all$MutName)
     s_errors = NA
     rho_errors = NA
     sigma_errors = NA
     missing_data = NA
     
     if(input$show_muts=="All Known Mutants"){
       if(input$disp_nas==TRUE){
         for (i in 1:nrow(m_all)){
           m=m_all[i,]
           s_errors = if(is.nan(m$s))c(s_errors,as.character(mut_names[i]))else s_errors
           rho_errors = if(is.nan(m$rho))c(rho_errors,as.character(mut_names[i]))else rho_errors
           sigma_errors = if(is.nan(m$sigma))c(sigma_errors,as.character(mut_names[i]))else sigma_errors}}
       else
         for (i in 1:nrow(m_all)){
           m=m_all[i,]
           if(is.na(m$s)|is.na(m$rho)|is.na(m$sigma)){missing_data = c(missing_data, as.character(mut_names[i]))} 
         }
       missing_data = missing_data[!is.na(missing_data)]
       miss = missing_data
       miss = unique(miss)
       miss = miss[!is.na(miss)]
       len_skipped <- length(miss)
       
       skipped = paste(miss,collapse=", ")
     }else
       if(input$show_muts=="Single Known Mutant"){
         len_skipped = 0
         index <- match(input$mut_select,m_all$MutName)
         s <- m_all$s[index]
         rho <-m_all$rho[index]
         sigma <-m_all$sigma[index]
         
         s_errors = if(is.nan(s))c(s_errors,as.character(mut_names[index])) else s_errors
         rho_errors = if(is.nan(rho))c(rho_errors,as.character(mut_names[index])) else rho_errors
         sigma_errors = if(is.nan(sigma))c(sigma_errors,as.character(mut_names[index])) else sigma_errors
         
       }else{          len_skipped = 0
       }
     
     s_errors <- s_errors[!is.na(s_errors)]
     sigma_errors <- sigma_errors[!is.na(sigma_errors)]
     rho_errors <- rho_errors[!is.na(rho_errors)]
     s<- if(length(s_errors)>0) paste("Assuming relative fitness for the following mutant(s): ",paste(s_errors,collapse = ", "),sep='\n')else NULL
     sigma<-if(length(sigma_errors)>0)paste("Assuming fractional slope change for the following mutant(s): ",paste(sigma_errors,collapse = ", "),sep='\n')else NULL
     rho<-if(length(rho_errors)>0) paste("Assuming fold-change of IC50 for the following mutant(s): ",paste(rho_errors,collapse = ", "),sep='\n')else NULL
     
     skips <-if(len_skipped>0)paste("The following mutants were removed due to incomplete data: ",skipped,sep = '\n')else NULL
     
     fit_fix <- paste(skips(),collapse = ", ")
     len_fitfix <- length(skips())
     too_fit <- if(len_fitfix>0)paste0("Limiting relative fitness of the following mutants: ", fit_fix)else NULL
     # Combine warnings and omissions:
     if(input$disp_nas==TRUE){paste(sigma,s,rho,too_fit,sep = '\n')}
     else paste(skips,too_fit,sep = '\n')
     
     
   })
   
   
   
   
   # Selection Window Plot
   output$warning3 <- renderText({
     m_all=mutValues()
     mut_names=as.character(m_all$MutName)
     s_errors = NA
     rho_errors = NA
     sigma_errors = NA
     missing_data = NA
     
     if(input$show_muts=="All Known Mutants"){
       if(input$disp_nas==TRUE){
         for (i in 1:nrow(m_all)){
           m=m_all[i,]
           s_errors = if(is.nan(m$s))c(s_errors,as.character(mut_names[i]))else s_errors
           rho_errors = if(is.nan(m$rho))c(rho_errors,as.character(mut_names[i]))else rho_errors
           sigma_errors = if(is.nan(m$sigma))c(sigma_errors,as.character(mut_names[i]))else sigma_errors}}
       else
         for (i in 1:nrow(m_all)){
           m=m_all[i,]
           if(is.na(m$s)|is.na(m$rho)|is.na(m$sigma)){missing_data = c(missing_data, as.character(mut_names[i]))} 
         }
       missing_data = missing_data[!is.na(missing_data)]
       # miss = missing_data
       # miss = unique(miss)
       # miss = miss[!is.na(miss)]
       len_skipped <- length(missing_data)
       
       skipped = paste(missing_data,collapse=", ")
     }else
       if(input$show_muts=="Single Known Mutant"){
         len_skipped = 0
         
         index <- match(input$mut_select,m_all$MutName)
         s <- m_all$s[index]
         rho <-m_all$rho[index]
         sigma <-m_all$sigma[index]
         
         s_errors = if(is.nan(s))c(s_errors,as.character(mut_names[index])) else s_errors
         rho_errors = if(is.nan(rho))c(rho_errors,as.character(mut_names[index])) else rho_errors
         sigma_errors = if(is.nan(sigma))c(sigma_errors,as.character(mut_names[index])) else sigma_errors
         
       }else{          len_skipped = 0
       }
     
     s_errors <- s_errors[!is.na(s_errors)]
     sigma_errors <- sigma_errors[!is.na(sigma_errors)]
     rho_errors <- rho_errors[!is.na(rho_errors)]
     s<- if(length(s_errors)>0) paste("Assuming relative fitness for the following mutant(s): ",paste(s_errors,collapse = ", "),sep='\n')else NULL
     sigma<-if(length(sigma_errors)>0)paste("Assuming fractional slope change for the following mutant(s): ",paste(sigma_errors,collapse = ", "),sep='\n')else NULL
     rho<-if(length(rho_errors)>0) paste("Assuming fold-change of IC50 for the following mutant(s): ",paste(rho_errors,collapse = ", "),sep='\n')else NULL
     
     
     # Omissions from Selection window
     
     
     fit_fix <- paste(skips(),collapse = ", ")
     len_fitfix <- length(skips())
     
     r0skip <- fitskips()
     r0_errors <- if(r0skip==1)paste0("The baseline WT fitness is too high to be suppressed by this drug.") else NULL 
     too_fit <- if(len_fitfix>0)paste0("Limiting relative fitness of the following mutants: ", fit_fix)else NULL
     
     omits <- if(len_skipped>0)paste0("The following mutants were omitted due to missing parameters: ", skipped)else NULL
     
     
     # Combine warnings and omissions:
     if(input$show_muts=="All Known Mutants"){
       if(input$disp_nas==TRUE){paste(sigma,s,rho,too_fit,r0_errors,sep = '\n')}else{
         paste(omits,too_fit,r0_errors,sep = '\n')
       }}else{
         if(input$show_muts=="Single Known Mutant"){
           paste(sigma,s,rho,too_fit,r0_errors,sep = '\n')
         }}
     
     
   })
   
   
   
   
   
   
   
    
    
    
    
    
    
    
    
    
    
    
    
     
})