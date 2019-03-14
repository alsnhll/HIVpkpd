rm(list = ls())
# setwd('~/Documents/R/HIVpkpd') #set working directory
library(plotly)
library(ggvis)
library(deSolve)
library(dplyr)
library(rootSolve)
library(visNetwork)
source('functions/getMutMatrix.R')

drugdata <- read.table("data/drug_params.txt",header=TRUE,row.names=1)
drugdata$ke <- log(2)/drugdata$halflife
mparams = read.table("data/mut_params.txt",header=TRUE)


# create sample inputs:
doses <- 4
drugs = "ETV"
R0 = 10
model = "Multiple Dose"
units = "µM"
show_muts = "All Known Mutants"
window_units = "days"
backmut = FALSE
directmulti = FALSE
directmultmult = FALSE
selfmut = FALSE
maxfit = 0.95
s=0.9
rho = 2
sigma = 0
# window_units = "concs"

input = data.frame(doses,drugs,model,R0,window_units,units,show_muts,backmut,directmulti,directmultmult,selfmut,maxfit,s,rho,sigma)
# generate data from inputs
drugValues <- drugdata[input$drugs,]
d <- drugdata["ETV",]
m_all <- mparams[mparams$DrugName=='ETV',]





# Paste whatever code you'd like to test below:




 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Generate Q matrix~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import all mutant names
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










#~~~~~~~~~~~~~~~~~~ Mutation Network Visualization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create colors for nodes based on infectivity at Cmax
conc_ind <- matrix(0,nrow = nrow(m_all)+1,ncol = 1)
k = 2
for(i in 1:length(mut_names)){
  m <- m_all[i,]
  #adjust missing or unrealistic data
  m$s <- if(is.nan(m$s))1/input$s else m$s
  m$s <- if(m$s <=1) 1/input$maxfit else m$s
  m$rho <- if(is.nan(m$rho))input$rho else m$rho
  m$sigma <- if(is.nan(m$sigma)) input$sigma else m$sigma
  
  conc_ind[k,] = d$IC50*m$rho*((input$R0/m$s-1))^(1/(d$m*(1+m$sigma))) # calculate concentration at which R0_mut = 1
  k = k+1
}
conc_ind[1,] <- d$IC50*((input$R0-1))^(1/(d$m)) # add WT
rownames(conc_ind) = c('WT',mut_names)

#Create vector of infectivity including newly added mutants (using assumed parameters):
conc_ind_all = matrix(0,nrow = nrow(Q),ncol = 1)
rownames(conc_ind_all) = rownames(Q)
for(i in 1:length(conc_ind)){
  index_muts <- which(rownames(conc_ind_all) %in% rownames(conc_ind)[i]) # find index of known mutant 
  conc_ind_all[index_muts,] <- conc_ind[i] # assign its fitness
}

# add assumed fitness for intermediate mutants based on input values
s <- if(1/input$s <=1) 1/input$maxfit else 1/input$s
rho <- input$rho
sigma <- input$sigma
conc_ind_all[-which(rownames(conc_ind_all)%in%rownames(conc_ind))] <- d$IC50*rho*((input$R0/s-1))^(1/(d$m*(1+sigma)))

conc_ind_all <- if(input$units=="µM") conc_ind_all else if(input$units=='µg/ml') conc_ind_all*d$mw/(10^3) else conc_ind_all*d$mw

#Create a function to generate a continuous color palette from white to red
bwPal <- colorRampPalette(c('blue','white','red')) # bwPal(x) creates a list of x colors that range from white to red (i.e. bwPal(2) gives #FFFFF and #FF0000 which are white and red)

#Use function to define colors of known mutant strains
color <- as.matrix(bwPal(nrow(Q))[as.numeric(factor(conc_ind_all))]) # assigns colors for mutants based on fitness 
rownames(color) = rownames(conc_ind_all)


# create NODES data frame
id <- c(1:nrow(Q))
label <- rownames(Q)

nodes <- data.frame(id,label,color = color,title = paste0("Minimum inhibitory <br>","concentration (",input$units,"): <br>",round(conc_ind_all,3)), stringsAsFactors = FALSE)





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

















  #Create sorted color dataframe to use in legend
  fitness_sort <- sort(conc_ind_all,decreasing = T,index.return = TRUE)
fitness_adj <- as.matrix(fitness_sort$x)
rownames(fitness_adj) <- rownames(conc_ind_all)[fitness_sort$ix]
temp_names <- rownames(conc_ind_all)[-c(1,which(rownames(conc_ind_all) %in% m_all$MutName))] #identifies names of mutants added during Q matrix formation (e.g. intermediates)
rownames(fitness_adj)[which(rownames(fitness_adj) %in% temp_names)] = 'Added Mutants'



# # Function to plot color bar
# 
# lut <- colorRampPalette(c("blue", "red"))(100)
# min <- round(min(conc_ind_all),digits = 2)
# max <- round(max(conc_ind_all),digits =2)
# nticks = 10
# title = paste("Viral fitness (R0) at Cmax")
# ticks=round(seq(min, max, len=nticks),digits = 3)
# scale = (length(lut)-1)/(max-min)
# 
# max <- if(max <=1)1.5 else max # make sure R0=1 is included even if all are below
# 
# # dev.new(width=1.75, height=5)
# par(mar=c(4,7,2,1)) 
# plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
# axis(2, c(1,fitness_adj), las=1,labels = c('R0 = 1',rep('',length(fitness_adj))))
# for (i in 1:(length(lut)-1)) {
#   y = (i-1)/scale + min
#   rect(0,y,10,y+1/scale, col=lut[i], border=NA)
# }

