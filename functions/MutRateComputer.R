MutRateComputer <- function(this_mut){
  
this_mut1 <- unlist(strsplit(this_mut,''))

#load data
aamap <- read.table("data/revgeneticcode.txt",header = TRUE,row.names = 1) #map from amino acid letter to codon
mut_matrix_DNA <- read.table("data/mut_matrix_DNA.txt",header = TRUE,row.names = 1) #matrix table of mutation rates from each starting to final nucleotide for HIV
base_compo_DNA <- read.table("data/base_compo_DNA.txt",header = TRUE) #table with base composition frequencies of HIV genome

# find all starting codons
first_codons <- t(aamap[,this_mut1[1]])
first_codons <- as.vector(first_codons[!is.na(first_codons)])

first_codon <- matrix(NA,nrow = length(first_codons),ncol = 3)
for(i in 1:length(first_codons)){
  first_codon[i,] <- unlist(strsplit(first_codons[i],''))
}
first_codon <- noquote(first_codon)

#find all ending codons
last_codons <- t(aamap[,this_mut1[length(this_mut1)]])
last_codons <- as.vector(last_codons[!is.na(last_codons)])

last_codon <- matrix(NA,nrow = length(last_codons),ncol = 3)
for(i in 1:length(last_codons)){
  last_codon[i,] <- unlist(strsplit(last_codons[i],''))
}
last_codon <- noquote(last_codon)





single_nt_subs <- 0 #switch - will become 1 if a single nucleotide substitution is possible

weight_each <- matrix(0,nrow = dim(first_codon)[1],ncol = 1) #frequency of each codon in HIV genome based on bp frequency
mut_val_per_codon <- matrix(0,nrow = dim(first_codon)[1],ncol = 1) # stores sum of rates at which each original codon can become one of the possible mutant codons
comp_codons <- matrix(0,nrow = 2, ncol = 3)

for(g in 1:dim(first_codon)[1]){ #for each codon that codes for original AA
  
  #get frequency of each starting codon by multiplying frequency of each bp in codon
  
  weight_each[g] <- base_compo_DNA[,first_codon[g,1]]*base_compo_DNA[,first_codon[g,2]]*base_compo_DNA[,first_codon[g,3]]
  
  mut_val = matrix(0,nrow = dim(last_codon)[1],1)
  
  for(g2 in 1:dim(last_codon)[1]){ # for each codon that codes for mutant AA
    
      #find nucleotide differences between mutant and origianl 
    comp_codons[1,] <- first_codon[g,]
    comp_codons[2,] <- last_codon[g2,]
    
    sub1 <- length(unique(comp_codons[,1]))-1
    sub2 <- length(unique(comp_codons[,2]))-1
    sub3 <- length(unique(comp_codons[,3]))-1
    subs <- sub1 +sub2+sub3
    
    subs_ind <- which(c(sub1,sub2,sub3) >0)# index the position of nucleotide changes
    
    mut_rate <- 1
    
    for(i in 1:length(subs_ind)){
      mut_rate <- mut_rate*mut_matrix_DNA[first_codon[g,subs_ind[i]],last_codon[g2,subs_ind[i]]]
      }
    #calculate the rate at which this particular nucleotide change happens
    mut_val[g2] <- mut_rate
    
    if(subs ==1){ #means codon has only 1 nt difference
      single_nt_subs = 1
    }
  }
  
  #record the total mutation rate for this starting codon
  mut_val_per_codon[g] = sum(mut_val)
}


#calculate total rate of mutation from original to mutant amino acid, by weighing the rate from each starting codon
#by the prevalence of that starting codon in HIV genome (based on nucleotide composition values)

weight_norm <- weight_each/sum(weight_each) #get relative prevalence of each starting codon
epsilon <- sum(weight_norm*mut_val_per_codon)


epsilon


}