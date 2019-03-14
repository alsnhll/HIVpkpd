MissingMutationFinder <- function(this_mut){

this_mut <- unlist(strsplit(this_mut,''))


#load data
aamap <- read.table("data/revgeneticcode.txt",header = TRUE,row.names = 1) #map from amino acid letter to codon
codonmap <- read.table("data/geneticcode.txt",header = TRUE,row.names = 1) #map from codon to amino acid letter (stop codons changed from * to NA)
mut_matrix_DNA <- read.table("data/mut_matrix_DNA.txt",header = TRUE,row.names = 1) #matrix tabel of mutation rates from each starting to final nucleotide for HIV

# find all starting codons
first_codons <- t(aamap[,this_mut[1]])
first_codons <- as.vector(first_codons[!is.na(first_codons)])

first_codon <- matrix(NA,nrow = length(first_codons),ncol = 3)
for(i in 1:length(first_codons)){
  first_codon[i,] <- unlist(strsplit(first_codons[i],''))
}
first_codon <- noquote(first_codon)

#find all ending codons
last_codons <- t(aamap[,this_mut[length(this_mut)]])
last_codons <- as.vector(last_codons[!is.na(last_codons)])

last_codon <- matrix(NA,nrow = length(last_codons),ncol = 3)
for(i in 1:length(last_codons)){
  last_codon[i,] <- unlist(strsplit(last_codons[i],''))
}
last_codon <- noquote(last_codon)

#loci <- noquote(paste(this_mut[2:4]))

#if >3 difference, give all intermediate states
intermediates_2step <- NA
intermediates_aa_2step <- NA
intermediates_3step <- c(NA,NA)
intermediates_aa_3step <- c(NA,NA)
num_subs <- matrix(0,nrow = dim(first_codon)[1],ncol = dim(last_codon)[1])
comp_codons <- matrix(0,nrow = 2, ncol = 3)

for(g in 1:dim(first_codon)[1]){ # for each codon that codes for original AA
  for(g2 in 1:dim(last_codon)[1]){ # for each codon that codes for mutant AA 
    
    #find nucleotide differences between mutant and origianl 
    comp_codons[1,] <- first_codon[g,]
    comp_codons[2,] <- last_codon[g2,]

    sub1 <- length(unique(comp_codons[,1]))-1
    sub2 <- length(unique(comp_codons[,2]))-1
    sub3 <- length(unique(comp_codons[,3]))-1
    subs <- sub1 +sub2+sub3
    
    num_subs[g,g2] <- subs
    
    if(subs == 2){ # codon has 2 nt difference
      
      
      
    #find the nucleotide that changed between mutant and original
      subs_ind <- which(c(sub1,sub2,sub3) >0)# index the position of nucleotide changes
      
      int1 <- first_codon[g,]
      int1[subs_ind[1]] <- last_codon[g2,subs_ind[1]]
      int1 <- paste(int1[1],int1[2],int1[3],sep = "")
      int2 <- first_codon[g,]
      int2[subs_ind[2]] <- last_codon[g2,subs_ind[2]]
      int2 <- paste(int2[1],int2[2],int2[3],sep = "")
      
      intermediates_2step <- rbind(intermediates_2step,int1,int2)
      
      #change back to amino acid
      intermediates_aa_2step <- rbind(intermediates_aa_2step,as.character(codonmap[rownames(codonmap)==int1,]),as.character(codonmap[rownames(codonmap)==int2,]))
    }
  }
}
  

if(min(num_subs)==1){
  # print_line = print("no missing mutations")
  intermediates <- NA
}else if(min(num_subs)==2){
  
  #find unique amino acid intermediates
  intermediates_aa_2step_u <- unique(intermediates_aa_2step)
  intermediates <- intermediates_aa_2step_u[!is.na(intermediates_aa_2step_u)] # stop codons were changed to NA's
  
}else if(min(num_subs) ==3){
  print_line = print("requires 3 codon changes!")
  
  for(g in 1:dim(first_codon)[1]){ # for each codon that codes for original AA
    for(g2 in 1:dim(last_codon)[1]){ # for each codon that codes for mutant AA   
      
      int1 <- paste(last_codon[g2,1],first_codon[g,2],first_codon[g,3],sep = "")
      int2 <- paste(first_codon[g,1],last_codon[g2,2],first_codon[g,3],sep = "")
      int3 <- paste(first_codon[g,1],first_codon[g,2],last_codon[g2,3],sep = "")
      
      int4 <- paste(last_codon[g2,1],last_codon[g2,2],first_codon[g,3],sep = "")
      int5 <- paste(last_codon[g2,1],first_codon[g,2],last_codon[g2,3],sep = "")
      int6 <- paste(first_codon[g,1],last_codon[g2,2],last_codon[g2,3],sep = "")
      
    intermediates_3step <- rbind(intermediates_3step,c(int1,int4),c(int1,int5),c(int2, int4),c(int2,int6),c(int3, int5),c(int3,int6))
    
    # change back to AA
    intermediates_aa_3step <- rbind(intermediates_aa_3step,c(as.character(codonmap[rownames(codonmap)==int1,]),as.character(codonmap[rownames(codonmap)==int4,])),
                                    c(as.character(codonmap[rownames(codonmap)==int1,]),as.character(codonmap[rownames(codonmap)==int5,])),c(as.character(codonmap[rownames(codonmap)==int2,]),as.character(codonmap[rownames(codonmap)==int4,])),
                                    c(as.character(codonmap[rownames(codonmap)==int2,]),as.character(codonmap[rownames(codonmap)==int6,])),c(as.character(codonmap[rownames(codonmap)==int3,]),as.character(codonmap[rownames(codonmap)==int5,])),
                                    c(as.character(codonmap[rownames(codonmap)==int3,]),as.character(codonmap[rownames(codonmap)==int6,])))
    
    intermediates_aa_3step_u <- unique(intermediates_aa_3step)
    x <- complete.cases(intermediates_aa_3step_u)
    intermediates <- intermediates_aa_3step_u[x,]
    }
  }
}

intermediates <- as.matrix(sort(intermediates))
}