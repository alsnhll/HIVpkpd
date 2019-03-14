getMutMatrix <- function(mut_names_all,back_mutation_on,direct_multi_hit,direct_multi_multi){

# import MutRateComputer and MissingMutationFinder functions
source('functions/MissingMutationFinder.R')
source('functions/MutRateComputer.R')
library(stringr)
library(data.table)

# create function to allow R to index like Matlab (wont run if start > end in the index start:end)
lim_index <- function(start,end) {
  if(end<start) return(integer(0))
  seq(start, end)
}

all_mut_list <- NA # single and intermediate multi-hit mutations that must be added
num_aa_change <- NA # number of AA substitutions required for all single and intermediate mutations

for(q in 1:length(mut_names_all)){ # for each genotype
  split_output <- t(unlist(strsplit(mut_names_all[q],'_'))) #split indivudal mutations into cells
  
  num_hits <- length(split_output)
  
  #add single AA mutations to list
  for(i in 1:length(split_output)){
    all_mut_list <- rbind(all_mut_list,split_output[i])
    num_aa_change <- rbind(num_aa_change,1)
  }
 
  if(num_hits > 1){
    
    #add the final multi-hit mutation to the list of all mutations
    #change / to _ so mutation name can also be variable name in mutation matrix
    mutcombine <- paste(split_output,sep = "_",collapse = "_")
    all_mut_list <- rbind(all_mut_list,mutcombine)
    num_aa_change <- rbind(num_aa_change,num_hits)
    
    #need to find all intermediates if it s a multi-hit (>2 step) mutation
    #find all k-step intermediates
    for(k in 2:num_hits-1){
      if(k >= 2){
      
        int_muts <- t(combn(split_output,k)) #concatenate. will always put subsets in order of appearance in full set
        
        for(j in 1:dim(int_muts)[1]){
          intmutcombine <- if(dim(int_muts)[2] ==2)paste(int_muts[j,1],int_muts[j,2],sep = "_")else if(dim(int_muts)[2]==3) paste(int_muts[j,1],int_muts[j,2],int_muts[j,3],sep = "_")
          all_mut_list <- rbind(all_mut_list,intmutcombine)
          num_aa_change<- rbind(num_aa_change,k)
        }
      }
    }
  } 
}



#find unique mutations
all_mut_list <- all_mut_list[!is.na(all_mut_list)]
num_aa_change <- num_aa_change[!is.na(num_aa_change)]

index <- which(!duplicated(all_mut_list))
all_mut_list <- all_mut_list[which(!duplicated(all_mut_list))]
num_aa_change <-num_aa_change[index]

#sort by number of hits
ind_sort <- sort(num_aa_change,decreasing = F, index.return = T)

num_aa_change <- ind_sort$x
all_mut_list <- all_mut_list[ind_sort$ix]

#separate out all single AA changes
ind_saa <- max(which(num_aa_change %in% 1))
single_aa_list <- all_mut_list[1:ind_saa]

#sort mutli-AA changes:
num_aa_change <- num_aa_change[order(num_aa_change,all_mut_list)]
all_mut_list <- all_mut_list[order(num_aa_change,all_mut_list)]


#sort single AA mutations by position in genome 
#(multi-AA already sorted this way if entered this way)

temp_list <- matrix(NA,nrow = length(single_aa_list),ncol = 1)

for(i in 1:ind_saa){
  #remove first letter, so that sorting occurs by the number
  singleAAchange <- unlist(strsplit(single_aa_list[i],''))
  locus <- paste(singleAAchange[2:(length(singleAAchange)-1)],collapse = '')
  temp_list[i] <- noquote(paste(locus,singleAAchange[length(singleAAchange)],sep = "",collapse = ''))
}

ind_sort2 <- sort(temp_list,decreasing = F, index.return = T)
single_aa_list <- single_aa_list[ind_sort2$ix]
all_mut_list[1:ind_saa] <- single_aa_list

#create mutation rate table
Q <- matrix(0,nrow = length(all_mut_list)+1,ncol = length(all_mut_list)+1)
colnames(Q) <- c('WT',all_mut_list)
rownames(Q) <- c('WT',all_mut_list)
Q <-as.table(Q)

#Fill in mutation rates for each single AA change
#if a single AA change requires 2 bp changes, add intermediates and their rates

ind_this_mut <- 2

for(i in 1:length(single_aa_list)){ #for each single aa
  
  this_mut <- single_aa_list[i]
  
  u <- MutRateComputer(this_mut)
  intermediates <- MissingMutationFinder(this_mut)
  
  if(is.na(intermediates[1])==TRUE){ #if single bp change
   
    Q['WT',this_mut] <- u
     
    if(back_mutation_on == 1){ 
      this_mut1 <- unlist(strsplit(this_mut,''))
      first_codon <- this_mut1[1]
      last_codon <- this_mut1[length(this_mut1)]
      locus <- paste(this_mut1[2:(length(this_mut1)-1)],collapse = '')
      
      rev_mut <- paste(last_codon,locus,first_codon,sep = "")
      
      ub <- MutRateComputer(rev_mut)
      Q[this_mut,'WT'] <- ub}
      
  }else{ #required intermediate states
    
    #if one intermediate state is needed
    if(dim(intermediates)[2]==1){
      
      #for each intermediate, calculate mutation rate and add it to Q
      for(j in 1:dim(intermediates)[1]){
        this_mut1 <- unlist(strsplit(this_mut,''))
        codon_plus_locus <-paste(this_mut1[1:(length(this_mut1)-1)],collapse="")
        int_mut <- paste(codon_plus_locus,intermediates[j], sep = "")
      
        #if this variable name doesn't already exist, create it
        if(is.na(match(int_mut,colnames(Q))) == TRUE){
          #augment Q
          #add column
          int_col <- matrix(0,nrow = dim(Q)[1],ncol = 1)
          colnames(int_col) <- int_mut
          Q <- cbind(Q,int_col)
          col_names <- colnames(Q)
          Q <- Q[,c(col_names[lim_index(1,ind_this_mut)],int_mut,col_names[lim_index((ind_this_mut+1),(length(col_names)-1))])] #move to correct location
          
          #add row
          int_row <- matrix(0,nrow = 1,ncol = dim(Q)[2])
          rownames(int_row) <- int_mut
          Q <- rbind(Q,int_row)
          row_names <- rownames(Q)
          Q <- Q[c(row_names[lim_index(1,ind_this_mut)],int_mut,row_names[lim_index((ind_this_mut+1),(length(row_names)-1))]),] #move to correct location
          
          ind_this_mut <- ind_this_mut +1
        }
        
        #enter mutation rates
        u1 <- MutRateComputer(int_mut)
        Q['WT',int_mut] <- u1
        
        if(back_mutation_on ==1){
          int_mut1 <- unlist(strsplit(int_mut,''))
          locus <- paste(int_mut1[2:(length(int_mut1)-1)],collapse = '')
          first_codon <- int_mut1[1]
          last_codon <- int_mut1[length(int_mut1)]
          rev_mut1 <- paste(last_codon,locus,first_codon,sep = "")
          
          u1b <- MutRateComputer(rev_mut1)
          Q[int_mut,'WT'] <- u1b
        }
        
        int_mut1 <- unlist(strsplit(int_mut,''))
        locus <- paste(int_mut1[2:(length(int_mut1)-1)],collapse = '')
        second_change <- paste(intermediates[j],locus,this_mut1[length(int_mut1)],sep="")
        u2 <- MutRateComputer(second_change)
        Q[int_mut,this_mut] <- u2
        
        if(back_mutation_on ==1){
          int_mut2 <- unlist(strsplit(second_change,''))
          locus <- paste(int_mut2[2:(length(int_mut2)-1)],collapse = "")
          first_codon <- int_mut2[1]
          last_codon <- int_mut2[length(int_mut2)]
          rev_mut2 <- paste(last_codon,locus,first_codon,sep = "")
          
          u2b <- MutRateComputer(rev_mut2)
          Q[this_mut,int_mut] <- u2b
        }
        # ind_this_mut <- ind_this_mut +1
      }
     
      if(direct_multi_multi ==1){
       
        #fill in mutation rates directly from one intermediate to another
        for(j in 1:dim(intermediates)[1]){
          for(k in 1:dim(intermediates)[1]){
            if(k != j){
              this_mut1 <- unlist(strsplit(this_mut,''))
              int_mut1 <- paste(c(this_mut1[1:length(this_mut1)-1],intermediates[j]), collapse = "")
              int_mut2 <- paste(c(this_mut1[1:length(this_mut1)-1],intermediates[k]), collapse = "")
              direct_int_int <- paste(c(intermediates[j],this_mut1[2:(length(this_mut1)-1)],intermediates[k]),collapse = "")
              u <- MutRateComputer(direct_int_int)
              Q[int_mut1,int_mut2] <- u
              
              direct_int_int_rev <- paste(c(intermediates[k],this_mut1[2:(length(this_mut1)-1)],intermediates[j]),collapse =  "")
              ub <- MutRateComputer(direct_int_int_rev)
              Q[int_mut2,int_mut1] <- ub
            }
          }
        }
      }
    }else{
      print("Error: haven't coded case where multiple intermediates are needed yet")
    }
  }
  ind_this_mut <- ind_this_mut + 1
}

#fill in mutation rates for multi-hit mutations based on single hit rates

while(ind_this_mut <= dim(Q)[1]){ # for all multi_hit mutations
  name_mut <- colnames(Q)[ind_this_mut]
  split_output <- t(unlist(strsplit(name_mut,'_'))) #split individual mutations
  
  num_hits <- length(split_output)
  
  #find all subsets with one less mutation and enter mutation rate from them
  
  for(k in 1:num_hits){
    
    #remove one mutation
    last_mut <- split_output[k]
    
    
    # #join the others
    name_one_less <- split_output[! split_output %in% last_mut]
    name_one_less <- paste(name_one_less[1:(length(name_one_less))],sep = "_",collapse = "_") 
    
    #enter mutation rate
    u <- MutRateComputer(last_mut)
    intermediates <- MissingMutationFinder(last_mut)
    if(is.na(intermediates[1])==TRUE){
      Q[name_one_less,name_mut] <-u
    }else{
      
      #if one intermediate state is needed
      if(dim(intermediates)[2]==1){
        
        ind_temp <- ind_this_mut
  
        # for each intermediate, calculate mutation rate to/from and add it to Q
        for(j in 1:dim(intermediates)[1]){
          last_mut1 <- unlist(strsplit(last_mut,''))
          last_mut1 <- paste(last_mut1[1:(length(last_mut1)-1)],sep='',collapse = '')
  
          int_mut1 <- paste(last_mut1,intermediates[j],sep = '',collapse = '') #one loci only
          
          mutcombine <- c(split_output[lim_index(1,k-1)],int_mut1,split_output[lim_index(k+1,num_hits)])
          mutcombine <- mutcombine[!is.na(mutcombine)]
          int_mut_all <- paste(mutcombine,sep = "_",collapse = "_") 
          
          #if this variable name doesn't already exist, create it
          if(is.na(match(int_mut_all,colnames(Q))) == TRUE){
            #augment Q
            #add column
            int_col <- matrix(0,nrow = nrow(Q),ncol = 1)
            colnames(int_col) <- int_mut_all
            Q <- cbind(Q,int_col)
            col_names <- colnames(Q)
            Q <- Q[,c(col_names[lim_index(1,ind_temp)],int_mut_all,col_names[lim_index((ind_temp+1),(length(col_names)-1))])] #move to correct location
            
            #add row
            int_row <- matrix(0,nrow = 1,ncol = ncol(Q))
            rownames(int_row) <- int_mut_all
            Q <- rbind(Q,int_row)
            row_names <- rownames(Q)
            Q <- Q[c(row_names[lim_index(1,ind_temp)],int_mut_all,row_names[lim_index((ind_temp+1),(length(row_names)-1))]),]
            
            ind_temp <- ind_temp +1
          }  
          
          #enter mutation rates
          u1 <- MutRateComputer(int_mut1)
          Q[name_one_less,int_mut_all] <- u1
          
          if(back_mutation_on ==1){
            int_intmut1 <- unlist(strsplit(int_mut1,''))
            rev_mut1 <- paste(c(int_intmut1[length(int_intmut1)],int_intmut1[2:length(int_intmut1)-1],int_intmut1[1]),sep = '',collapse = '')
            u1b <- MutRateComputer(rev_mut1)
            Q[int_mut_all,name_one_less] <- u1b
          }
          last_mut2 <- unlist(strsplit(last_mut,''))
          second_change <- paste(c(intermediates[j],last_mut2[2:length(last_mut2)]),sep = "",collapse = '')
          u2 <- MutRateComputer(second_change)
          Q[int_mut_all,name_mut] <- u2
          
          if(back_mutation_on ==1){
            second_change2 <- unlist(strsplit(second_change,''))
            rev_mut2 <- paste(c(second_change2[length(second_change2)],second_change2[2:(length(second_change2)-1)],second_change2[1]),sep = "",collapse = '')
            u2b <-MutRateComputer(rev_mut2)
            Q[name_mut,int_mut_all] <- u2b
          }
          # ind_temp <- ind_temp +1
        }
      }
    }
  }
  
  ind_this_mut=ind_this_mut+1
}


#add in rates at which mutations don't occur - diagonal entries
for(i in 1:ncol(Q)){
  Q[i,i] <- 1
}

#augment Q by adding in direct rates of multi-hit mutations from WT

if(direct_multi_hit ==1){

  for(i in 2:ncol(Q)){

    if(Q['WT',i]==0){

      this_mut <- colnames(Q)[i]
      split_output <- t(unlist(strsplit(this_mut,'_'))) #split individual mutations into cells
      num_hits <- length(split_output)
      u_vec <- matrix(0,nrow = num_hits,ncol = 1)

      for(k in 1:num_hits){

        u_vec[k] <- MutRateComputer(split_output[k])

      }

      Q['WT',i] <- prod(u_vec)

      if(back_mutation_on ==1){
        ub_vec <- matrix(0,nrow = num_hits,ncol = 1) #backwards mutation

        for(k in 1:num_hits){
          this_mut <- unlist(strsplit(split_output[k],''))
          locus <- paste(this_mut[2:(length(this_mut)-1)],collapse = '')
          rev_mut <- paste(c(this_mut[length(this_mut)],locus,this_mut[1]),sep = "",collapse = '')

          ub_vec[k] <- MutRateComputer(rev_mut)

        }

        Q[i,'WT'] <- prod(ub_vec)

      }
    }
  }
}

#allow direct mutation-to-mutation transitions
#back mutation must be on to allow this

if(direct_multi_multi==1){

  if(back_mutation_on ==1){

    for(i in 2:nrow(Q)){
      for(j in 2:ncol(Q)){

        if(Q[i,j]==0 & i!=j){

          u_tot = 1
          start_gen <- rownames(Q)[i]
          end_gen <- colnames(Q)[j]
          split_start <- t(unlist(strsplit(start_gen,'_')))
          split_end <- t(unlist(strsplit(end_gen,'_')))

          #record loci (numeric) values only
          start_loci <- matrix(0,nrow = length(split_start),ncol = 1)
          end_loci <- matrix(0,nrow = length(split_end),ncol = 1)
          for(k in 1:length(split_start)){
            split_piece <- unlist(strsplit(split_start[k],''))
            start_loci[k] <- paste(split_piece[2:(length(split_piece)-1)],collapse = '')
          }
          for(k in 1:length(split_end)){
            split_piece <- unlist(strsplit(split_end[k],''))
            end_loci[k] <-paste(split_piece[2:(length(split_piece)-1)],collapse = '')
          }

          been_end <- matrix(0,nrow = length(split_end),ncol = 1)

          #for each original allele
          for(k in 1:length(split_start)){
            #if this loci is in the end genotype
            ind <- match(start_loci[k],end_loci)
            if(is.na(ind)==TRUE){#loci not in end genotype
              #record mutation rate directly from WT
              u_tot <- u_tot * Q[split_start[k],'WT']
            }else #loci in end genotype
              if(Q[split_start[k],split_end[ind]]==0){ #means multi-hit mutation required at this location
                Q[split_start[k],split_end[ind]] <- Q[split_start[k],'WT']*Q['WT',split_end[ind]]
                u_tot <- u_tot *Q[split_start[k],split_end[ind]]
              }else{
                u_tot <- u_tot*Q[split_start[k],split_end[ind]]
              }
            been_end[ind] <- 1
          }

          #for final alleles that are WT in the original genotype
          for(k in 1:length(split_end)){
            if(been_end[k]==0){
              u_tot <- u_tot*Q['WT',split_end[k]]
            }
          }

          Q[i,j] <- u_tot
        }
      }
    }
  }
  else{
    print("Error: cannot have direct_multi_multi on without having back_mutations on")
  }

}

as.table(Q)
}









