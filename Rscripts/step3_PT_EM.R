###arguments of function: preprocessed data, matrix of p_ij, all species ordered based on count of reads, contigs: nb of reads that make it
### metaMix ( blast.output, ##matrix of pij
##            weights, ##if some contigs have many reads,
##            temporary count per species)


##########################-------------------------------------------- MAIN ---------------------------------------------------------------------------------######################################
library(Rmpi)
library(Matrix)
library(gtools)
library(DIRECT)

####need to use gibbs function
source("EM.R")

print(warnings())
StartTime<-Sys.time()
noChains<-12

sieve <- function(n) {
  n <- as.integer(n)
  if(n > 1e6) stop("n too large")
  primes <- rep(TRUE, n)
  primes[1] <- FALSE
  last.prime <- 2L
  for(i in last.prime:floor(sqrt(n)))
    {
      primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
      last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
    }
          which(primes)
}


list.integers <- sieve(1000)

node.ids <- list.integers[ 1:noChains ]
### for object prime id number need to know how many objects I need to split my abundances/species into
# "noChains" slaves will start. Specify both here and on "-pe" option to qsub command

mpi.spawn.Rslaves(nslaves = noChains)  #number of slaves to spawn, should be equal to individual chains

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

######### need pij info. "preprocess.RData" was produced by step1_preprocessMCMC.R
load("EM_species.RData")

rm(output500Tent)
print(pijSparseUnknown[1,])
unknown<-10^-20
pijSparseUnknown<-cBind(pijSparseUnknown, unknown)
print(pijSparseUnknown[1,])


##########################---------------------------------------SINGLE CHAIN FUNCTION -----------------------------------------------------------------------------------------------#######################

singleChain <- function(TotalIter, exchangeInterval){

  if (file.exists('~/.Rprofile')) source('~/.Rprofile')
  print(.libPaths())
  
  library(Matrix)
  library(DIRECT)
  library(gtools)
   
  message("These are the warnings:")
  print(warnings())

  pUnknown<-1e-20
    
  # Each slave gets its own copy of ind and chain based on mpi process rank
  ind <- mpi.comm.rank()
  print(ind)

  ### Create matrix that will hold the estimDator of the log-likelihood
  estimNew<-matrix(0, nrow=TotalIter, ncol=1)

  estimNew[1,]<- sum(readWeights[,"weight"])*log(pUnknown)*temper[ind]    ### logL when all reads are unknown
print(estimNew[1,])
  
  presentSpecies<-"unknown"
  abundUsedSp<-c("unknown"=1)
  
  ### Create 2 lists. One that holds species names and one with their abundances. Each list has 2 elements. Element1: present species  and Element2: tentative species considered in current iteration.
  speciestoUse<-list("presentSp"=presentSpecies, "tentativeSp"=NULL)
  abundUsedSpecies<-list("presentSp"= abundUsedSp, "tentativeSp"=NULL)
  
  message("\nThis is the temperature for this slave")
  print(temper[ind])
  
### create matrix to hold info on which species was added or removed, which species were present and logL , through iterations
  record<-matrix(0, ncol=(5+lenSp), nrow=(TotalIter))
  record<-as.data.frame(record)
  colnames(record)<-c("Iter", "Move", "Candidate Species", allSpecies[,1],  "unknown", "logL")
  record[,1]=1:(TotalIter)

  record[1, presentSpecies]<-1
  record[1,(5+lenSp)]<-estimNew[1,]             ###last colum is logL
  
  swaps.attempted<-0
  swaps.accepted<-0

  oddFlag<-0

################ ----------------------------------- Begin MCMC (within single chain)  ------------------------------------------------------------------------------------------############  
  for (i in  2:TotalIter)  {

    cat('\n',i,'\n')                                 ### temporary - debugging purposes

#### Create 3 functions to use repeatedly in adding/removing species steps.
######################### 1a. For add step, from which species can I sample my candidate species
    potentialSpecies<-function() {                    
      toChooseFrom <-  allSpecies[,1][!(allSpecies[,1] %in% speciestoUse[[1]])]   ### Species remaining after omittig "present species"
      toChooseFrom<-as.character(toChooseFrom)
      potentialSp<- allSpecies[allSpecies$taxonID %in% toChooseFrom,]
      resultPS<-list("potentialSp"=potentialSp, "toChooseFrom"=toChooseFrom)
      return(resultPS)
    }

    ######################## 1b. For add step, sample candidate organism to add and create object for tentative species
    moveAdd<-function() {
      randSp <- sample(as.character(potentialSp[,1]),1, prob=potentialSp[,2])   ### choose random species, weights need not sum to one. Here weights are based on intial read counts.         
      cat('Add species', randSp, '\n')                     ### temporary - debugging purposes
      speciestoUse[[2]]<-c(randSp, speciestoUse[[1]])      ### tentative present species, to use in Gibbs below.
      record[i, 2]<- "Add"                     ### record info on species we are adding
      record[i,3]<- randSp

      result<-list("record"=record, "speciestoUse"=speciestoUse)      
      return(result)
    }


    ######################## 2. For remove step, sample candidate organism to remove and create object for tentative species
    moveRemove<-function() {
      toremoveFrom <- speciestoUse[[1]][!(speciestoUse[[1]]%in%"unknown")]         ### species set from which to remove (i.e present ones bar unknown)
      toremoveFrom<-as.character(toremoveFrom)          
      removeProb<-(1/abundUsedSpecies[[1]][toremoveFrom])/sum(1/abundUsedSpecies[[1]][toremoveFrom])       ###### sampling weight inversely proportional to assigned read counts.
     ###Flatten removal probabilities
      percentiles<-quantile(removeProb,  probs=c(0.2, 0.8))
      removeProb[which(removeProb  >= percentiles["80%"])] <- percentiles["80%"]
      removeProb[which(removeProb  <= percentiles["20%"])] <- percentiles["20%"]

    ### sample candidate species to remove
      randSp<-sample(toremoveFrom, 1, prob=removeProb )   ### weights need not sum to one
      cat('Remove species', randSp, '\n')                 ### temporary - debugging purposes
      speciestoUse[[2]] <- speciestoUse[[1]][!(speciestoUse[[1]] %in% randSp)]    ### tentative present species, to use in Gibbs below.
      record[i, 2]<- "Remove"                     ### record info on species we are adding
      record[i,3]<- randSp
      result<-list("record"=record, "speciestoUse"=speciestoUse)      
      return(result)
    }

    addStep <-0.5      ####proability of doing add step

    
      
              ######################################################## Add species ########################################      
    if (runif(1)<= addStep) {                             
      resultPS<-potentialSpecies()
      potentialSp<-as.data.frame(resultPS$potentialSp)
      toChooseFrom<-as.character(resultPS$toChooseFrom)
      
      if (length(toChooseFrom) != 0L){  
        result<-moveAdd()
        speciestoUse<-result$speciestoUse
        record<-result$record
        
        
      } else {                                 ################## ### If pool of species to add empty, go to remove step  
        message("No more species to choose from, all are kept as present.")
        result<-moveRemove()
        speciestoUse<-result$speciestoUse
        record<-result$record        
        addStep<-0             ### so do remove step with prob=1
      }  

    }  ###close if runif(1)<=addStep

      
    else {    ################################################### Remove species #################################################
      toremoveFrom <- speciestoUse[[1]][!(speciestoUse[[1]]%in%"unknown")]         ### species set from which to remove (i.e present ones bar unknown)
      toremoveFrom<-as.character(toremoveFrom)

         #### First check that more than 1 species are present, so don't risk of remaining only with X bin and wasting iteration.
      if (length(toremoveFrom) > 1) {
        result<-moveRemove()
        speciestoUse<-result$speciestoUse
        record<-result$record        
        
      }  else if  (length(speciestoUse[[1]])==1) {                         #### speciestoUse[[1]]==1 when only unknown bin is there / We can only add
        message("We cannot remove, we only have unknown bin")
        potentialSp<-allSpecies
        result<-moveAdd()
        speciestoUse<-result$speciestoUse
#        print(speciestoUse)
        record<-result$record        
        addStep <- 1        ### so do add step with prob=1
#        print(addStep)
        
      } else {                                                ### we are adding species instead. Removing would leave us just with X-bin again.
        message("\nCannot remove a species, so instead we add\n")   ### temporary
        resultPS<-potentialSpecies()
        potentialSp<-as.data.frame(resultPS$potentialSp)
        toChooseFrom<-as.character(resultPS$toChooseFrom)          
        result<-moveAdd()
        speciestoUse<-result$speciestoUse
        record<-result$record        
        addStep <- 1
        
      }        
    } ###close else (remove species)


    noSpecies<-length(speciestoUse[[2]])
    print(noSpecies)
    hyperP<-rep(1, noSpecies)
    startW<-rdirichlet(1, hyperP)
   
    output100Tent<-EM(pij=pijSparseUnknown, iter=EMiter, species=speciestoUse[[2]], abund=startW)  ### EM function

    estimator <- output100Tent$logL[EMiter,2] + (noSpecies* penalty) #### penalise likelihood                                  

    message("EM took: ", output100Tent$RunningTime)

  
######## parameters mean
    
    mean1<-as.numeric(output100Tent$abundances[EMiter,2:(noSpecies+1)])
    names(mean1)<-names(output100Tent$abundances[EMiter,2:(noSpecies+1)])
                        
    abundUsedSpecies[[2]]<-mean1

    estimNew[i,]<-estimator*temper[ind]                          ###tempered likelihood
    message('\n', estimNew[i,] , ' ',  estimNew[i-1,])                        ### temporary print - which values am I comparing?

    gc()
    
    ####### flag of adding/removing
    if (record[i,2]=="Add") {
      stepAdd<-TRUE
    } else {stepAdd<-FALSE}

    removeStep <- 1 - addStep    ##remove step


    if (stepAdd) {

      candidateAdd<-allSpecies[allSpecies[,1]==(as.character(record[i,3])),2]
      removeProba<-(1/abundUsedSpecies[[2]])/sum(1/abundUsedSpecies[[2]])       ###### sampling weight inversely proportional to size. Tentative Species

###Flattening the sampling probabilities
      percentiles<-quantile(removeProba,  probs=c(0.2, 0.8))
      removeProba[which(removeProba  >= percentiles["80%"])] <- percentiles["80%"]
      removeProba[which(removeProba  <= percentiles["20%"])] <- percentiles["20%"]
      candidateRemove<-removeProba[as.character(record[i,3])]

#      print(removeStep)
#      print(addStep)
      print(candidateRemove) ### gia kapoio logo twra einai NA
      print(candidateAdd)
      if (runif(1) < min( 1, exp(estimNew[i] - estimNew[i-1] + log(removeStep) - log(addStep) + log(candidateRemove) - log(candidateAdd) )))  {     ### accept "add species"  with prob min{1, P(D|i)*P(removeSpecies)*P(remove Specific Species)/P(D|i-1)*P(addSpecies)*P(add specific species) 

        ######## Accept move #######
        estimNew[i,]<-estimNew[i,]                                  ### if move is accepted, record new logL
        speciestoUse[[1]]<-speciestoUse[[2]]                        ### if move is accepted, tentative species becomes present species.
        abundUsedSpecies[[1]]<-abundUsedSpecies[[2]]                        ###       -->>-->>--   , abundances of new set of species are kept  
        cat('Present species become:', speciestoUse[[1]], '\n')
#        message("And their abundances: ", abundUsedSpecies)

      }  else {                                                    ######### Reject Move ########
        estimNew[i,]<-estimNew[i-1,]                               ### if move is rejected, record previous logL
        speciestoUse[[1]]<-speciestoUse[[1]]                       ### if move is rejected, keep present species as it is.
        abundUsedSpecies[[1]]<-abundUsedSpecies[[1]]
        cat('Present species remain:', speciestoUse[[1]], '\n')    ### temporary
 #       message("And their abundances: ", abundUsedSpecies)

      }

    }    ### close  if(stepAdd==TRUE)

    
    else {   ########################################  type of move was to remove species   ##############
      removeProba<-(1/abundUsedSpecies[[1]])/sum(1/abundUsedSpecies[[1]])   ###### sampling weight inversely proportional to size. Do this for present Species
      candidateRemove<-removeProba[as.character(record[i,3])]
      candidateAdd<-allSpecies[allSpecies[,1]==(as.character(record[i,3])),2]

      
      if (runif(1) < min( 1, exp(estimNew[i] - estimNew[i-1] + log(addStep) - log(removeStep) + log(candidateAdd) - log(candidateRemove) ))) {############# Accept "remove species"        
        estimNew[i,]<-estimNew[i,]
        speciestoUse[[1]]<-speciestoUse[[2]]                       ### tentative species becomes present species.
        abundUsedSpecies[[1]]<-abundUsedSpecies[[2]]
        cat('Present species become:', speciestoUse[[1]], '\n')    
        
      }  else {                                                    ######## Reject Move ###########
        estimNew[i,]<-estimNew[i-1,]
        speciestoUse[[1]]<-speciestoUse[[1]]
        abundUsedSpecies[[1]]<-abundUsedSpecies[[1]]
        cat('Present species remain:', speciestoUse[[1]], '\n')          
      }

    }   ###close else (i.e (stepAdd==FALSE))
    
################################# exchange between slaves
       
    if( i%%exchangeInterval == 0 ){   ### every nth (2nd for now) iteration
      message("\n\nTime to attempt an exchange")      
      oddFlag<-oddFlag+1
 #     message("This is the odd flag", oddFlag)        
      swap<-0

      estim.current<-estimNew[i,]/temper[ind]    ######### need untempered logL

      
########################################### CREATE prime tags for object to send around
      allowedLength<-175
      Nsubobjects<-round(length(abundUsedSpecies[[1]])/allowedLength)+1
      object.ids <- list.integers[ (noChains+1):(noChains + 4 + Nsubobjects) ]    ### 4 objcts for logL, swap message, untempered , species PLUS as many as necessary for abundances      

      

      if (ind%%2 == oddFlag%%2) {  ###when oddFlag zero , the following code concerns even-numbered slaves. For oddFlag 1, it concerns odd-numbered slaves.                 
        ind.partner<-ind+1

        if (0<ind.partner && ind.partner<(noChains+1)){           
          estim.partner<-mpi.recv.Robj(ind.partner, tag=object.ids[1]*node.ids[ind.partner])  #### receive the logdensity of above  partner
          message("I received the untempered: ", estim.partner)
          swaps.attempted<-swaps.attempted+1          
          lalpha<-(estim.partner - estim.current)*(temper[ind] - temper[ind.partner] )
          message("This is the acceptance probability: ", min(1,exp(lalpha)))
          
###I will check whether the exchange should happen in one of the partners
          if (runif(1)< min(1, exp(lalpha))) {    ############# exp((chain2 - chain1)*(T1 - T2))
            message("I exchanged the values")               
            swap<-1         
          }   ## end of if runif(1)<lalpha M-H step

          else {message("I didn't exchange the values")}
          
          mpi.send.Robj(obj=swap,dest=ind.partner,tag=object.ids[2]*node.ids[ind])
          message("I send message swap: ", swap)
        }     ### end of   if (0<ind.partner && ind.partner<(noChains+1)){          

        if(swap==1){
          swaps.accepted<-swaps.accepted+1

          mpi.send.Robj(obj=estim.current,dest=ind.partner,tag=object.ids[3]*node.ids[ind])
          message("I just sent the untempered: ", estim.current)
          species.swap<-speciestoUse[[1]]
          mpi.send.Robj(species.swap,dest=ind.partner,tag=object.ids[4]*node.ids[ind])
#          message("Succesfully sent species.swap")
          speciestoUse[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[4]*node.ids[ind.partner])
#          message("Succesfully received species")


          
          abund.swap<-abundUsedSpecies[[1]]
#            message("About to print the abund.swap, its size is: ", object.size(abund.swap))
                                        #print(abund.swap)
#            message("About to send abund.swap", ind.partner, ' ', ind)
          mpi.send.Robj(abund.swap,dest=ind.partner,tag=object.ids[5]*node.ids[ind])
#          message("Succesfully sent abund.swap")
          abundUsedSpecies[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[5]*node.ids[ind.partner])
#            message("Succesfully received abundances")
          

          message("I received the following abundances from: ", ind.partner, " here ", ind)
          message("The new logdensity will be the one of the neighbor but tempered with my temperature: ", estim.partner*temper[ind])
          estimNew[i,]<-estim.partner * temper[ind]
        }     ### end of if(swap==1)         
        
      } else {  ##### ###when oddFlag zero , the following code concerns odd-numbered slaves. For oddFlag 1, it concerns even-numbered slaves. I.e say what partners should do
        ind.partner<-ind-1;           
        if(0<ind.partner && ind.partner<(noChains+1)){
          mpi.send.Robj(obj=estim.current,dest=ind.partner,tag=object.ids[1]*node.ids[ind])
          message("I just sent the untempered: ", estim.current)
          swap<-mpi.recv.Robj(ind.partner,tag=object.ids[2]*node.ids[ind.partner])
          message("I received the swap message: ", swap)
          swaps.attempted<-swaps.attempted+1
        }  ####end of  if(0<ind.partner && ind.partner<(noChains+1)){
          
                  
        if(swap==1){
          swaps.accepted<-swaps.accepted+1
#          message("So far I have than many accepted swaps: ", swaps.accepted)
#          message("I am about to receive the untempered: ")
          estim.partner<-mpi.recv.Robj(ind.partner, tag=object.ids[3]*node.ids[ind.partner] )  #### receive the logdensity of above  partner
          message("I received the untempered: ", estim.partner)


          species.swap<-speciestoUse[[1]]        

          speciestoUse[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[4]*node.ids[ind.partner])           
#          message("Succesfully received species")
          mpi.send.Robj(species.swap,dest=ind.partner,tag=object.ids[4]*node.ids[ind])
#          message("Succesfully sent species.swap")

          
          abund.swap<-abundUsedSpecies[[1]]
#          message("About to print the abund.swap, its size is: ", object.size(abund.swap))
                                        #print(abund.swap)
          abundUsedSpecies[[1]]<-mpi.recv.Robj(ind.partner,tag=object.ids[5]*node.ids[ind.partner])
#          message("Succesfully received abundances")
#          message("About to send abund.swap", ind.partner, ' ', ind)
          mpi.send.Robj(abund.swap,dest=ind.partner,tag=object.ids[5]*node.ids[ind])
#          message("Succesfully sent abund.swap")
          
          
          message("I received the following species from: ", ind.partner, " here ", ind)
          message("The new logdensity will be the one of the neighbor but tempered with my temperature: ", estim.partner*temper[ind])
          estimNew[i,]<-estim.partner * temper[ind]
          
        }
      }
    }

    gc()
############ Use indicator variable for species presence per iteration
    record[i, speciestoUse[[1]]]<-ifelse(speciestoUse[[1]] %in% colnames(record), 1, 0)
    record[i,(5+lenSp)]<-estimNew[i,]             ###last colum is logL


    
  }    ###end of for loop for internal iterations

  resultSC<-list("estimNew"=estimNew, "record"=record, "usedSp"=speciestoUse[[1]], "abundUsedSp" = abundUsedSpecies[[1]], "swaps.attempt"=swaps.attempted, "swaps.accept"=swaps.accepted)                                       
  return(resultSC)

   
}    ##end of singleChain function


##########################################----------------------- MAIN ---------------------------

exchangeInterval<-1                             ###leave chains run in parallel for that many iterations
ExternIter<-1870                              ### make chains communicate 1 times
TotalIter<-exchangeInterval * ExternIter
print(TotalIter)

#### broadcast necessary objects/functions to slaves
mpi.bcast.Robj2slave(exchangeInterval)
mpi.bcast.Robj2slave(ExternIter)
mpi.bcast.Robj2slave(TotalIter)

#mpi.bcast.Robj2slave(object.ids)
mpi.bcast.Robj2slave(list.integers)
mpi.bcast.Robj2slave(node.ids)

#### These are the objects that are loaded with the load("preprocess.RData") at beginning of script
mpi.bcast.Robj2slave(orderedSpecies)
mpi.bcast.Robj2slave(readWeights)
mpi.bcast.Robj2slave(pijSparseUnknown)

############################################### Penalty Definition

total.reads<-sum(readWeights[,"weight"])
read.support<-40
penalty.pij<-matrix(data=0, ncol=2, nrow=total.reads)
LLunknown<-total.reads*log(unknown)
penalty.pij[,1]<-unknown
penalty.pij[1:read.support,2] <- 1/284332   ###median gene/genome length
penalty.abund<-c((total.reads-read.support)/total.reads , read.support/total.reads)
LLunknownPlus1<-colSums(log(penalty.pij %*% penalty.abund))
penalty<- LLunknown-LLunknownPlus1




mpi.bcast.Robj2slave(penalty)

mpi.bcast.Robj2slave(EM)

######### 
allSpecies<-orderedSpecies[,c("taxonID", "samplingWeight")]
mpi.bcast.Robj2slave(allSpecies)

lenSp<-nrow(allSpecies)
mpi.bcast.Robj2slave(lenSp)

########## EM iterations
EMiter<-10

#mpi.bcast.Robj2slave(BurnIn)
mpi.bcast.Robj2slave(EMiter)

mpi.bcast.Robj2slave(noChains)


##Tempering Vector --Power Decay
temper<-vector()
K<-0.001
a<-3/2
for (n in 2:noChains){
      temper[1]<-1
            temper[n]<-(temper[n-1]-K)^a
    }


for (i in 1:noChains) {names(temper)[i]<- paste("slave",i, sep="")}  ###names

mpi.bcast.Robj2slave(temper)


### flag for adding /removing species
stepAdd<-vector(mode = "logical")

###send function to slaves
mpi.bcast.Robj2slave(singleChain)



##################--------------------- Call the function in all children & collect results #################

### Now call the single chain function
  message('Now calling the single chains')
  result<-mpi.remote.exec(singleChain(TotalIter, exchangeInterval))
  message('Done calling the single chains')
  
  


EndTime<-Sys.time()

duration<-EndTime-StartTime

rm(list= ls()[!ls() %in% c("result", "duration")])


gc()

save.image("PT_EM.RData")

mpi.close.Rslaves(dellog=FALSE)
mpi.quit()



