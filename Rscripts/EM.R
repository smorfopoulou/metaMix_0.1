########################################################## EM algorithm #############################################################################0###
#Function takes 4 arguments:
#pij: input data with generative prob for each read and species 
#iter: number of  EMiterations
#species: taxon ids of species that are assumed present in the sample
#abund: Intial values for the mixing weights
############################################################################################################################################################

EM = function(pij, iter, species, abund)
{
  time1 <- Sys.time() 

####### Matrices to record the mixing weights and species assignments through iterations
  abundances<-matrix(0, ncol=(length(species)+1), nrow=iter)
  abundances[,1]=1:iter
  logLs<-matrix(0, ncol=2, nrow=iter)
  logLs[,1]=1:iter

  
####### Create progress bar
#  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  
#### select  the species you want and at the same time rearrange columns of pij.temp to be same with species vector, so multiplication work

  pij.temp<-pijSparseUnknown[,species]  
  timeA <- Sys.time() 

####### Begin iterations
  for (i in 1:iter) {
    w<-as.numeric(abund)
 
#####Do at the same time STEP1 and STEP2
######## STEP1.Calculate  zij {zij=pijwj/ЬЇзг_j(pijwj)} - Expectation Step
    zij<-(t(t(pij.temp)*w))
    
    sumProd<-rowSums(zij)
    zijFinal <- zij / sumProd
    rowWeights <- as.numeric(readWeights[ dimnames(zijFinal)[[1]] ,"weight"])
    zijWeighted<-zijFinal*rowWeights
    
    
######## STEP2. Compute wj - Maximization step 

    abund<- colSums(zijWeighted)/sum(colSums(zijWeighted))
    names(abund)<-colnames(zijWeighted)

### record log-likelihood
    logL<-colSums(log((pij.temp%*%abund))*readWeights[rownames(pij.temp),"weight"])


###### record output    
    abundances[i,2:(length(species)+1)]<-w   
    logLs[i,2]<-logL

#### update progress bar
 #   setTxtProgressBar(pb, i)

    
} ###end of iterat

###close progress bar
#  close(pb)
  
### Processing time
  time2<-Sys.time()
  timeDiff<-time2-time1
  iterTimeDiff<-time2-timeA

  abundances=data.frame(abundances)
  names(abundances)<-c("Iter", colnames(zij))
  
  result<-list("abundances"=abundances, "logL"=logLs, "RunningTime"=timeDiff)
  return(result)
}
