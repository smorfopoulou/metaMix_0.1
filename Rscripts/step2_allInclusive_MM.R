library(Matrix)
library(gtools)
library(DIRECT)

####need to use EM function
source("EM.R")


######### need pij info. "preprocess.RData" was produced by step1_preprocessMCMC.R
load("preprocess.RData")

EMiter<-500

noSpecies<-length(colnames(pijSparseUnknown))
print(noSpecies)
hyperP<-rep(1, noSpecies)
startW<-rdirichlet(1, hyperP)

output100Tent<-EM(pij=pijSparseUnknown, iter=EMiter, species=colnames(pijSparseUnknown), abund=startW)  ### EM function
approxSpecies0<-names(which(round(colMeans(output100Tent$abundances[EMiter,])*sum(readWeights[,"weight"]))>0))
approxSpecies0<-approxSpecies0[-1]

approxPij<-pijSparseUnknown[, approxSpecies0]
approxSpecies.with.counts<-round(colMeans(output100Tent$abundances[EMiter,2:length(colnames(output100Tent$abundances))])*sum(readWeights[,"weight"]))
ordered.approx.species<-cbind(approxSpecies.with.counts, approxSpecies.with.counts/sum(approxSpecies.with.counts))
colnames(ordered.approx.species)<-c( "countReads", "samplingWeight")
ordered.approx.species<-  data.frame("taxonID"=rownames(ordered.approx.species), ordered.approx.species, stringsAsFactors=FALSE)
orderedSpecies<-ordered.approx.species[order(-ordered.approx.species[,2]) , ]                      #### order them by read count
orderedSpecies<-orderedSpecies[which(orderedSpecies$countReads>=1),]    ###potential species are the ones that have at least one read assigned to them 

orderedSpecies<- orderedSpecies[-which(orderedSpecies$taxonID=="unknown"),]
approxSpecies<-orderedSpecies$taxonID
pijSparseUnknown<-pijSparseUnknown[,approxSpecies]


## ###Flattening the sampling probabilities
percentiles<-quantile(orderedSpecies$samplingWeight,  probs=c(0.2, 0.8))
orderedSpecies$samplingWeight[which(orderedSpecies$samplingWeight  >= percentiles["80%"])] <- percentiles["80%"]
orderedSpecies$samplingWeight[which(orderedSpecies$samplingWeight  <= percentiles["20%"])] <- percentiles["20%"]



### remove objects
rm(list= ls()[!ls() %in% c("output100Tent", "approxSpecies", "pijSparseUnknown", "approxPij", "orderedSpecies", "readWeights", "approxSpecies0")])

gc()

   
save.image("EM_species.RData")

