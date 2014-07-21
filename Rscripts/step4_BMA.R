library(Matrix)
library(gtools)
library(ggplot2)

load("PT_EM.RData")
load("EM_species.RData")

nIter<-2000
unknown<-10^-20
pijSparseUnknown<-cBind(pijSparseUnknown, unknown)

source("Gibbs.R")


dyn.load('src/fast_multinom_weight.so')

fast.rmultinom.weight <- function(proba.matrix, z.matrix, seed, weights) {
    return( .Call("C_rmultinom_weight", proba.matrix, z.matrix, as.integer(seed), weights) )
  }


StartTime<-Sys.time()

message("Associate taxonIDs with scientific names: reading \"names.dmp\" is going to take some time")

taxonNames<-read.table(pipe("awk -F\"|\" 'BEGIN {OFS=\"\t\"} {if ($4~/scientific/) print $1, $2 }' /home/ucbtsmo/scratchDir/sofia/taxonomy/names.dmp | sed \"s/\t\t*/\t/g\" | awk  'BEGIN {OFS=FS=\"\t\"} {print $1,$2}'"), sep="\t", stringsAsFactors=FALSE, quote="\"")

colnames(taxonNames)<-c("taxonID","scientName" )


ofInterest<-result$slave1$record[round(nIter/5):nIter,4:(ncol(result$slave1$record)-1)]
present.probabilities<- apply(ofInterest, MARGIN=2, function(x) sum(x)/length(x))
poster.prob.all<-present.probabilities[present.probabilities>0]

poster.prob<-present.probabilities[present.probabilities>=0.9]

poster.probM<-as.data.frame(poster.prob)
poster.probM$taxonID<-rownames(poster.probM)
poster.prob<-merge(taxonNames, poster.probM, by.x="taxonID", by.y="taxonID", all.y=T)
poster.prob[which(poster.prob[,"taxonID"]=="unknown"),"scientName"]<-"unknown"
poster.prob<-poster.prob[order(-poster.prob[,"poster.prob"]),]


pdf("barplot_09_20.pdf")
par(mar=c(5.1,14.1,4.1,2.1))


 barplot(poster.prob[,"poster.prob"], names=poster.prob[,"scientName"], horiz=T,  las=2, cex.names=0.5, beside=T, col=c(rep("navy",  length(which(poster.prob[,"poster.prob"]==1))), topo.colors(nrow(poster.prob) -  length(which(poster.prob[,"poster.prob"]==1)))))
     

dev.off()

finalSpecies<-poster.prob[,"taxonID"]
noSpecies<-length(finalSpecies)

###parameters for gibbs
hyperP<-rep(1, noSpecies)
startW<-rdirichlet(1, hyperP)

BurnIn<-nIter/20
GibbsCycles<-nIter/10

message("Running final longer chain")
output1000<-gibbsSparseFinal(pij=pijSparseUnknown, iter=GibbsCycles, species=finalSpecies, abund=startW,  hyperParam=hyperP)
finalAssignments<-matrix(output1000$assignments[GibbsCycles,], ncol=1, dimnames=list(colnames(output1000$assignments[GibbsCycles,])))
finalAssignmentsDF <- data.frame(taxonID=rownames(finalAssignments), finalAssignments=unlist(finalAssignments))
finalAssignmentsDF<- finalAssignmentsDF[which(finalAssignmentsDF$taxonID!="Iter"),]


presentSpecies<-merge(taxonNames, finalAssignmentsDF, by.x="taxonID", by.y="taxonID", all.y=T)
presentSpecies[presentSpecies[,"taxonID"]=="unknown",][,"scientName"]<-"unknown"
presentSpecies<-presentSpecies[order(-presentSpecies[,"finalAssignments"]),]

presentSpecies.allInfo<-merge(presentSpecies, poster.probM, by.x="taxonID", by.y="taxonID", all.y=T)
presentSpecies.allInfo<-presentSpecies.allInfo[order(-presentSpecies.allInfo[,"finalAssignments"]),]



message("Results in \"presentSpecies_assignedReads_allInfo_09.tab\"")
write.table(presentSpecies.allInfo, "presentSpecies_assignedReads_allInfo_09_20.tab", sep="\t")



###Classification Probability

noSpecies<-nrow(presentSpecies)
mean1000<-output1000$abundances[GibbsCycles,2:(noSpecies+1)]
zij<-output1000$pijs %*% diag(mean1000)
sumProd<-rowSums(zij)
zijFinal <- as.matrix(zij / sumProd)
colnames(zijFinal)<-colnames(output1000$pijs)
zijFinal<-zijFinal[,presentSpecies[,"taxonID"]]



assignedReads<-list()
classProb<-list()
for (i in presentSpecies[,"taxonID"]){
  print(i)
  assignedReads[[i]]<-rownames(output1000$assignedReads[output1000$assignedReads[,i]>0,])
  classProb[[i]]<-zijFinal[assignedReads[[i]], i]
}


scientNames<-vector()
for (i in names(classProb)) {
  scN<-presentSpecies[presentSpecies[,"taxonID"]==i, "scientName"]
  scientNames<-append(scientNames, scN)
}
names(classProb)<-scientNames


pdf("histograms_cdf_09_20.pdf")
for (i in names(classProb)) {
  temp<-data.frame(read=names(classProb[[i]]), prob=classProb[[i]], stringsAsFactors=F)

  ploti<-ggplot(temp, aes(x=prob)) + stat_bin(aes(y=..count../sum(..count..), fill = ..count../sum(..count..))) + stat_ecdf() + labs(list(title = paste(nrow(temp),"reads assigned to ", i), x = "Classification probability", y = "Percentage of reads")) + guides(fill=guide_legend(title="Percentage of reads"))  + scale_x_continuous(breaks=c(0,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) ) + theme(plot.title = element_text(size = 12))
  print(ploti)
}
dev.off()





EndTime<-Sys.time()

duration<-EndTime-StartTime
print(EndTime - StartTime)

rm(list= ls()[!ls() %in% c("result", "pijSparseUnknown",   "presentSpecies.allInfo", "duration", "output1000", "zijFinal", "assignedReads", "classProb", "poster.prob", "poster.prob.all", "taxonNames")])


gc()


save.image("presentSpecies_barplot_09.RData")

 
