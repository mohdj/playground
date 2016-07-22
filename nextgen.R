sampleSet <- c("one","two","three","four","five","six")
dat <- getSimulatedBigData(N=10,sampleSet=sampleSet,numOfSingleChoiceCols=3,numOfMultiChoiceCols=2,
                           possibleNumOfItemsPerValue=1:3)

colTypes=c(x1="SingleChoice",x2="SingleChoice",x3="SingleChoice",x4="MultiChoice",x5="MultiChoice")
weight=c(x1=-5,x2=-4,x3=4,x4=5,x5=3)

groupSize <- 10
weight <- weight[weight!=0]
colTypes <- colTypes[names(weight)]
numberOfGroup <- floor(nrow(dat)/groupSize)
nr <- nrow(dat)
clusInfo <- list(weight=weight,colTypes=colTypes,numberOfGroup=numberOfGroup,
                 idealClusterSize=groupSize,datColNames=names(dat),distanceType="SimilarityBased")
dtm <- setupDocumentTermFrequencyForTextCols(dat,clusInfo)
clusInfo$columnTermFrequencyColumnMapping <- getColumnTermFrequencyColumnMapping(dtm,clusInfo)
#clusInfo$columnTermFrequencyColumnMappingIdx <- getColumnTermFrequencyColumnMappingIdx(dtm,clusInfo)
#clusInfo <- getAvgNumOfItemForMultiChoiceColumn(dat,clusInfo)

# singleChoiceColIdx <- which(clusInfo$colTypes %in% c("SingleChoice"))
# singleChoiceSimColIdx <- which(clusInfo$colTypes %in% c("SingleChoice") & clusInfo$weight > 0)
# singleChoiceDissimColIdx <- which(clusInfo$colTypes %in% c("SingleChoice") & clusInfo$weight < 0)
MultiChoiceSimColIdx <- which(clusInfo$colTypes %in% c("MultiChoice","PlainText") & clusInfo$weight > 0)
MultiChoiceDissimColIdx <- which(clusInfo$colTypes %in% c("MultiChoice","PlainText") & clusInfo$weight < 0)
MultiChoiceColIdx <- which(clusInfo$colTypes %in% c("MultiChoice","PlainText"))

matsc <- as.matrix(dtm[,which(clusInfo$colTypes %in% c("SingleChoice"))])
singleChoiceSimColIdx <- which(colnames(matsc) %in% names(weight)[weight>0])
singleChoiceDissimColIdx <- which(colnames(matsc) %in% names(weight)[weight<0])

MultiChoiceTermColNames <- unlist(clusInfo$columnTermFrequencyColumnMapping)
matmcterm <- as.matrix(dtm[,MultiChoiceTermColNames])
columnTermFrequencyColumnMappingIdx <- list()
for(i in which(clusInfo$colTypes %in% c("MultiChoice","PlainText"))){
  nm <- names(clusInfo$colTypes)[i]
  columnTermFrequencyColumnMappingIdx[[i]] <- grep(pattern = paste0("__",nm),x = colnames(matmcterm))
}

dist <- matrix(nrow = nr,ncol = nr)

system.time({for(i in 1:(nr-2)){
  tmpdist <- matrix(nrow = (nr-i),ncol = length(colTypes))
  x <- matsc[(i+1):nr,]
  y <- matsc[1:(nr-i),]
  comp <- (x==y)*1
  tmpdist[,singleChoiceSimColIdx] <- 1-comp[,singleChoiceSimColIdx]
  tmpdist[,singleChoiceDissimColIdx] <- comp[,singleChoiceDissimColIdx]
  
  for(j in MultiChoiceColIdx){
    xj <- matmcterm[(i+1):nr,columnTermFrequencyColumnMappingIdx[[j]]]
    yj <- matmcterm[1:(nr-i),columnTermFrequencyColumnMappingIdx[[j]]]
    plusxy <- xj+yj
    multiplyxy <- xj*yj
    similar <- rowSums(multiplyxy)
    totalPresent <- rowSums(plusxy) - 2*similar
    jacSim <- similar/totalPresent # Jacard similarity
    themeSim <- 1 - similar/6
    if(weight[j]>0){
      tmpdist[,j] <- pmin(themeSim,jacSim)*0.7 + 0.3*jacSim
    } else{
      tmpdist[,j] <- jacSim
    }
  }
  #tmpdist <- tmpdist*abs(weight[colnames(tmpdist)])
  #print(i)
}
})