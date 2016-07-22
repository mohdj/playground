library(dplyr)
library(tm)
library(stringr)
library(data.table)

# Generate Synthetic data ------------------
generateRandomSample <- function(sampleSetToGenerateFrom,numOfObsToGenerate){
  # Like below
  # dat <- data.frame(x1=sample(c(0,1),numOfObs,T),x2=sample(c(1,2,3),numOfObs,T))
  strcode <- paste0("data.frame(",paste(paste0(names(sampleSetToGenerateFrom),"=sample(sampleSetToGenerateFrom[[",1:length(sampleSetToGenerateFrom),"]],numOfObsToGenerate,T)"),
                             collapse=","),")")
  out <- eval(parse(text=strcode))
}
generateCompatibleVector <- function(x,weight,numOfObsToGenerate,sampleSetToGenerateFrom){
  my_sample <- function(x,size){
    if(length(x)>1){
      out <- sample(x,size,T)
    }
    else {
      out <- rep(x,size)
    }
    return(out)
  }
  colNames <- names(sampleSetToGenerateFrom)
  strcode <- paste0("data.frame(",paste0(paste0(colNames,"=my_sample(sampleSetToGenerateFrom[[",1:length(sampleSetToGenerateFrom),"]],numOfObsToGenerate)"),collapse=","),")")
  out <- eval(parse(text=strcode))
  return(out)
}

generateCollabData <- function(sampleSetToGenerateFrom,weight,numOfCluster,numOfRandomObs,numOfObservation){
  out <- list()
  
  if(numOfCluster > 0){
    clusterDat <- generateRandomSample(sampleSetToGenerateFrom,numOfCluster)
    if(length(numOfObservation)==1) numOfObservation <- rep(x = numOfObservation,numOfCluster)
    out <- list()
    for(i in 1:nrow(clusterDat)){
      # Get a sample set from which to pull (for equality condition it should be one element, for other it can be more but lesser by 1)
      tmpSampeSet <- list()
      x <- clusterDat[i,]
      for(j in 1:length(sampleSetToGenerateFrom)){
        xVal <- x[[names(sampleSetToGenerateFrom)[j]]]
        if(weight[j] < 0){
          vals <- sampleSetToGenerateFrom[[j]]
          #vals <- vals[which(vals!=xVal)]
          tmpSampeSet[[names(sampleSetToGenerateFrom)[j]]] <- vals
        } else {
          tmpSampeSet[[names(sampleSetToGenerateFrom)[j]]] <- xVal
        }
      }
      out[[i]] <- generateCompatibleVector(x,weight,numOfObservation[i],tmpSampeSet)
    }
    out <- rbind_all(out)
    out <- rbind_list(clusterDat,out)
  }
  
  if(numOfRandomObs > 0){
    tmp <- generateRandomSample(sampleSetToGenerateFrom,numOfRandomObs)
    out <- rbind_list(out,tmp)
  }
  
  return(out)
}

# Perform Cluster Validation ------------
getClusterError <- function(clusInfo){
  avgclusDist <- 0
  for(clusNo in 1:length(clusInfo$memberClusterMap)){
    memberIdxInClus <- clusInfo$memberClusterMap[[clusNo]]
    avgclusDist <- avgclusDist+getAvgClusterDistance(memberIdxInClus,clusInfo)*length(memberIdxInClus)
  }
  return(avgclusDist)
}

getClusteringPerformance <- function(dat,clusInfo){
  clusterLabels <- unique(clusInfo$clusterNumbers)
  
  colNamesToConsider <- names(clusInfo$weight != 0)
  clusInfo$colTypes <- clusInfo$colTypes[colNamesToConsider]
  factorOrTimeZoneColNames <- names(clusInfo$colTypes)[which(clusInfo$colTypes %in% c("SingleChoice","TimeZone"))]
  
  keyVal <- data.frame(clusterNumber=numeric(),colName=character(),Proportion=numeric())
  for(clusterNum in clusterLabels){
    datInClus <- dat[which(clusInfo$clusterNumbers==clusterNum),]
    keyVal <- rbind(keyVal,data.frame(clusterNumber=clusterNum,colName="GroupSize",Proportion=nrow(datInClus)))
    
    for(nm in factorOrTimeZoneColNames){
      tmptable <- table(datInClus[[nm]])
      if(length(tmptable)>0){ # that means there are some elements, all elements are not NA
        memberProportion <- data.frame(colValue=character(),value=numeric())
        for(i in 1:length(tmptable)){
          if(length(names(tmptable)[i])==0) {
            browser()
          }
          val <- round(as.numeric(tmptable[i])/nrow(datInClus),digits = 2)
          if(is.na(val)) val <- 0
          memberProportion <- rbind(memberProportion,data.frame(colValue=names(tmptable)[i],value=val))
          
        }
        proportion <- paste(paste(as.character(memberProportion$colValue),memberProportion$value,sep=":"),collapse=",")
        keyVal <- rbind(keyVal,data.frame(clusterNumber=clusterNum,colName=nm,Proportion=proportion))
      }
    }
    
    listColNames <- names(clusInfo$colTypes)[which(clusInfo$colTypes %in% c("MultiChoice","PlainText"))]
    for(nm in listColNames){
      proportions <- round(sort(colMeans(datInClus[,clusInfo$columnTermFrequencyColumnMapping[[nm]]]),decreasing = T),digits = 2)
      proportions <- proportions[which(proportions>0)]
      names(proportions) <- gsub(pattern = paste0("__",nm),replacement = "",x = names(proportions))
      proportions <- paste(paste(names(proportions),proportions,sep=":"),collapse=",")
      keyVal <- rbind(keyVal,data.frame(clusterNumber=clusterNum,colName=nm,Proportion=proportions))
    }
    
  }
  #keyVal$clusterNumber <- paste0("Cluster_",keyVal$clusterNumber)
  #keyVal$value <- round(keyVal$value,digits = 2)
  #keyVal$value[is.na(keyVal$value)] <- 0
  #keyVal <- tidyr::spread(data = keyVal,key = clusterNumber,value = Proportion)
  keyVal <- keyVal %>% select(ColumnName=colName,TeamSet=clusterNumber,Proportion) %>% arrange(TeamSet)
  return(keyVal)
}

# Determine distances ---------------
populateDistanceMatrix <- function(dat,distFunc,clusInfo){
  numOfRow <- nrow(dat)
  datMat <- as.matrix(dat)
  mydist <- matrix(nrow=numOfRow,ncol=numOfRow)
  for(i in 1:numOfRow){
    for(j in i:numOfRow){
      val <- distFunc(datMat[i,],datMat[j,],clusInfo)
      mydist[i,j] <- val
      mydist[j,i] <- val 
    }
  }
  return(mydist)
}

getCosineSim <- function(x,y){
  c <- sum(x*y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
  return(c)
}

getDist <- function(x,y,clusInfo){
  dists <- numeric()
  weight <- clusInfo$weight
  colTypes <- clusInfo$colTypes
  columnTermFrequencyColumnMapping <- clusInfo$columnTermFrequencyColumnMapping
  importantTagDetail <- clusInfo$importantTagDetail
  
  colNamesToConsider <- names(weight != 0)
  colTypes <- colTypes[colNamesToConsider]
  
  # calculate distance of SingleChoice columns
  factorCols <- names(colTypes)[which(colTypes=="SingleChoice")]
  comp <- as.numeric(x[factorCols]==y[factorCols])
  names(comp) <- factorCols
  factorCols4Similarity <- names(which(weight[factorCols]>0))
  factorCols4DisSimilarity <- names(which(weight[factorCols]<0))
  simDist <- 1-comp[factorCols4Similarity]
  dissimDist <- comp[factorCols4DisSimilarity]
  dists <- c(simDist,dissimDist)
  # missing value treatment - assign max distance if now known
  dists[is.na(dists)] <- 1
  
  # calculate distance for timezone column
  timzoneColName <- names(colTypes)[which(colTypes=="TimeZone")]
  if(length(timzoneColName)>0){
    timediff <- abs(as.numeric(x[timzoneColName])-as.numeric(y[timzoneColName]))
    if(timediff==0) dists[timzoneColName] <- 0
    else if(timediff<=1) dists[timzoneColName]=0.2
    else if(timediff<=2) dists[timzoneColName] = 0.9
    else if(timediff > 2) dists[timzoneColName] = 1  
  }
  
  # calculate distance of list columns
  listCols <- names(colTypes)[which(colTypes %in% c("MultiChoice","PlainText"))]
  for(i in listCols){
    if(weight[i]>0){
      colnamestmp <- columnTermFrequencyColumnMapping[[i]]
      simDist <- getSimilarityBasedDistance(x[colnamestmp],y[colnamestmp],clusInfo$avgNumOfItemForMultiChoiceColumn[[i]],clusInfo)
      jacdist <- getJaccardDistance(x[colnamestmp],y[colnamestmp])
      simdist <- min(simDist,jacdist)
      dists[i] <- 0.7*simdist + 0.3*jacdist
    } else if(weight[i]<0){
      colnamestmp <- columnTermFrequencyColumnMapping[[i]]
      dists[i] <- 1-getJaccardDistance(x[colnamestmp],y[colnamestmp])
    }
  }
  
  xydist <- sum(dists*abs(weight[names(dists)]))
  # also try weighted harmonic mean or geometric mean, because mean may try to bias certian things
  return(xydist)
}

getSimilarityBasedDistance <- function(x,y,maxNumOfSimilarity,clusInfo){
  x <- x>0
  y <- y>0
  cfm <- table(x,y)
  if(nrow(cfm)==2 & ncol(cfm)==2){
    #maxNumOfSimilarity <- min(length(x),5)
    dist <- 1-cfm["TRUE","TRUE"]/maxNumOfSimilarity
    dist <- ifelse(dist<0,0,dist)
  } else {
    dist <- 1
  }
  return(dist)
}

getJaccardDistance <- function(x,y){
  x <- as.numeric(x>0)
  y <- as.numeric(y>0)
  # refer wikipedia "jaccard index" for notations
  cfm <- table(x,y)
  if(nrow(cfm)==2 & ncol(cfm)==2){
    jaccardDist <- (cfm[1,2]+cfm[2,1])/(cfm[1,2]+cfm[2,1]+cfm[2,2])
  } else {
    jaccardDist <- 1
  }
  return(jaccardDist)
}

# generate Clusters --------------
generateCluster <- function(dat,clusInfo,typeOfClustering="hclust"){
  clusInfo$typeOfClustering <- typeOfClustering
  mydist <- populateDistanceMatrix(dat,getDist,clusInfo)
  clusInfo[["distMatrix"]] <- mydist
  if(typeOfClustering=="hclust"){
    clus2 <- cluster::agnes(x = as.dist(mydist),diss = T)
    clusterNumbers <- cutree(clus2, k=clusInfo$numberOfGroup) # cut tree into 5 clusters
  }else if(typeOfClustering=="pam"){
    clus1 <- cluster::pam(x = as.dist(mydist),diss = T,k = clusInfo$numberOfGroup) 
    clusterNumbers <- clus1$clustering
  }
  clusInfo[["clusterNumbers"]] <- clusterNumbers
  clusterIndices <- list()
  for(clusNo in unique(clusterNumbers)){
    clusterIndices[[clusNo]] <- which(clusterNumbers==clusNo)
  }
  clusInfo[["memberClusterMap"]] <- clusterIndices
  
  clusInfo[["clusterPerformance"]] <- getClusteringPerformance(dat,clusInfo)
  clusInfo[["ClusterError"]] <- getClusterError(clusInfo)
  
  return(clusInfo)
}

setupDocumentTermFrequencyForTextCols <- function(dat,clusInfo){
  listSeparator <- ","
  colNames <- names(clusInfo$colTypes)
  # List from a set
  listfromsetidx <- colNames[which(clusInfo$colTypes=="MultiChoice")]
  for(i in listfromsetidx){
    dat[[i]][is.na(dat[[i]])] <- "" # make NAs as empty strings
    tmp <- gsub(pattern = " ",replacement = "",x = dat[[i]])
    tmp <- gsub(pattern = listSeparator,replacement = " ",x = tmp)
    askdoc <- VectorSource(tmp)
    askdoc <- VCorpus(askdoc)
    askdoc <- tm_map(askdoc, tolower)   # *Converting to lowercase:*    
    askdoc <- tm_map(askdoc, stripWhitespace)   # *Stripping whitespace   
    askdoc <- tm_map(askdoc, PlainTextDocument)  
    dtm <- DocumentTermMatrix(askdoc)  
    dtm <- as.data.frame(as.matrix(dtm))
    
    termColumnSums <- colSums(dtm,na.rm = T)
    termColumnSums <- termColumnSums[termColumnSums>1]
    dtm <- dtm[,names(termColumnSums)]
    
    row.names(dtm) <- NULL
    names(dtm) <- paste0(names(dtm),"__",i)
    dat <- cbind(dat,dtm)
  }
  # Plain text
  plainTextidx <- colNames[which(clusInfo$colTypes=="PlainText")]
  for(i in plainTextidx){
    dat[[i]][is.na(dat[[i]])] <- "" # make NAs as empty strings
    #tmp <- gsub(pattern = " ",replacement = "",x = dat[[i]])
    #tmp <- gsub(pattern = listSeparator,replacement = " ",x = tmp)
    askdoc <- VectorSource(dat[[i]])
    askdoc <- VCorpus(askdoc)
    askdoc <- tm_map(askdoc, removePunctuation)
    askdoc <- tm_map(askdoc, removeNumbers)
    askdoc <- tm_map(askdoc, tolower)   # *Converting to lowercase:*    
    askdoc <- tm_map(askdoc, removeWords, getStopWordsCollabAI())
    askdoc <- tm_map(askdoc, stemDocument)
    askdoc <- tm_map(askdoc, stripWhitespace)   # *Stripping whitespace   
    askdoc <- tm_map(askdoc, PlainTextDocument)  
    dtm <- DocumentTermMatrix(askdoc)  
    #dtm <- removeSparseTerms(dtm, 0.1)
    dtm <- as.data.frame(as.matrix(dtm))
    dtm <- reduceDocumentTermMatrixDataFrame(dtm)
    row.names(dtm) <- NULL
    names(dtm) <- paste0(names(dtm),"__",i)
    dat <- cbind(dat,dtm)
  }
  return(dat)
}

getColumnTermFrequencyColumnMapping <- function(dat,clusInfo){
  listCols <- names(clusInfo$colTypes)[which(clusInfo$colTypes %in% c("MultiChoice","PlainText"))]
  columnTermFrequencyColumnMapping <- list()
  for(i in listCols){
    columnTermFrequencyColumnMapping[[i]] <- grep(pattern = paste0("__",i),x = names(dat),value = T)
  }
  return(columnTermFrequencyColumnMapping)
}

getAvgNumOfItemForMultiChoiceColumn <- function(dat,clusInfo){
  avgNumOfItem <- list()
  for(nm in names(clusInfo$columnTermFrequencyColumnMapping)){
    avgNumOfItem[[nm]] <- round(mean(rowSums(dat[,clusInfo$columnTermFrequencyColumnMapping[[nm]]]),na.rm = T),digits = 0)
  }
  clusInfo$avgNumOfItemForMultiChoiceColumn <- avgNumOfItem
  return(clusInfo)
}

createGroups <- function(dat,weight,colTypes,groupSize,importantTags="",assignTagWeightByFreq=T){
  numberOfGroup <- floor(nrow(dat)/groupSize)
  weight <- weight[weight!=0]
  typeOfClustering <- "hclust"
  colTypes <- colTypes[names(weight)]
  weight <- getAlgoWeights(weight)
  #names(dat) <- gsub(pattern = "\\.",replacement = "",x = names(dat))
  clusInfo <- list(weight=weight,colTypes=colTypes,numberOfGroup=numberOfGroup,
                   idealClusterSize=groupSize,datColNames=names(dat),distanceType="SimilarityBased")
  dat <- setupTimeZoneColumn(dat,clusInfo)
  dat <- setupDocumentTermFrequencyForTextCols(dat,clusInfo)
  clusInfo$columnTermFrequencyColumnMapping <- getColumnTermFrequencyColumnMapping(dat,clusInfo)
  clusInfo <- getAvgNumOfItemForMultiChoiceColumn(dat,clusInfo)
  clusInfo <- generateCluster(dat = dat,clusInfo=clusInfo,typeOfClustering = typeOfClustering)
  #write.csv(clusInfo$clusterPerformance,"perfBeforeClassBalance.csv")
  clusInfo <- balanceClusterSize(clusInfo)
  #write.csv(clusInfo$clusterPerformance,"perfAfterClassBalance.csv")
  clusInfo <- performMovementsToImproveClusters(clusInfo)
  clusInfo$clusterPerformance <- getClusteringPerformance(dat,clusInfo)
  return(clusInfo)
}

setupTimeZoneColumn <- function(dat,clusInfo){
  timeZoneMap <- data.frame(timezone=c('A','B','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','India')
                            ,utctz=c(1:12,-1:-12,0,5.5),stringsAsFactors = F)
  timzoneColName <- names(clusInfo$colTypes)[which(clusInfo$colTypes=="TimeZone")]
  if(length(timzoneColName)>0){
    dat[[timzoneColName]] <- as.character(dat[[timzoneColName]])
    codestr <- paste0("dat %>% inner_join(timeZoneMap,by=c('",timzoneColName,"'='timezone'))")
    dat <- eval(parse(text=codestr))
    dat[[timzoneColName]] <- dat$utctz
    dat$utctz <- NULL  
  }
  return(dat)
}

getAlgoWeights <- function(sliderWeight){
  # concept of weight distribution
  # starts at 25, in next doubles (50), than triples (150), than quadruples (600), than 5 times (3000)
  # this happens on both sides, but extremely high importance is given to -5, since that means hard cutoff
  #algoWeightMap <- data.frame(sliderWeight=-5:5,algoWeight=c(-20000,-600,-150,-50,-25,0,25,50,150,600,3000))
  algoWeightMap <- data.frame(sliderWeight=-5:5,algoWeight=c(-10000,-100,-10,-1,-0.1,0,0.1,1,10,100,1000))
  tmp <- data.frame(sliderWeight=sliderWeight) %>% inner_join(algoWeightMap,by="sliderWeight")
  weight <- tmp$algoWeight
  names(weight) <- names(sliderWeight)
  return(weight)
}

getStopWordsCollabAI <- function(){
  words <- c("learn","like","play","love","around","better","make","think","able","best")
  words <- c(words,stopwords("english"))
  return(words)
}

reduceDocumentTermMatrixDataFrame <- function(dat){
  termInMoreThan1Doc <- colSums(dat)
  termInMoreThan1Doc <- termInMoreThan1Doc[termInMoreThan1Doc>1]
  dat <- dat[,names(termInMoreThan1Doc)]
  return(dat)
}

# Misc --------------
identifyColumnTypes <- function(dat){
  colnm <- names(dat)
  colType <- list()
  sampsize <- ifelse(nrow(dat)>100,100,nrow(dat))
  for(nm in colnm){
    x <- sample(x = dat[[nm]],size = sampsize)
    comp <- sum(grepl(pattern = ",",x = x))
    if(comp>0) colType[[nm]] <- "MultiChoice"
    else colType[[nm]] <- "SingleChoice"
  }
  return(colType)
}