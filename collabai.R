library(dplyr)
library(tm)
library(stringr)
library(data.table)

# Member Movement - Functions to be used later ----------------------
getAvgClusterDistance <- function(memberIdxInClus,clusInfo){
  if(length(memberIdxInClus)==1) 
    return(0)
  tmpdist <- clusInfo$distMatrix[memberIdxInClus,memberIdxInClus]
  for(idx in 1:length(memberIdxInClus)){
    tmpdist[idx,idx] <- NA
  }
  avgDist <- mean(tmpdist,na.rm = T)
  return(avgDist)
}

getAvgClusterDistanceForAllClusters <- function(clusInfo,clusNos=NULL){
  if(is.null(clusNos)) clusNos <- 1:length(clusInfo$memberClusterMap)
  out <- sapply(clusNos,function(clusno,clusInfo){
    return(getAvgClusterDistance(clusInfo$memberClusterMap[[clusno]],clusInfo))
  },clusInfo)
  return(out)
}

getVarianceIncAfterRemoval <- function(memberIdxToBeRemoved,clusNo,clusInfo){
  memberIdxInClus <- clusInfo$memberClusterMap[[clusNo]]
  origAvgDist <- getAvgClusterDistance(memberIdxInClus,clusInfo)
  afterRemovalAvgDist <- getAvgClusterDistance(memberIdxInClus[which(memberIdxInClus != memberIdxToBeRemoved)],clusInfo)
  return(afterRemovalAvgDist-origAvgDist)
}

getVarianceIncAfterAdding <- function(memberIdxToBeAdded,clusNo,clusInfo){
  memberIdxInClus <- clusInfo$memberClusterMap[[clusNo]]
  origAvgDist <- getAvgClusterDistance(memberIdxInClus,clusInfo)
  afterAddingAvgDist <- getAvgClusterDistance(c(memberIdxInClus,memberIdxToBeAdded),clusInfo)
  return(afterAddingAvgDist-origAvgDist)
}

getVarianceIncreaseForMemberMovement <- function(memberIdx,clusterToBeAddedTo,clusInfo){
  clusterToBeRemovedFrom <- clusInfo$clusterNumbers[memberIdx]
  varianceIncForRemoval <- getVarianceIncAfterRemoval(memberIdx,clusterToBeRemovedFrom,clusInfo)
  varianceIncForAdding <- getVarianceIncAfterAdding(memberIdx,clusterToBeAddedTo,clusInfo)
  varianceInc <- varianceIncForRemoval+varianceIncForAdding
}

determineTheClusterToMoveTo <- function(memberIdx,undersizeCluster,clusInfo){
  #for the given member
  avgdist <- numeric()
  for(clusNo in undersizeCluster){
    tmp <- getAvgDistanceFromOtherClusterMembers(memberIdx,clusInfo$distMatrix,clusInfo$memberClusterMap[[clusNo]])
    avgdist <- c(avgdist,tmp)
  }
  return(undersizeCluster[which.min(avgdist)])
}

fullGreedyMemberMovement <- function(clusterIndices,distMatrix){
  #   Otherwise you can go for greedy approach as of now whereby you find avgdistance of a member with other
  #   cluster members and choose (clusterSize-N) top members sorted by desc distance.
  #   Than go one by one in moving these items to undersize clusters
  #   Here again you can go via greedy way (one by one from the top)
}

GreedyListButOptimalMovement <- function(){
  #   Otherwise you can go for greedy approach as of now whereby you find avgdistance of a member with other
  #   cluster members and choose (clusterSize-N) top members sorted by desc distance.
  # followed by exhaustive search  (all the time consider all possible movements of all available members)
}

fullOptimalMemberMovement <- function(){
  # Consider all the members of all the oversize clusters and all their possible movements to all 
  # undersize clusters, choose the movement which has least increase in overall variance of both the
  # affected cluster
}

determineTheBestTransition_optimal <- function(memberIndicesToBeConsidered,undersizeCluster,clusInfo){
  #for all the preselected members find the avg distance each element has wrt to other clusters 
  # and choose the one with minimum distance
  besttransition <- data.frame(memberIdx=integer(),clusterNo=integer(),avgdist=numeric())
  for(clusNo in undersizeCluster){
    avgdists <- sapply(X = memberIndicesToBeConsidered, FUN = getVarianceIncreaseForMemberMovement,clusNo,clusInfo)
    besttransition <- rbind(besttransition,data.frame(memberIdx=memberIndicesToBeConsidered[which.min(avgdists)],clusterNo=clusNo,avgdist=min(avgdists)))
  }
  return(besttransition[which.min(besttransition$avgdist),])
}

fillUndersizeClusters_optimal <- function(clusInfo){
  # Otherwise you can go for greedy approach as of now whereby you find avgdistance of a member with other
  #   cluster members and choose (clusterSize-N) top members sorted by desc distance.
  
  # 1. Determine exraCapacityCluster - cluster which can donate members & 
  #     Determine undersize clusters
  #     Determine how many movements have to be made to undersize clusters in order to make them 
  #     adequately sized
  extraCapacitycluster <- numeric()
  undersizeCluster <- numeric()
  numberOfMovementsReqdForUndersizeClus <- 0
  for(clusNo in 1:length(clusInfo$memberClusterMap)){
    if(length(clusInfo$memberClusterMap[[clusNo]]) > clusInfo$idealClusterSize)
      extraCapacitycluster <- c(extraCapacitycluster,clusNo)
    else if(length(clusInfo$memberClusterMap[[clusNo]]) < clusInfo$idealClusterSize){
      undersizeCluster <- c(undersizeCluster,clusNo)
      numberOfMovementsReqdForUndersizeClus = numberOfMovementsReqdForUndersizeClus + clusInfo$idealClusterSize - length(clusInfo$memberClusterMap[[clusNo]])
    }
  }
  # 2. Determine the list of members which needs to be removed - full greedy way
  #     Get members from extra capacity cluster which have overall higher avg distances from their cluster members
  membersToBeMoved <- integer()
  for(clusNo in extraCapacitycluster){
    membersToBeMoved <- c(membersToBeMoved,clusInfo$memberClusterMap[[clusNo]])
  }
  
  
  # Iteratively move all the members one by one to undersize clusters
  while(numberOfMovementsReqdForUndersizeClus > 0){
    transition <- determineTheBestTransition_optimal(membersToBeMoved,undersizeCluster,clusInfo) 
    # Remove all the members of this cluster if this cluster no longer have capacity to donate
    membersOrigClusNo <- clusInfo$clusterNumbers[transition$memberIdx]
    origClusterSize <- length(clusInfo$memberClusterMap[[membersOrigClusNo]])
    if((origClusterSize-1)==clusInfo$idealClusterSize)
      membersToBeMoved <- membersToBeMoved[which(!(membersToBeMoved %in% clusInfo$memberClusterMap[[membersOrigClusNo]]))]
    # Update clusInfo
    clusInfo <- updateClusInfoForTransition(transition$memberIdx,transition$clusterNo,clusInfo)
    # Remove current cluster from undersize cluster if its adequately filled
    if(length(clusInfo$memberClusterMap[[transition$clusterNo]]) == clusInfo$idealClusterSize)
      undersizeCluster <- undersizeCluster[which(undersizeCluster != transition$clusterNo)]
    # Remove current member from the list of membersToBeMoved
    membersToBeMoved <- membersToBeMoved[which(membersToBeMoved != transition$memberIdx)]
    numberOfMovementsReqdForUndersizeClus <- numberOfMovementsReqdForUndersizeClus-1
  }
  
  return(clusInfo)
}

# Member Movement - Current cluster balancing code --------------------

getAvgDistanceFromOtherClusterMembers <- function(memberIdx,clusInfo,clusNo=NULL){
  if(is.null(clusNo)){
    clusterNoOfMember <- clusInfo$clusterNumbers[memberIdx]
    otherMemberIdxInClus <- clusInfo$memberClusterMap[[clusterNoOfMember]]
    otherMemberIdxInClus <- otherMemberIdxInClus[which(otherMemberIdxInClus != memberIdx)]  
  } else{
    otherMemberIdxInClus <- clusInfo$memberClusterMap[[clusNo]]
  }
  tmpdist <- clusInfo$distMatrix[memberIdx,otherMemberIdxInClus]
  avgdist <- mean(tmpdist)
  return(avgdist)
}

fillUndersizeClusters <- function(clusInfo){
  # Otherwise you can go for greedy approach as of now whereby you find avgdistance of a member with other
  #   cluster members and choose (clusterSize-N) top members sorted by desc distance.
  
  # 1. Determine exraCapacityCluster - cluster which can donate members & 
  #     Determine undersize clusters
  #     Determine how many movements have to be made to undersize clusters in order to make them 
  #     adequately sized
  extraCapacitycluster <- numeric()
  undersizeCluster <- numeric()
  numberOfMovementsReqdForUndersizeClus <- 0
  for(clusNo in 1:length(clusInfo$memberClusterMap)){
    if(length(clusInfo$memberClusterMap[[clusNo]]) > clusInfo$idealClusterSize)
      extraCapacitycluster <- c(extraCapacitycluster,clusNo)
    else if(length(clusInfo$memberClusterMap[[clusNo]]) < clusInfo$idealClusterSize){
      undersizeCluster <- c(undersizeCluster,clusNo)
      numberOfMovementsReqdForUndersizeClus = numberOfMovementsReqdForUndersizeClus + clusInfo$idealClusterSize - length(clusInfo$memberClusterMap[[clusNo]])
    }
  }
  # 2. Determine the list of members which needs to be removed - full greedy way
  #     Get members from extra capacity cluster which have overall higher avg distances from their cluster members
  membersWhichCanBeRemoved <- list()
  for(clusNo in extraCapacitycluster){
    avgdists <- sapply(X = clusInfo$memberClusterMap[[clusNo]], FUN = getAvgDistanceFromOtherClusterMembers,clusInfo)
    avgdists <- data.frame(memberIdx=clusInfo$memberClusterMap[[clusNo]],avgdist=avgdists)
    avgdists <- avgdists %>% arrange(desc(avgdist))
    avgdists <- avgdists[1:(nrow(avgdists)-clusInfo$idealClusterSize),]
    membersWhichCanBeRemoved[[clusNo]] <- avgdists
  }
  membersWhichCanBeRemoved <- rbind_all(membersWhichCanBeRemoved)
  membersWhichCanBeRemoved <- membersWhichCanBeRemoved %>% arrange(desc(avgdist))
  membersToBeMoved <- membersWhichCanBeRemoved$memberIdx[1:numberOfMovementsReqdForUndersizeClus]
  
  # Iteratively move all the members one by one to undersize clusters
  while(length(membersToBeMoved) > 0 & length(undersizeCluster) > 0){
    transition <- determineTheBestTransition(membersToBeMoved,undersizeCluster,clusInfo)  
    clusInfo <- updateClusInfoForTransition(transition$memberIdx,transition$clusterNo,clusInfo)
    # Remove current cluster from undersize cluster if its adequately filled
    if(length(clusInfo$memberClusterMap[[transition$clusterNo]]) == clusInfo$idealClusterSize)
      undersizeCluster <- undersizeCluster[which(undersizeCluster != transition$clusterNo)]
    # Remove current member from the list of membersToBeMoved
    membersToBeMoved <- membersToBeMoved[which(membersToBeMoved != transition$memberIdx)]
  }
  
  return(clusInfo)
  
}

determineTheBestTransition <- function(memberIndicesToBeConsidered,undersizeCluster,clusInfo){
  #for all the preselected members find the avg distance each element has wrt to other clusters 
  # and choose the one with minimum distance
  besttransition <- data.frame(memberIdx=integer(),clusterNo=integer(),avgdist=numeric())
  for(clusNo in undersizeCluster){
    avgdists <- sapply(X = memberIndicesToBeConsidered, FUN = getAvgDistanceFromOtherClusterMembers,clusInfo,clusNo)
    besttransition <- rbind(besttransition,data.frame(memberIdx=memberIndicesToBeConsidered[which.min(avgdists)],clusterNo=clusNo,avgdist=min(avgdists)))
  }
  return(besttransition[which.min(besttransition$avgdist),])
}

updateClusInfoForTransition <- function(memberIdxTobeMoved,clusterToBeMovedTo,clusInfo){
  # Only these two needs to be updated - clusterNumbers, memberClusterMap
  originalClusNo <- clusInfo$clusterNumbers[memberIdxTobeMoved]
  origClusMembers <- clusInfo$memberClusterMap[[originalClusNo]]
  clusInfo$memberClusterMap[[originalClusNo]] <- origClusMembers[which(origClusMembers!=memberIdxTobeMoved)]
  clusInfo$memberClusterMap[[clusterToBeMovedTo]] <- c(clusInfo$memberClusterMap[[clusterToBeMovedTo]],memberIdxTobeMoved)
  clusInfo$clusterNumbers[memberIdxTobeMoved] <- clusterToBeMovedTo
  return(clusInfo)
}

removeFromOversizedCluster <- function(clusInfo){
  # Identify oversize & possibleSinkCluster cluster
  oversizeCluster <- numeric()
  possibleSinkCluster <- numeric()
  for(clusNo in 1:length(clusInfo$memberClusterMap)){
    if(length(clusInfo$memberClusterMap[[clusNo]]) > clusInfo$idealClusterSize+1)
      oversizeCluster <- c(oversizeCluster,clusNo)
    else if(length(clusInfo$memberClusterMap[[clusNo]]) == clusInfo$idealClusterSize)
      possibleSinkCluster <- c(possibleSinkCluster,clusNo)
  }
  
  # Determine members which should be removed - currently usual greedy way
  membersToBeMoved <- integer()
  for(clusNo in oversizeCluster){
    avgdists <- sapply(X = clusInfo$memberClusterMap[[clusNo]], FUN = getAvgDistanceFromOtherClusterMembers,clusInfo)
    avgdists <- data.frame(memberIdx=clusInfo$memberClusterMap[[clusNo]],avgdist=avgdists)
    avgdists <- avgdists %>% arrange(desc(avgdist))
    membersToBeMoved <- c(membersToBeMoved,avgdists$memberIdx[1:(nrow(avgdists)-clusInfo$idealClusterSize-1)])
  }
  
  # Iteratively move all the members one by one from oversize cluster to possibleSinkClusters
  while(length(membersToBeMoved) > 0 & length(possibleSinkCluster)>0){
    transition <- determineTheBestTransition(membersToBeMoved,possibleSinkCluster,clusInfo)  
    clusInfo <- updateClusInfoForTransition(transition$memberIdx,transition$clusterNo,clusInfo)
    # Remove current cluster from possible sink cluster list, since no more capacity
    possibleSinkCluster <- possibleSinkCluster[which(possibleSinkCluster != transition$clusterNo)]
    # Remove current member from the list of membersToBeMoved
    membersToBeMoved <- membersToBeMoved[which(membersToBeMoved != transition$memberIdx)]
    
  }
  
  return(clusInfo)
}

balanceClusterSize <- function(clusInfo){
  #clusInfo <- fillUndersizeClusters(clusInfo) 
  clusInfo <- fillUndersizeClusters_optimal(clusInfo) 
  clusInfo <- removeFromOversizedCluster(clusInfo)
  return(clusInfo)
}

# Randomized Improvement --------------
getAllMemberAvgDistancesFromGivenClusters <- function(clusNumbers,clusInfo){
  alldists <- as.data.frame(matrix(nrow = length(clusInfo$clusterNumbers),ncol = 1+length(clusNumbers)))
  names(alldists) <- c("memberIdx",paste0("c",clusNumbers))
  alldists$memberIdx <- 1:length(clusInfo$clusterNumbers)
  for(clusNo in clusNumbers){
    alldists[[paste0("c",clusNo)]] <- sapply(alldists$memberIdx,FUN = getAvgDistanceFromOtherClusterMembers,clusInfo,clusNo)
  }
  return(alldists)
}

performMemberSwap <- function(memberId1, memberId2,clusInfo){
  rightClusNo <- clusInfo$clusterNumbers[memberId1]
  leftClusNo <- clusInfo$clusterNumbers[memberId2]
  rightClusMemberIdx <- clusInfo$memberClusterMap[[rightClusNo]]
  leftClusMemberIdx <- clusInfo$memberClusterMap[[leftClusNo]]
  clusInfo$memberClusterMap[[rightClusNo]] <- c(memberId2,rightClusMemberIdx[which(rightClusMemberIdx!=memberId1)])
  clusInfo$memberClusterMap[[leftClusNo]] <- c(memberId1,leftClusMemberIdx[which(leftClusMemberIdx!=memberId2)])
  
  clusInfo$clusterNumbers[memberId1] <- leftClusNo
  clusInfo$clusterNumbers[memberId2] <- rightClusNo
  return(clusInfo)  
}

performMovementsToImproveClusters <- function(clusInfo,maxNumberOfIterations=NULL){
  if(is.null(maxNumberOfIterations)){
    numOfClus <- length(clusInfo$memberClusterMap)
    maxNumberOfIterations <- 5*numOfClus^numOfClus
    maxNumberOfIterations <- min(maxNumberOfIterations,1000000) # max 1M time which takes 4 mins
    maxNumberOfIterations <- max(maxNumberOfIterations,10000) #atleast 10k times since anyways if there isnt any improvement this will be stopped
  }
  clusterLabels <- 1:length(clusInfo$memberClusterMap)
  clusdists <- getAllMemberAvgDistancesFromGivenClusters(clusterLabels,clusInfo)
  avgDistAllClusterPrev <- mean(getAvgClusterDistanceForAllClusters(clusInfo))
  for(i in 1:maxNumberOfIterations){
    if(i%%1000==0){ # every 1000 iteration check if there is any improvement else stop this.
      avgDistAllClusterCurr <- mean(getAvgClusterDistanceForAllClusters(clusInfo))
      if(avgDistAllClusterCurr==avgDistAllClusterPrev) break
      else avgDistAllClusterPrev <- avgDistAllClusterCurr
    }
    #randomly pick to clusters to consider for swapping
    clustersToSwap <- sample(clusterLabels,size = 2)
    rightClusterNo <- clustersToSwap[1]
    leftClusterNo <- clustersToSwap[2]
    rightClusterMember <- clusInfo$memberClusterMap[[rightClusterNo]]
    leftClusterMember <- clusInfo$memberClusterMap[[leftClusterNo]]
    rightClusterMemberDistDecrease <- clusdists[rightClusterMember,(rightClusterNo+1)]-clusdists[rightClusterMember,(leftClusterNo+1)]
    leftClusterMemberDistDecrease <- clusdists[leftClusterMember,(leftClusterNo+1)]-clusdists[leftClusterMember,(rightClusterNo+1)]
    if(max(rightClusterMemberDistDecrease)+max(leftClusterMemberDistDecrease)>0){
      #swap will happen
      rightMemberToSwap <- rightClusterMember[which.max(rightClusterMemberDistDecrease)]
      leftMemberToSwap <- leftClusterMember[which.max(leftClusterMemberDistDecrease)]
      avgErrorBefore <- mean(getAvgClusterDistanceForAllClusters(clusInfo,clustersToSwap))
      tmp <- performMemberSwap(rightMemberToSwap,leftMemberToSwap,clusInfo)
      avgErrorAfter <- mean(getAvgClusterDistanceForAllClusters(tmp,clustersToSwap))
      if(avgErrorAfter<avgErrorBefore) clusInfo <- tmp
    }
  }
  return(clusInfo)
}

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