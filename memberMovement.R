# Functions to be used later ----------------------
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


# Current cluster balancing code --------------------

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

