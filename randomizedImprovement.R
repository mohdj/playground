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

#  clusInfo <- readRDS("CollaborationAI/clusinfotmp.rdata")
# # clusInfo <- performMovementsToImproveClusters(clusInfo)
# 
# 
# system.time({kk <- performMovementsToImproveClusters(clusInfo,maxNumberOfIterations = 100000)})
