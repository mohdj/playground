source("CollaborationAI/rscripts/generateSet.R")
source("CollaborationAI/rscripts/memberMovement.R")

# Big data clustering function -----------------
prepareDataForKmeans <- function(dat,clusInfo){
  newdat <- data.frame(tmpid=1:nrow(dat))
  weight <- clusInfo$bigDataClusWeight
  weight <- weight[weight>0] # take only cols with similarity criteria
  dat <- dat[,names(weight)]
  weight <- weight/sum(clusInfo$bigDataClusWeight)
  colTypes <- colTypes[names(weight)]
  listSeparator <- ","
  
  singleChoiceidx <- names(colTypes)[which(colTypes=="SingleChoice")]
  for(i in singleChoiceidx){
    tmp <- dat[[i]]
    tmp[is.na(tmp)] <- "" # make NAs as empty strings
    tmp <- paste0("x",tmp) # this is done to convert all numeric to string
    tmp <- gsub(pattern = " ",replacement = "",x = dat[[i]])
    askdoc <- VectorSource(tmp)
    askdoc <- VCorpus(askdoc)
    askdoc <- tm_map(askdoc, tolower)   # *Converting to lowercase:*    
    askdoc <- tm_map(askdoc, PlainTextDocument)  
    dtm <- DocumentTermMatrix(askdoc)  
    dtm <- as.data.frame(as.matrix(dtm))
    row.names(dtm) <- NULL
    names(dtm) <- paste0(names(dtm),"__",i)
    dtm <- dtm*0.5 # since only two cols can have dissimilarity at max, so max distance is 1
    dtm <- dtm*weight[[i]] # multiply it by weight
    newdat <- cbind(newdat,dtm)
  }
  # List from a set
  listfromsetidx <- names(colTypes)[which(colTypes=="MultiChoice")]
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
    row.names(dtm) <- NULL
    names(dtm) <- paste0(names(dtm),"__",i)
    dtm <- dtm*(1/ncol(dtm)) # distributing the weight between them evenly
    dtm <- dtm*weight[[i]] # multiply it by weight
    newdat <- cbind(newdat,dtm)
  }
  # Plain text
  plainTextidx <- names(colTypes)[which(colTypes=="PlainText")]
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
    dtm <- dtm*(1/ncol(dtm)) # distributing the weight between them evenly
    dtm <- dtm*weight[[i]] # multiply it by weight
    newdat <- cbind(newdat,dtm)
  }
  newdat$tmpid <- NULL
  return(newdat)
}

performBigDataClustering <- function(dat,clusInfo){
  out <- list()
  dtm <- prepareDataForKmeans(dat,clusInfo)
  clusInfo$columnTermFrequencyColumnMapping <- getColumnTermFrequencyColumnMapping(dtm,clusInfo)
  clus.fit <- clara(dtm,k=bigDataClusNumberOfGroup,metric = "manhattan",stand = F,correct.d = T)
  clusInfo$clusterNumbers <- clus.fit$clustering
  bigDataClusNumbers <- unique(clusInfo$clusterNumbers)
  currentClusCounter <- 0
  newdat <- list()
  for(bdClusNo in bigDataClusNumbers){
    tmp <- dat[which(clusInfo$clusterNumbers==bdClusNo),]
    localClusInfo <- createGroups(dat = tmp,weight = clusInfo$weight,colTypes = clusInfo$colTypes,groupSize = clusInfo$idealClusterSize)
    clusterNumbers <- localClusInfo$clusterNumbers + currentClusCounter
    tmp$cluster <- clusterNumbers
    newdat[[bdClusNo]] <- tmp
    currentClusCounter <- currentClusCounter + localClusInfo$numberOfGroup
    print(bdClusNo)
  }
  newdat <- rbind_all(newdat)
  
  clusInfo$clusterNumbers <- newdat$cluster
  clusterIndices <- list()
  for(clusNo in unique(clusInfo$clusterNumbers)){
    clusterIndices[[clusNo]] <- which(clusInfo$clusterNumbers==clusNo)
  }
  clusInfo[["memberClusterMap"]] <- clusterIndices
  
  allTermColName <- character()
  for(nm in names(clusInfo$columnTermFrequencyColumnMapping)){
    allTermColName <- c(allTermColName,clusInfo$columnTermFrequencyColumnMapping[[nm]])  
  }
  tmp <- dtm[,allTermColName]
  tmp <- as.data.frame(tmp>0)
  tmp <- lapply(tmp,as.numeric)
  tmp <- as.data.frame(tmp)
  tmp <- cbind(newdat,tmp)
  perf <- getClusteringPerformance(tmp,clusInfo)
  out[["ClusteringResult"]] <- newdat
  out[["ClusteringPerformance"]] <- perf
  return(out)
}

getSimulatedBigData <- function(N,sampleSet,numOfSingleChoiceCols,numOfMultiChoiceCols,possibleNumOfItemsPerValue){
  # Generate single choice columns
  sampset <- eval(parse(text=paste0("c('",paste(paste0("a",sampleSet),collapse="','"),"')")))
  codestr <- paste0("data.frame(",paste(paste0("x",1:numOfSingleChoiceCols,"=sample(sampset,N,T)"),collapse=","),",stringsAsFactors=F)")
  dat <- eval(parse(text=codestr))
  #generate multi choice columns
  numOfCols <- numOfSingleChoiceCols+numOfMultiChoiceCols
  for(j in (numOfSingleChoiceCols+1):numOfCols){
    newcol <- list()
    for(i in 1:N){
      newcol[[i]] <- paste(paste0("a",sample(x = sampleSet,size = sample(x = possibleNumOfItemsPerValue,size = 1))),collapse=",")
    }
    newcol <- sapply(newcol,c)  
    dat[[paste0("x",j)]] <- newcol
  }  
  return(dat)
}

# Testing Validation of Code ----------------
sampleSet <- c("one","two","three","four","five","six")
dat <- getSimulatedBigData(N=1000,sampleSet=sampleSet,numOfSingleChoiceCols=3,numOfMultiChoiceCols=1,
                           possibleNumOfItemsPerValue=1:3)
bigDataClusWeight <- c(x1=5,x2=3,x4=2) # Final
bigDataClusNumberOfGroup <- 10
weight <- c(x1=4,x2=3,x3=-5,x4=2)
groupSize <- 7
colTypes <- c(x1="SingleChoice",x2="SingleChoice",x3="SingleChoice",x4="MultiChoice")
clusInfo <- list(bigDataClusWeight=bigDataClusWeight,bigDataClusNumberOfGroup=bigDataClusNumberOfGroup,
                 colTypes=colTypes,idealClusterSize=groupSize,datColNames=names(dat),weight=weight)

out <- performBigDataClustering(dat,clusInfo)

# Whats works and Next Steps -------------

# Whats works
# - Kmeans beatifully creates excellant similarity groups (single choice, multiple choice) also
# - taking into consideration relative weights - thats superb
# - ability to use kmeans with manhattan distance with factors, text was excellant, by converting 
# - everything to term frequency and than normalising the weights with 0.5 for single choice, 1/n (n=setsize)
# - for multi choice was a mind blowing idea and its beatifully working

# Next Step
# - Currently there are some violations of separation constraint
# - also there are possibility to get better solution by trying some movement across different clusters
# - may not be too much potential but still
# - Idea is to identify cluster & cluster members which are violating separation constraints & other high 
# - priority constraints and find the cluster with which it can be swapped (swapping is easier)