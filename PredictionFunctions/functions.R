##Functions needed for the prediction of epi-genetic age

#### Read COV files
readCovFiles <- function(inputFolder, backupInformation, sampleSelection=NULL){
  sitesToConsiderFull = unique(c(backupInformation[,1],backupInformation[,2]))
  scFiles = list.files(inputFolder,full.names = T,recursive = T)
  scFiles = scFiles[grep(scFiles,pattern = ".cov$|.cov.gz$")]
  
  if(!is.null(sampleSelection)){
    scFiles=scFiles[grep(scFiles,pattern = sampleSelection)]
  }
  
 
  if(length(scFiles)==0){
    print("Error no files to be read in.")
  }
  
  
  options(warn = 2)
  scMeth = NULL
  for(f in scFiles){
    print(f)
    scInfo <- read.delim(f,as.is=T,header=F,colClasses = c("character","double","double","double","double","double"))
    
    scInfo[,1] = paste(scInfo[,1],scInfo[,2],sep=":")
    
    ##subset to interesting sites.
    scInfo = scInfo[which(scInfo[,1] %in% sitesToConsiderFull),]
    
    #Make into 0-1
    scInfo[,4] = scInfo[,4]/100
    
    #Read-depth.
    scInfo[,5] = scInfo[,5]+scInfo[,6]
    print(paste("Average read-depth clock sites:",mean(scInfo[,5])))
    
    #Subset to methylation and row id only.
    scInfo = scInfo[,c(1,4)]
    
    partsName = strsplit(f,"/")[[1]]
    colnames(scInfo) <- c("ID",paste("MethRate_", partsName[length(partsName)],sep = ""))
    
    if(is.null(scMeth)){
      scMeth = scInfo
    } else {
      scMeth = merge(scMeth, scInfo,by = "ID",all = T)
    }
  }
  options(warn = 1)
  
  #Force (missing) primary sites to be added as NA. Code is expecting these are at least present.
  toAdd = unique(backupInformation[,1])
  toAdd = toAdd[-which(toAdd %in% scMeth$ID)]
  if(length(toAdd)!=0){
    extraMatDat = matrix(data = NA,nrow = length(toAdd),ncol = ncol(scMeth))
    colnames(extraMatDat) = colnames(scMeth)
    extraMatDat[,which(colnames(extraMatDat)=="ID")] = toAdd
    extraMatDat = as.data.frame(extraMatDat)
    for(i in 2:ncol(extraMatDat)){
      extraMatDat[,i] = as.numeric(extraMatDat[,i])
    }
    scMeth = rbind(scMeth,extraMatDat)
  }
  
  rownames(scMeth) = scMeth$ID
  scMeth = scMeth[,-which(colnames(scMeth)=="ID")]
  scMethMat = as.matrix(scMeth)
  rm(scMeth)
  return(scMethMat)
}

#### Read subsetted COV files
readScCovFiles <- function(inputFolder, backupInformation){
  sitesToConsiderFull = unique(c(backupInformation[,1],backupInformation[,2]))
  scFiles = list.files(inputFolder,full.names = T,recursive = T)
  scFiles = scFiles[grep(scFiles,pattern = ".cov$|.cov.gz$")]
  
  if(length(scFiles)==0){
    print("Error no files to be read in.")
  }
  
  options(warn = 2)
  scMeth = NULL
  for(f in scFiles){
    print(f)
    scInfo <- read.delim(f,as.is=T,header=F,colClasses = c("character","double","double"))
    
    scInfo[,1] = paste(scInfo[,1],scInfo[,2],sep=":")
    
    ##subset to interesting sites.
    scInfo = scInfo[which(scInfo[,1] %in% sitesToConsiderFull),]
    
    #Make into 0-1
    #scInfo[,4] = scInfo[,4]/100
    
    print("Average read-depth clock sites: 1")
    
    #Subset to methylation and row id only.
    scInfo = scInfo[,c(1,3)]
    
    partsName = strsplit(f,"/")[[1]]
    colnames(scInfo) <- c("ID",paste("MethRate_", partsName[length(partsName)],sep = ""))
    
    if(is.null(scMeth)){
      scMeth = scInfo
    } else {
      scMeth = merge(scMeth, scInfo,by = "ID",all = T)
    }
  }
  options(warn = 1)
  
  #Force (missing) primary sites to be added as NA. Code is expecting these are at least present.
  toAdd = unique(backupInformation[,1])
  toAdd = toAdd[-which(toAdd %in% scMeth$ID)]
  if(length(toAdd)!=0){
    extraMatDat = matrix(data = NA,nrow = length(toAdd),ncol = ncol(scMeth))
    colnames(extraMatDat) = colnames(scMeth)
    extraMatDat[,which(colnames(extraMatDat)=="ID")] = toAdd
    extraMatDat = as.data.frame(extraMatDat)
    for(i in 2:ncol(extraMatDat)){
      extraMatDat[,i] = as.numeric(extraMatDat[,i])
    }
    scMeth = rbind(scMeth,extraMatDat)
  }
  
  rownames(scMeth) = scMeth$ID
  scMeth = scMeth[,-which(colnames(scMeth)=="ID")]
  scMethMat = as.matrix(scMeth)
  rm(scMeth)
  return(scMethMat)
}

processArrayData <- function(inputFile, backupInformation, sampleSelection=NULL){
  sitesToConsiderFull = unique(c(backupInformation[,1],backupInformation[,2]))
  methMat = read.delim(inputFile,as.is=T,check.names = F,sep=",")
  methMat$chr = gsub(pattern = "chr",replacement = "",methMat$chr)
  methMat["ID"] = paste(methMat$chr,methMat$start,sep = ":")
  methMat = methMat[,c(ncol(methMat),5:(ncol(methMat)-1))]
  methMat= methMat[which(methMat$ID%in%sitesToConsiderFull),]
  
  #Force (missing) primary sites to be added as NA. Code is expecting these are at least present.
  toAdd = unique(backupInformation[,1])
  toAdd = toAdd[-which(toAdd %in% methMat$ID)]
  if(length(toAdd)!=0){
    extraMatDat = matrix(data = NA,nrow = length(toAdd),ncol = ncol(methMat))
    colnames(extraMatDat) = colnames(methMat)
    extraMatDat[,which(colnames(extraMatDat)=="ID")] = toAdd
    extraMatDat = as.data.frame(extraMatDat)
    for(i in 2:ncol(extraMatDat)){
      extraMatDat[,i] = as.numeric(extraMatDat[,i])
    }
    methMat = rbind(methMat,extraMatDat)
  }
  
  rownames(methMat) = methMat$ID
  methMat = methMat[,-which(colnames(extraMatDat)=="ID")]
  methMat = as.matrix(methMat)
  return(methMat)
}

##Prediction
predictAges <- function(scMethMat, backupInformation, expectedMethMat){
  ##Test how the prediction would work on bulk [Not assuming 0-1 values only.]
  sitesToConsider = unique(backupInformation[,1])
  methData_validation_sel = scMethMat[which(rownames(scMethMat) %in% (sitesToConsider)),]
  methData_validation_sel = methData_validation_sel[order(rownames(methData_validation_sel)),]
  expectedMethMat_sel = expectedMethMat[which(rownames(expectedMethMat) %in% rownames(methData_validation_sel)),]
  if(!all(rownames(methData_validation_sel)==rownames(expectedMethMat_sel))){
    print("something wrong with the input matrices")
    stop();
  }
  
  #To store prediction information.
  predictionMatrix <- matrix(NA,ncol=2, nrow=ncol(methData_validation_sel))
  colnames(predictionMatrix) <- c("sitesUsed","predictedAge")
  ##Quick correlation and MLE based method (for speed up). We do 2 MLE methods one weighted one unweighted, for correlations we can't use read-depth.
  for(sc in 1:ncol(methData_validation_sel)){
    ##statics.
    ageProbability = rep(0,ncol(expectedMethMat_sel))
    maxDistance = 0
    
    #Actual data.
    relMeth = which(!is.na(methData_validation_sel[,sc]))
    replacementValues = NULL
    
    ##Try and add backup info.
    if(length(relMeth)<length(sitesToConsider)){
      #Get infomation for replacement values.
      replacementValues = backupInformation[which(backupInformation[,1]%in% rownames(methData_validation_sel)[(which(is.na(methData_validation_sel[,sc])))]),]
      replacementValues = replacementValues[which(replacementValues[,2]!="-"),]
      if(!is.null(nrow(replacementValues))) {
        if(nrow(replacementValues)==0){
          replacementValues = NULL
        } else {
          
          ##select interesting sites for this sample.
          replacementValues = replacementValues[which(replacementValues[,2] %in% rownames(scMethMat)),]
          
          #S1 drop backup sites that are missing.
          replacementValues = replacementValues[which(!(replacementValues[,2] %in% names(which(is.na(scMethMat[which(rownames(scMethMat)%in%replacementValues[,2]),sc]))))),]
          #This one can be length 0 (non-safe).
          if(is.na(replacementValues[1,1])){
            # print(paste("s2",sc,replacementValues[1]))
            replacementValues = NULL
          } else if(is.null(dim(replacementValues))){
            #only one entry so per definition unique.
            replacementValues = replacementValues[2]
          } else {
            #Select 1 site per missing value.
            replacementValues = replacementValues[which(!duplicated(replacementValues[,1])),]
            #This one can't be length 0 (safe).
            if(is.null(dim(replacementValues))){
              replacementValues = replacementValues[2]
            } else {
              #Making sure each backup is only introduced once.
              replacementValues = unique(replacementValues[,2])
            }
          }
        }
      }
      
      ##
    }

    if((length(relMeth)+length(replacementValues))<5){
      print(c(rownames(methData_validation_sel)[relMeth],replacementValues))
      next();
    }
    
    methData_validation_sel_t = NULL
    expectedMethMat_sel_t = NULL    
    if((is.null(replacementValues))){
      methData_validation_sel_t = methData_validation_sel[relMeth,sc]
      expectedMethMat_sel_t = expectedMethMat_sel[relMeth,]
    } else {
      ##Need to add the backup sites here.
      methData_validation_sel_t = rbind(methData_validation_sel,scMethMat[which(rownames(scMethMat)%in%replacementValues),])
      expectedMethMat_sel_t = rbind(expectedMethMat_sel,expectedMethMat[which(rownames(expectedMethMat)%in%replacementValues),])
      
      ##Fix issue with adding one row (no row names in that case.)
      if(length(replacementValues)==1){
        rownames(methData_validation_sel_t)[nrow(methData_validation_sel_t)] = as.character(replacementValues)
        rownames(expectedMethMat_sel_t)[nrow(expectedMethMat_sel_t)] = as.character(replacementValues)
      }
      methData_validation_sel_t = methData_validation_sel_t[order(rownames(methData_validation_sel_t)),]
      expectedMethMat_sel_t = expectedMethMat_sel_t[order(rownames(expectedMethMat_sel_t)),]

      relMeth = which(!is.na(methData_validation_sel_t[,sc]))
      
      methData_validation_sel_t = methData_validation_sel_t[relMeth,sc]
      expectedMethMat_sel_t = expectedMethMat_sel_t[relMeth,]
    }
    
    if(!all(rownames(methData_validation_sel_t)==rownames(expectedMethMat_sel_t))){
      print("something wrong with the input matrices.")
      stop();
    }
    
    predictionMatrix[sc,1] <- length(relMeth)
    ageProbability = colSums(log(1-abs(expectedMethMat_sel_t - methData_validation_sel_t)))
    plot(ageProbability)
    predictionMatrix[sc,2] = floor(median(as.numeric(names(ageProbability)[which(ageProbability == max(ageProbability))])))
    rm(methData_validation_sel_t,expectedMethMat_sel_t)
    
  }
  rownames(predictionMatrix) = colnames(methData_validation_sel)
  return(predictionMatrix)
}

##Prediction and simulated expected age prediction on the same sites.
predictAgesAndCalculateExpectedGivenAge <- function(scMethMat, backupInformation, expectedMethMat, expectedAges, nSimulations, plotHist=T){
  #scMethMat = inputMethMatrix; expectedMethMat = expectedMethMatrix; expectedAges = ageInfo;
  
  sitesToConsider = unique(backupInformation[,1])
  methData_validation_sel = scMethMat[which(rownames(scMethMat) %in% (sitesToConsider)),]
  methData_validation_sel = methData_validation_sel[order(rownames(methData_validation_sel)),]
  expectedMethMat_sel = expectedMethMat[which(rownames(expectedMethMat) %in% rownames(methData_validation_sel)),]
  
  ##Order colnames and first column
  methData_validation_sel = methData_validation_sel[,order(colnames(methData_validation_sel))]
  expectedAges = expectedAges[order(expectedAges[,1]),]
  expectedAges = expectedAges[which(expectedAges[,1] %in% colnames(methData_validation_sel)),]
  
  if(!all(rownames(methData_validation_sel)==rownames(expectedMethMat_sel))){
    print("something wrong with the input matrices")
    stop();
  }
  
  if(!all(colnames(scMethMat) == expectedAges[,1])){
    print("something wrong with the age information provided for the input methylation matrix.")
    stop();
  }
  
  #To store prediction information.
  predictionMatrix <- matrix(NA,ncol=9, nrow=ncol(methData_validation_sel))
  colnames(predictionMatrix) <- c("actualAge","sitesUsed","readDepth","predictedAge","ageDeviation","medianRandomAge","IQR","IQR_Ratio","FDR")
  ##Quick correlation and MLE based method (for speed up). We do 2 MLE methods one weighted one unweighted, for correlations we can't use read-depth.
  for(sc in 1:ncol(methData_validation_sel)){
    print(expectedAges[sc,])
    ##statics.
    ageProbability = rep(0,ncol(expectedMethMat_sel))
    maxDistance = 0
    
    #Actual data.
    relMeth = which(!is.na(methData_validation_sel[,sc]))
    replacementValues = NULL
    
    ##Try and add backup info.
    if(length(relMeth)<length(sitesToConsider)){
      #Get infomation for replacement values.
      replacementValues = backupInformation[which(backupInformation[,1]%in% rownames(methData_validation_sel)[(which(is.na(methData_validation_sel[,sc])))]),]
      replacementValues = replacementValues[which(replacementValues[,2]!="-"),]
      if(!is.null(nrow(replacementValues))) {
        if(nrow(replacementValues)==0){
          replacementValues = NULL
        } else {
          
          ##select interesting sites for this sample.
          replacementValues = replacementValues[which(replacementValues[,2] %in% rownames(scMethMat)),]
          
          #S1 drop backup sites that are missing.
          replacementValues = replacementValues[which(!(replacementValues[,2] %in% names(which(is.na(scMethMat[which(rownames(scMethMat)%in%replacementValues[,2]),sc]))))),]
          #This one can be length 0 (non-safe).
          if(is.na(replacementValues[1,1])){
            # print(paste("s2",sc,replacementValues[1]))
            replacementValues = NULL
          } else if(is.null(dim(replacementValues))){
            #only one entry so per definition unique.
            replacementValues = replacementValues[2]
          } else {
            #Select 1 site per missing value.
            replacementValues = replacementValues[which(!duplicated(replacementValues[,1])),]
            #This one can't be length 0 (safe).
            if(is.null(dim(replacementValues))){
              replacementValues = replacementValues[2]
            } else {
              #Making sure each backup is only introduced once.
              replacementValues = unique(replacementValues[,2])
            }
          }
        }
      }
      
      ##
    }
    
    if((length(relMeth)+length(replacementValues))<5){
      next();
    }
    
    methData_validation_sel_t = NULL
    expectedMethMat_sel_t = NULL    
    if((is.null(replacementValues))){
      methData_validation_sel_t = methData_validation_sel[relMeth,sc]
      expectedMethMat_sel_t = expectedMethMat_sel[relMeth,]
    } else {
      ##Need to add the backup sites here.
      methData_validation_sel_t = rbind(methData_validation_sel,scMethMat[which(rownames(scMethMat)%in%replacementValues),])
      expectedMethMat_sel_t = rbind(expectedMethMat_sel,expectedMethMat[which(rownames(expectedMethMat)%in%replacementValues),])
      
      ##Fix issue with adding one row (no row names in that case.)
      if(length(replacementValues)==1){
        rownames(methData_validation_sel_t)[nrow(methData_validation_sel_t)] = as.character(replacementValues)
        rownames(expectedMethMat_sel_t)[nrow(expectedMethMat_sel_t)] = as.character(replacementValues)
      }
      methData_validation_sel_t = methData_validation_sel_t[order(rownames(methData_validation_sel_t)),]
      expectedMethMat_sel_t = expectedMethMat_sel_t[order(rownames(expectedMethMat_sel_t)),]
      
      relMeth = which(!is.na(methData_validation_sel_t[,sc]))
      
      methData_validation_sel_t = methData_validation_sel_t[relMeth,sc]
      expectedMethMat_sel_t = expectedMethMat_sel_t[relMeth,]
    }
    
    if(!all(rownames(methData_validation_sel_t)==rownames(expectedMethMat_sel_t))){
      print("something wrong with the input matrices.")
      stop();
    }
    
    predictionMatrix[sc,1] <- expectedAges[sc,2]
    predictionMatrix[sc,2] <- length(relMeth)
    predictionMatrix[sc,3] <- expectedAges[sc,3]
    ageProbability = colSums(log(1-abs(expectedMethMat_sel_t - methData_validation_sel_t)))
    predictionMatrix[sc,4] = floor(median(as.numeric(names(ageProbability)[which(ageProbability == max(ageProbability))])))
    predictionMatrix[sc,5] = abs(predictionMatrix[sc,1]-predictionMatrix[sc,4])
    
    ##Start predictions on random methylation values, on the right sides and forming the expected bulk methylation profile for a given age.
    colOfInterest = which(colnames(expectedMethMat)==expectedAges[sc,2])
    
    expectedMethMat_sel_r = as.matrix(expectedMethMat_sel_t[,colOfInterest],ncol=1)
    rownames(expectedMethMat_sel_r) = rownames(expectedMethMat_sel_t)
    
    expectedRandomScData = matrix(NA,ncol=nSimulations,nrow=nrow(expectedMethMat_sel_r))
    rownames(expectedRandomScData)= rownames(expectedMethMat_sel_r)
    
    if(!all(rownames(expectedRandomScData) == rownames(expectedMethMat_sel_r))){
      #Rownames don't match
      stop();
    }
    
    for(x in 1:expectedAges[sc,3]){
      for(rId in 1:nrow(expectedRandomScData)){
        methylationvalues <- runif(nSimulations, 0, 1)
        methVexp = expectedMethMat_sel_r[rId,1]
        methVexp = quantile(methylationvalues,(1-methVexp))
        if(all(is.na(expectedRandomScData[rId,]))){
          methylationvalues[which(methylationvalues<=methVexp)]=0
          methylationvalues[which(methylationvalues>methVexp)]=1
          expectedRandomScData[rId,] = methylationvalues 
        } else {
          if((mean(expectedRandomScData[rId,])/(x-1))>=expectedMethMat_sel_r[rId,1]){
            methylationvalues[which(methylationvalues<=methVexp)]=0
            methylationvalues[which(methylationvalues>methVexp)]=1
          } else {
            methylationvalues[which(methylationvalues<methVexp)]=0
            methylationvalues[which(methylationvalues>=methVexp)]=1
          }
          expectedRandomScData[rId,] = expectedRandomScData[rId,]+methylationvalues 
        }
      }
    }
    expectedRandomScData = expectedRandomScData/expectedAges[sc,3]
    #print(mean(rowMeans(expectedRandomScData)-expectedMethMat_sel_r[,1]))
    
    ##Age predict.
    predictedAgesRandom = NULL
    
    for(cId in 1:nSimulations){
      ageProbabilityR = (1-abs(expectedMethMat_sel_t - expectedRandomScData[,cId]))
      ageProbabilityR = colSums(log(1-abs(expectedMethMat_sel_t - expectedRandomScData[,cId])))
      predictedAgesRandom = c(predictedAgesRandom,as.numeric(colnames(expectedMethMat_sel_t)[which(ageProbabilityR == max(ageProbabilityR))]))
    }
    
    if(plotHist){
      hist(predictedAgesRandom,xlim=c(min(as.numeric(colnames(expectedMethMat))),max(as.numeric(colnames(expectedMethMat)))))
      abline(v=predictionMatrix[sc,4],col="red")
      abline(v=median(predictedAgesRandom),col="blue")
    }
    
    predictionMatrix[sc,6] =  median(predictedAgesRandom)
    predictionMatrix[sc,7] =  IQR(predictedAgesRandom)

    predictionMatrix[sc,8] =  (predictionMatrix[sc,4]-predictionMatrix[sc,6])/predictionMatrix[sc,7]
    if(is.na(predictionMatrix[sc,8])){
      predictionMatrix[sc,8] = NA
      predictionMatrix[sc,9] = 1
    } else if(predictionMatrix[sc,8]>0){
      predictionMatrix[sc,9] = length(which(predictedAgesRandom>=predictionMatrix[sc,4]))/nSimulations
    } else {
      predictionMatrix[sc,9] = length(which(predictedAgesRandom<=predictionMatrix[sc,4]))/nSimulations
    }
    #Making sure we don't under-estimate FDR, rounding up to the minimal possible at this level of permutations.
    if(predictionMatrix[sc,9]==0){
      predictionMatrix[sc,9] = 0.9/nSimulations
    }
    rm(methData_validation_sel_t,expectedMethMat_sel_t,expectedRandomScData)
  }
  
  rownames(predictionMatrix) = colnames(methData_validation_sel)
  return(predictionMatrix)
}

##Simulating and dropping sites X times
dropSites <- function(scMethMat, nRep = 20, minSites=5, maxSites = 750 ) {
  BabrahamSimExtended = NULL
  
  for(sampleN in 1:ncol(scMethMat)){
    sampleName = colnames(scMethMat)[sampleN]
    
    #Here we go to only the observed sites to start.
    toSampleFrom = which(!is.na(scMethMat[,sampleN]))
    
    if(length(toSampleFrom)<maxSites){
      maxSites = length(toSampleFrom)  
    }
    
    steps = (maxSites - minSites+1)
    #print(paste(sampleName,maxSites))
    simExtended = matrix(NA,nrow=nrow(scMethMat),ncol=(steps*nRep))
    colCounter = 0;
    colnames(simExtended) = rep("",ncol(simExtended))
    rownames(simExtended) = rownames(scMethMat)
    for(siteN in c(minSites:maxSites)){
      if(siteN!=maxSites){
        for(repN in 1:nRep){
          colCounter = colCounter+1
          sampleInfo = as.numeric(sample(toSampleFrom,siteN))
          simExtended[sampleInfo,colCounter] = scMethMat[sampleInfo,sampleN]
          colnames(simExtended)[colCounter] = paste(sampleName,"#Observed:",siteN,"#Replication:",repN,sep="")
        }
      } else {
        colCounter = colCounter+1
        simExtended[toSampleFrom,colCounter] = scMethMat[toSampleFrom,sampleN]
        colnames(simExtended)[colCounter] = paste(sampleName,"#Observed:",siteN,"#All",sep="")
        ##All sites.
      }
    }
    BabrahamSimExtended = cbind(BabrahamSimExtended,simExtended)
  }
  return(BabrahamSimExtended)
}
