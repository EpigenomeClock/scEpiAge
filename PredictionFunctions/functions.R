##Functions needed for the prediction of epi-genetic age

#### Read COV files
readCovFiles <- function(inputFolder, backupInformation){
  sitesToConsiderFull = unique(c(backupInformation[,1],backupInformation[,2]))
  scFiles = list.files(inputFolder,full.names = T)
  scMeth = NULL
  for(f in scFiles){
    scInfo <- read.delim(f,as.is=T,header=F,colClasses = c("character","double","double","double","double","double"))
    
    scInfo[,1] = paste(scInfo[,1],scInfo[,2],sep=":")
    #Proper bismark file with all info.
    scInfo[,4] = scInfo[,4]/100
    scInfo[,5] = scInfo[,5]+scInfo[,6]
    scInfo = scInfo[,c(1,4)]
    
    scInfo = scInfo[which(scInfo[,1] %in% sitesToConsiderFull),]
    partsName = strsplit(f,"/")[[1]]
    colnames(scInfo) <- c("ID",paste("MethRate_", partsName[length(partsName)],sep = ""))
    
    if(is.null(scMeth)){
      scMeth = scInfo
    } else {
      scMeth = merge(scMeth, scInfo,by = "ID",all = T)
    }
  }
  
  rownames(scMeth) = scMeth$ID
  scMeth = scMeth[,-1]
  scMethMat = as.matrix(scMeth)
  rm(scMeth)
  return(scMethMat)
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
    predictionMatrix[sc,1] <- length(relMeth)
    replacementValues = NULL
    
    if(length(relMeth)<length(sitesToConsider)){
      #Get infomation for replacement values.
      replacementValues = backupInformation[which(backupInformation[,1]%in% names(which(is.na(methData_validation_sel[,sc])))),]
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
    
    if(is.null(replacementValues)){
      for(cpg in (1:(length(relMeth)))){
        scores = log(1-abs(expectedMethMat_sel[relMeth[cpg],]- methData_validation_sel[relMeth[cpg],sc]))
        ageProbability = ageProbability+scores;
      }
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
      
      if(!all(rownames(methData_validation_sel_t)==rownames(expectedMethMat_sel_t))){
        print("something wrong with the input matrices, after placing backup values.")
        stop();
      }
      
      relMeth = which(!is.na(methData_validation_sel_t[,sc]))
      
      predictionMatrix[sc,1] <- length(relMeth)
      #MLE
      for(cpg in (1:(length(relMeth)))){
        #predictionMatrix
        #Works for high methylation.
        scores = log(1-abs(expectedMethMat_sel_t[relMeth[cpg],]- methData_validation_sel_t[relMeth[cpg],sc]))
        ageProbability = ageProbability+scores;
      }
      rm(methData_validation_sel_t,expectedMethMat_sel_t)
    }
    
    predictionMatrix[sc,2] = floor(median(as.numeric(names(ageProbability)[which(ageProbability == max(ageProbability))])))
    
  }
  rownames(predictionMatrix) = colnames(methData_validation_sel)
  return(predictionMatrix)
}