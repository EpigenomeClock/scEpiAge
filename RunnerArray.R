###### Settings ###### 
#setwd("C:/Documents and Settings/admm414r/Documents/GitHub/AgingClock_v2/")
setwd("C:/Documents and Settings/Administrator/Documents/GitHub/AgingClock_v2/")
##Tissue, options are: liver, lung, blood.
tissue="liver"

##File with info
fileInput = "./array_dummy_data.csv"

## output file with aging information.
outputFileName="./dummyDataOut.txt"

##If there are known age calculate confidence on the difference observed.
nSimulations = 5

##Age information
ageInfo = NULL

## ouput file name with extended information.
outputFileNameExtended = NA

plotHist=T

###################### 

## Actual code form here:
calcExtendedStats = F
if(!is.na(outputFileNameExtended) & !is.null(ageInfo)){
  calcExtendedStats=T
} else {
  print("Not given age information or output filename for the extended analysis")
}

###### Source ########
##source prediction functions.
source("./PredictionFunctions/functions.R")
###################### 

##### Load data ######
if(tissue=="liver"){
  expectedMethMatrix <- read.delim("./ExpectedMethylationMatrices/ExpectedMethMat_Liver.tsv",as.is=T,row.names=1,check.names = F)
  backupInformation <- read.delim("./SiteInformation/clockSites_Liver_BackUp.txt",as.is=T)
} else if(tissue=="lung"){
  expectedMethMatrix <- read.delim("./ExpectedMethylationMatrices/ExpectedMethMat_Lung.tsv",as.is=T,row.names=1,check.names = F)
  backupInformation <- read.delim("./SiteInformation/clockSites_Lung_BackUp.txt",as.is=T)
} else if(tissue=="blood"){
  expectedMethMatrix <- read.delim("./ExpectedMethylationMatrices/ExpectedMethMat_Blood.tsv",as.is=T,row.names=1,check.names = F)
  backupInformation <- read.delim("./SiteInformation/clockSites_Blood_BackUp.txt",as.is=T)
} else{
  print("No valid tissue selected")
  stop();
}


##Predict.
inputMethMatrix = processArrayData(fileInput,backupInformation);

if(calcExtendedStats){
  predictionOutVersusExpected = predictAgesAndCalculateExpectedGivenAge(inputMethMatrix, backupInformation, expectedMethMatrix, ageInfo, nSimulations, plotHist);
  write.table(predictionOutVersusExpected,outputFileNameExtended,sep="\t",quote=F)
} else {
  predictionOut = predictAges(inputMethMatrix, backupInformation, expectedMethMatrix);
  write.table(predictionOut,outputFileName,sep="\t",quote=F)  
}

