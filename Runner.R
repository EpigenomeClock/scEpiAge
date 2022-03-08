###### Settings ###### 
setwd("C:/Documents and Settings/admm414r/Documents/GitHub/AgingClock_v2/")
##Tissue, options are: liver, lung, blood.
tissue="liver"

##Folder with COV files.
#folderInput = "Gravina_2016/Bulk/"
folderInput = "Gravina_2016/SingleCell/"

## output file with aging information.
#outputFileName="./preditionGravinaBulk.txt"
outputFileName="./preditionGravinaSc.txt"
###################### 

###### Source ########
##source prediction functions.
source("./PredictionFunctions/functions.R")
###################### 

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

inputMethMatrix = readCovFiles(folderInput,backupInformation);


predictionOut = predictAges(inputMethMatrix, backupInformation, expectedMethMatrix);

write.table(predictionOut,outputFileName,sep="\t",quote=F)
