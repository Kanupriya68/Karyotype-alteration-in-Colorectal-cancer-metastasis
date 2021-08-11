#Read Hartwig liver copy number data
met_liver <- read.csv("C:/Research/Research_project/tcga/HartwigLiver.csv",header = TRUE,sep = ",")
#Read Hartwig lung copy number data
met_lung <- read.csv("C:/Research/Research_project/tcga/HartwigLung.csv",header = TRUE,sep = ",")

#Combine metastatic data
met_lung <- subset(met_lung,select = c(4:30))
metastatic_cndata <- cbind(met_liver,met_lung)
metastatic_cndata <- subset(metastatic_cndata, select = -c(1:3))

#Read TCGA liver copy number data
pri_liver <- read.csv("C:/Research/Research_project/TCGAbinned/liverOutput.csv", header = TRUE,sep = ",")
pri_liver <- subset(pri_liver,select = -c(1:5))
#Read TCGA lung copy number data
pri_lung <- read.csv("C:/Research/Research_project/TCGAbinned/Lungoutput.csv", header = TRUE,sep = ",")
pri_lung <- subset(pri_lung,select = -c(1:5))

#Combine primary data
primary_cndata <- cbind(pri_liver,pri_lung)

#Transpose metastatic data
t_meta_cnData <- t(metastatic_cndata)

#Name columns as per bin number
metaBinCols <- colnames(t_meta_cnData)

#Create empty list for gain and loss to append later
metaGainList <- list()
metaLossList <- list()

#Create for loop to calculate gain and loss proportions bin wise
for(i in metaBinCols){
  gain <- length(which(t_meta_cnData[,i] > 2))
  loss <- length(which(t_meta_cnData[,i] < 2))
  total <- length(t_meta_cnData[,i])
  
  gainRatio <- gain/total
  lossRatio <- loss/total
  metaGainList <- append(metaGainList, gainRatio)
  metaLossList <- append(metaLossList, lossRatio)
  
}

#Unlist the object created from the loop to extract the array
metaGainList <- unlist(metaGainList)
metaLossList <- unlist(metaLossList)



####################################################

#Transpose
t_primary_cnData <- t(primary_cndata)
#Name columns
PimaryBinCols <- colnames(t_primary_cnData)

#Create empty list for append
primaryGainList <- list()
primaryLossList <- list()

#Create for loop to calculate gain and loss proportions
for(i in PimaryBinCols){
  gain <- length(which(t_primary_cnData[,i] > 2))
  loss <- length(which(t_primary_cnData[,i] < 2))
  total <- length(t_primary_cnData[,i])
  
  gainRatio <- gain/total
  lossRatio <- loss/total
  primaryGainList <- append(primaryGainList, gainRatio)
  primaryLossList <- append(primaryLossList, lossRatio)
  
}

#Extract array from object created from above loop
primaryGainList <- unlist(primaryGainList)
primaryLossList <- unlist(primaryLossList)

#Perform t test to compare gains in primary and metastatic
t.test(primaryGainList,metaGainList)

#Perform t test to compare losses in primary and metastatic
t.test(primaryLossList,metaLossList)












