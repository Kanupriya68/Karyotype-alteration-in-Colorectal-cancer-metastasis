##Load libraries
y<-c("openxlsx", "tidyverse", "survminer", "ComplexHeatmap", "RColorBrewer", "ggbeeswarm", "MASS", "grid", "gtable", "gsubfn", "proto", "circlize", "xlsx", "RColorBrewer", "ComplexHeatmap", "grid", "data.table", "dplyr", "survminer", "ggpubr", "magrittr", "survival", "openxlsx", "reshape2", "reshape", "copynumber", "BiocGenerics", "parallel", "ggplot2", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
lapply(x, require, character.only = TRUE)
library(methods)

#pdf("Fig1A.pdf", width=10, height=10)

#Read metadata and create sampleID matrix
biopsySite <- read.csv(file = "C:/Research/Research_project/biopsy_metadata.csv", header = TRUE, sep =";")

#Read csv file containing patient sample copynumbers and chromosome info
chr_data <- read.csv(file = "C:/Research/Research_project/test/df.csv", header = TRUE, sep =",")
#Remove first column 'X'
chr_data <- subset(chr_data, select= -c(1,2,3))

#Transpose copynumber dataframe
transposed_chr_data <- t(chr_data)
transposed_chr_data <- as.data.frame(transposed_chr_data)

#Name the columns as per chromosome names
colnames_chr_data <- colnames(chr_data)
transposed_chr_data$sampleid <- colnames_chr_data
transposed_chr_data<-transposed_chr_data[,c(3113, 1:3112)]

#Remove first row containing chromosome numbers
#sampleCopyNumbers <- transposed_chr_data[-1,]
#Extract list of sampleIDs from sample matrix
#sampleId <- rownames(sampleCopyNumber)
#Add list of sampleIDs as column in sample matrix
#sampleCopyNumbers <- cbind(sampleCopyNumbers,sampleId)
#Move last column as first column
#sampleCopyNumbers <- sampleCopyNumbers[,c(ncol(sampleCopyNumbers),1:(ncol(sampleCopyNumbers)-1))]

#Convert into a dataframe 
#sampleCopyNumbers <- as.data.frame(sampleCopyNumbers)

##Merge sample matrix and biopsySite matrix
raw_df = merge(transposed_chr_data, biopsySite, by.x = "sampleid", by.y = "sampleId")
raw_df <- raw_df[,c(ncol(raw_df),1:(ncol(raw_df)-1))]

#Arrange as per biopsysite location
df = raw_df %>% arrange(biopsySite)
#Subset Lung and Liver metastatic site wise
newData <- df[c(34:387),]
#Remove first two columns
#newData <- subset(newData, select = -c(biopsySite,sampleid))
#Name columns
colnames_newData <- colnames(newData[-c(1,2)])

#Convert copynumber ranges into categories
#cnAbbr <- cut(newData[,c(1)] , breaks=c(-10,3,5,6,20),labels=c("loss","gain","amp","cnLOH"))

#Create function to convert numerical variables into categorical
my_function <- function(x){
  cnAbbr <- cut(x, breaks=c(-10,2,5,6,20),labels=c("loss","gain","amp","cnLOH"))
  return(cnAbbr)
}

#Apply above function to entire dataframe
newData[,colnames_newData]<- sapply(X= newData[,colnames_newData],
      FUN = my_function)

#Remove Y and X chromosome bins
newData_22 <- subset(newData,select = -c(2899:3114))
#Name columns of newData_23 
#newData_22 <- subset(newData_22, select = -c(1,2))
newData_22_colnames <- colnames(newData_22)


for (x in newData_22_colnames){
  catLength <- length(levels(as.factor(newData_22[[x]])))
  if (catLength < 2){
    newData_22_colnames <- newData_22_colnames[newData_22_colnames != x]

  }
  
}


#Create empty dataframe
fisherBins <- data.frame(binNumber=as.character(),fisherPvalue=as.numeric())

#Insert Pvalues bin number wise in above empty data frame
for (i in newData_22_colnames){
  
  #Create empty data frame and insert biopsySite list as first column
  testdf <- data.frame(biopsySite=newData$biopsySite)
  #Attach second column as bin (looping:attaching each bin one by one)
  testdf$bin <- newData_22[[i]]
  #Repeat for all six combinations of contingency table
  #Subset categories for contingency table
  testdf1 <- filter(testdf, bin == "amp"|bin == "gain")
  
  siteCatLength <- length(levels(as.factor(testdf1$biopsySite)))
  binCatLength <- length(levels(as.factor(testdf1$bin)))
  
  if (siteCatLength > 1 & binCatLength > 1){
    
    TAB <- table(testdf1$biopsySite,testdf1$bin)
    print(TAB)
    fisher_table <- fisher.test(TAB) #Run fisher test on contingency table
    print(fisher_table)
    fisher_pvalue <- fisher_table$p.value #Extract p-value column
    fisherBins[i, ] <- c(i, fisher_pvalue) #Append p-values row by row in data frame
    
  }
  else{
    
  }

}

#Convert p-value column and binNumber into numeric class
fisherBins$fisherPvalue <- as.numeric(fisherBins$fisherPvalue)
fisherBins$binNumber <- as.numeric(fisherBins$binNumber)

#Round off p-values upto three decimal place
fisherBins$fisherPvalue <- lapply(fisherBins$fisherPvalue, round, 3)


#Extract first three columns from chrSep data frame
fisherDf <- subset(chrSep,select = c(chromosome,start,stop))
fisherDf <- tibble::rowid_to_column(fisherDf, "binNumber")

#merge fisherDf and fisherBins to create final p-value table
final_fisher_df <- left_join(fisherDf,fisherBins, by="binNumber")

#Adjust p-values with bonferroni correction
final_fisher_df$fisherPvalue <- vapply(final_fisher_df$fisherPvalue , paste, collapse = ", ", character(1L))
#final_fisher_df$bonferroni <- p.adjust(final_fisher_df$fisherPvalue , method = "bonferroni")


#Remove X and Y chromosome from final_fisher_df
final_fisher_df <- final_fisher_df[-c(2897:3113),]

#Adjust p-values with FDR correction and add new column FDR (corrected P-values)
final_fisher_df$FDR <- p.adjust(final_fisher_df$fisherPvalue , method = "fdr")


#Round off p-values upto three decimal place
final_fisher_df$FDR <- lapply(final_fisher_df$FDR, round, 4)


#Add boneferroni result into final df
#final_fisher_df <- left_join(fisherDf,fisherBins, by="binNumber")

#Export data frame into a csv file
final_fisher_df$FDR <- vapply(final_fisher_df$FDR , paste, collapse = ", ", character(1L))
write.csv(final_fisher_df,file = "C:/Research/Research_project/fishertestdf/ampgain.csv")

#Read p-value CSVs
gainloss <- read.csv("C:/Research/Research_project/fishertestdf/gainloss.csv",sep = ";",header = TRUE)
cnLOHloss <- read.csv("C:/Research/Research_project/fishertestdf/cnLOHloss.csv",sep = ";",header = TRUE)
cnLOHgain <- read.csv("C:/Research/Research_project/fishertestdf/cnLOHgain.csv",sep = ";",header = TRUE)
cnLOHamp <- read.csv("C:/Research/Research_project/fishertestdf/cnLOHamp.csv",sep = ";",header = TRUE)
amploss <- read.csv("C:/Research/Research_project/fishertestdf/amploss.csv",sep = ";",header = TRUE)
ampgain <- read.csv("C:/Research/Research_project/fishertestdf/ampgain.csv",sep = ";",header = TRUE)

#Remove extra index column
#gainloss <- subset(gainloss, select = -c(X))
#cnLOHloss <- subset(cnLOHloss, select = -c(X))
#cnLOHgain <- subset(cnLOHgain, select = -c(X))
#cnLOHamp <- subset(cnLOHamp, select = -c(X))
#amploss <- subset(amploss, select = -c(X))
#ampgain <- subset(ampgain, select = -c(X))



pValue_df <-  left_join(gainloss,cnLOHloss, by = "binNumber")
pValue_df <-  left_join(pValue_df,cnLOHgain, by = "binNumber")
pValue_df <-  left_join(pValue_df,cnLOHamp, by = "binNumber")
pValue_df <-  left_join(pValue_df,amploss, by = "binNumber")
pValue_df <-  left_join(pValue_df,ampgain, by = "binNumber")

write.csv(pValue_df,file = "C:/Research/Research_project/fishertestdf/allPvalues.csv")



