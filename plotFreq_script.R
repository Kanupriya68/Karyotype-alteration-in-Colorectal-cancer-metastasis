#Read TCGA segments data (11084 unique patients data)
TCGA_segs <- read.csv("C:/Research/Research_project/TCGA_data/TCGA_mastercalls.abs_segtabs.fixed.txt",
                      sep = "\t")

#Read TCGA metadata with biopsy locations 
clinical <- read.csv("C:/Research/Research_project/TCGA_data/clinical.tsv",
                      sep = "\t")

#Check number of unique biopsy locations
#biopsy_locations <- unique(clinical$site_of_resection_or_biopsy)


#Arrange data frame as per biopsy location and filter out all Liver,Lung and colon samples
library(dplyr)
a <- clinical %>% arrange(site_of_resection_or_biopsy)%>%
  filter(site_of_resection_or_biopsy == "Lung, NOS" | site_of_resection_or_biopsy == "Lower lobe, lung" |
           site_of_resection_or_biopsy == "Upper lobe, lung" | site_of_resection_or_biopsy == "Middle lobe, lung"|
           site_of_resection_or_biopsy == "Liver"|site_of_resection_or_biopsy == "Colon, NOS"|
           site_of_resection_or_biopsy == "Sigmoid colon"|site_of_resection_or_biopsy == "Transverse colon"|
           site_of_resection_or_biopsy == "Descending colon"|site_of_resection_or_biopsy == "Ascending colon")


#Rename samples uniformly as per biopsy location site
b <- subset(a,select = c(2,116))    
b[b == "Lung, NOS"] <- "Lung"
b[b == "Upper lobe, lung"] <- "Lung"
b[b == "Lower lobe, lung"] <- "Lung"
b[b == "Middle lobe, lung"] <- "Lung"
b[b == "Colon, NOS"] <- "Colon"
b[b == "Sigmoid colon"] <- "Colon"
b[b == "Transverse colon"] <- "Colon"
b[b == "Descending colon"] <- "Colon"
b[b == "Ascending colon"] <- "Colon"

#Rearrange again
biopsy_locations <- b %>% arrange(site_of_resection_or_biopsy)
#Rename sample ID column
colnames(biopsy_locations)[1] <- "Sample"
colnames(biopsy_locations)[2] <- "Biopsy Location"

#Clean sample ID names to match with biopsy_locations data frame 
TCGA_segs$sample <- substr(TCGA_segs$Sample, 1, 12)
TCGA_segs <- TCGA_segs %>% relocate(sample, .after = Sample)
TCGA_segs <- subset(TCGA_segs, select = -c(1))
colnames(TCGA_segs)[1] <- "Sample"
#colnames(TCGA_segs)[9] <- "CopyNumber"

#Merge both dataframe to match as per biopsy location
df <- left_join(TCGA_segs,biopsy_locations, by = "Sample")
#relocate Biopsy location column next to Sample column
df <- df %>% relocate(`Biopsy Location`, .after = Sample)

#Remove unmatched sample id by removing rows with NA in biopst location column
df <- na.omit(df)

#Remove duplicated rows
df <- df[!duplicated(df),]

#write.csv(df, file = "C:/Research/Research_project/tcga/site.csv", row.names = FALSE)

#Select necessary column to run pcf
pcf_df <- subset(df, select = c(1:10))
pcf_df <- subset(pcf_df, select = -c(8,9))

#Create df for Liver
Liver_df <- pcf_df %>% filter(`Biopsy Location` == "Liver")



#Create df for Lung
Lung_df <- pcf_df %>% filter(`Biopsy Location` == "Lung")



#Create df for Colon
Colon_df <- pcf_df %>% filter(`Biopsy Location` == "Colon")



#write.csv(Colon_df,file = "C:/Research/Research_project/tcga/colon.csv",row.names = FALSE)

#write.csv(Liver_df,file = "C:/Research/Research_project/tcga/Liver.csv", row.names = FALSE)

#write.csv(Lung_df,file = "C:/Research/Research_project/tcga/Lung.csv", row.names = FALSE)

#Import Hartwig copynumber data 
Hartwig_data <- read.csv("C:/Research/Research_project/TCGA_data/dataCN.csv",header = TRUE, sep = ",")
HartwigPhenotype <- read.csv("C:/Research/Research_project/binnedCnDataframe/liver_Lung_patientIDs.csv",header = TRUE, sep = ";")

ID <- subset(Hartwig_data, select = -c(1:6))
tID <- t(ID)
Hdf <- cbind(sampleid = rownames(tID), tID)
rownames(Hdf) <- 1:nrow(Hdf)
colnames(Hdf)[1] <- "sampleId" 
Hdf <- as.data.frame(Hdf)
Hartwig_phen <- left_join(Hdf,HartwigPhenotype, by= "sampleId")
Hartwig_phen <- Hartwig_phen %>% relocate(biopsySite, .after= "sampleId")%>%
  arrange(biopsySite)

##########################
######Hartwig_Liver#######
##########################

Hartwig_Liver <- Hartwig_phen[1:327,]
Hartwig_Liver <- subset(Hartwig_Liver, select= -c(2))
Hartwig_Liver <- t(Hartwig_Liver)
colnames(Hartwig_Liver) <- Hartwig_Liver[1,]
Hartwig_Liver <- Hartwig_Liver[-1,]

append <- subset(Hartwig_data, select= c(3,4,5))
Hartwig_Liver <- cbind(append,Hartwig_Liver)

#write.csv(Hartwig_Liver,file = "C:/Research/Research_project/tcga/HartwigLiver.csv", row.names = FALSE)

#write.csv(Hartwig_Lung,file = "C:/Research/Research_project/tcga/HartwigLung.csv", row.names = FALSE)

Hartwig_Liver <- subset(Hartwig_Liver, select= -c(3))
Hartwig_na <- mutate(Hartwig_Liver, across(everything(), ~replace_na(.x, 2)))
Hartwig_Liver <- mutate_all(Hartwig_na, function(x) as.numeric(as.character(x)))

#Run pcf on Hartwig Liver df
Hartwig_Liver_pcf <- pcf(Hartwig_Liver, normalize = FALSE)


#Plot frequency of Hartwig liver df

plotFreq(Hartwig_Liver_pcf,thres.gain = c(2,5), thres.loss = 2, main="Frequency Plot of Hartwig : Phenotype- Liver , n=327")

#WGD correction
wgd_HartwigLiver <- read.csv("C:/Research/Research_project/Hartwig_WGD_Samples_Liver_Lung.csv",sep = ";",header = TRUE)

wgdLiversubset <- Hartwig_Liver[ ,wgd_HartwigLiver$LiverWGD_SampleIDs]

wgdcorrectedLiver <- Hartwig_Liver

meanHartwig_Liver <- data.frame(sampleID=character(),meanploidy=numeric())

#Calculate mean ploidy for Hartwig Liver
for (i in colnames(Hartwig_Liver[3:329])){
  average <- mean(Hartwig_Liver[[i]])
  df_newmean <- data.frame(i,average)
  names(df_newmean) <- c("sampleID","meanploidy")
  meanHartwig_Liver <- rbind(meanHartwig_Liver,df_newmean)
}

#Reduce whole genome doubled copy number values
for (i in wgd_HartwigLiver$LiverWGD_SampleIDs){
  as <- meanHartwig_Liver[meanHartwig_Liver$sampleID == i,]
  mean_ploidy <- as[[2]]
  if (mean_ploidy > 3 & mean_ploidy < 3.5){
    wgdcorrectedLiver[i] <- wgdcorrectedLiver[i]-1
  }
  else if (mean_ploidy >= 3.5){
    wgdcorrectedLiver[i] <- wgdcorrectedLiver[i]-2
  }
  
}

#Rerun pcf on wgd corrected df
wgdLiver_pcf <- pcf(wgdcorrectedLiver,normalize = FALSE)

plotFreq(wgdLiver_pcf,thres.gain = 2,thres.loss = 2,main="Frequency Plot of Hartwig : Liver, n=327")



#######################
####Hartwig_Lung#######
#######################

#Subset Lung sample from Hartwig main df
Hartwig_Lung <- Hartwig_phen[328:354,]
#Remove extra columns
Hartwig_Lung <- subset(Hartwig_Lung, select= -c(2))
#Transpose df to get sample IDs as column names
Hartwig_Lung <- t(Hartwig_Lung)
#Name columns as per sample ID
colnames(Hartwig_Lung) <- Hartwig_Lung[1,]
#Remove index row and column
Hartwig_Lung <- Hartwig_Lung[-1,]
#Add chromosome, start , end columns
Hartwig_Lung <- cbind(append,Hartwig_Lung)
#Remove extra column
Hartwig_Lung <- subset(Hartwig_Lung, select= -c(3))
#Replace all NAs with diploid value '2'
Hartwig_Lung_na <- mutate(Hartwig_Lung, across(everything(), ~replace_na(.x, 2)))
#Convert character into numeric class
Hartwig_Lung <- mutate_all(Hartwig_Lung_na, function(x) as.numeric(as.character(x)))
#Run pcf to get segments
Hartwig_Lung_pcf <- pcf(Hartwig_Lung, normalize = FALSE)
#Plot segments in frequency plot to visualize gain and loss percentage
plotFreq(Hartwig_Lung_pcf,thres.gain = c(2,5), thres.loss = 2, main="Frequency Plot of Hartwig : Phenotype- Lung , n=27")

#WGD correction
#wgd_HartwigLung <- read.csv("C:/Research/Research_project/Hartwig_WGD_Samples_Liver_Lung.csv",sep = ";",header = TRUE)
wgdHartwig_Lung <- c("CPCT02010474T","CPCT02040102T","CPCT02040188T","CPCT02040226T","CPCT02040264T","CPCT02050184T","CPCT02060259T","CPCT02070232T","CPCT02080084T","CPCT02160038T","CPCT02190007T","CPCT02210064T","CPCT02220009T","CPCT02330036T","CPCT02330043T","CPCT02380015T","CPCT02380016T","DRUP01050025T","DRUP01050031T","DRUP01330008T")
wgdLungsubset <- Hartwig_Lung[ ,wgdHartwig_Lung]

wgdcorrecteddf <- Hartwig_Lung

meanHartwig_Lung <- data.frame(sampleID=character(),meanploidy=numeric())

for (i in colnames(Hartwig_Lung[3:29])){
  average <- mean(Hartwig_Lung[[i]])
  df_newmean <- data.frame(i,average)
  names(df_newmean) <- c("sampleID","meanploidy")
  meanHartwig_Lung <- rbind(meanHartwig_Lung,df_newmean)
}

wgdcorrectedlung <- Hartwig_Lung

for (i in wgdHartwig_Lung){
  as <- meanHartwig_Lung[meanHartwig_Lung$sampleID == i,]
  mean_ploidy <- as[[2]]
  if (mean_ploidy > 3 & mean_ploidy < 3.5){
    wgdcorrectedlung[i] <- wgdcorrectedlung[i]-1
  }
  else if (mean_ploidy >= 3.5){
    wgdcorrectedlung[i] <- wgdcorrectedlung[i]-2
  }
  
}

#Rerun pcf on wgd corrected df
wgdLung_pcf <- pcf(wgdcorrectedlung,normalize = FALSE)

plotFreq(wgdLung_pcf,thres.gain = 2,thres.loss = 2,main="Frequency Plot of Hartwig : Lung, n=27")







############################
########TCGA_Colon##########
############################


#Read binned TCGA data (Liver, Lung, Colon)
Colon_binned <- read.csv("C:/Research/Research_project/TCGAbinned/Colonoutput.csv",sep = ",")
Colon_binned <- subset(Colon_binned, select = -c(1,2,5))

#Replace NAs with diploid
Colon_na <- mutate(Colon_binned, across(everything(), ~replace_na(.x, 2)))
Colon_pcf <- mutate_all(Colon_na, function(x) as.numeric(as.character(x)))

#Run pcf 
ColonFreq <- pcf(Colon_pcf, normalize = FALSE)

#Plot frequency 
plotFreq(ColonFreq,thres.gain = c(2,5),thres.loss = 2, main = "Frequency Plot of TCGA : Phenotype-Colon, n=290")
 
#WGD correction
wgd_Colon <- read.csv("C:/Research/Research_project/tcga/wgdTCGAColon.csv",sep = ",",header = TRUE)
wgd_Colon$Sample <- gsub(wgd_Colon$Sample,pattern = "-",replacement = ".")
wgdColonsubset <- Colon_pcf[ ,wgd_Colon$Sample]
#chromStart <- Colon_pcf[,c(1,2)]
#wgdColon_pcf <- cbind(chromStart,wgdColonsubset)
wgdcorrecteddf <- Colon_pcf
for (i in wgd_Colon$Sample){
  as <- wgd_Colon[wgd_Colon$Sample == i,]
  mean_ploidy <- as[[2]]
  if (mean_ploidy > 3 & mean_ploidy < 3.5){
    wgdcorrecteddf[i] <- wgdcorrecteddf[i]-1
  }
  else if (mean_ploidy >= 3.5){
    wgdcorrecteddf[i] <- wgdcorrecteddf[i]-2
  }
  
}

#Rerun pcf on wgd corrected df
wgdColon_pcf <- pcf(wgdcorrecteddf,normalize = FALSE)

plotFreq(wgdColon_pcf,thres.gain = 2,thres.loss = 2,main="Frequency Plot of TCGA : Colon, n=290")


#################################
#########TCGA_Liver##############
#################################

#Read binned TCGA data (Liver)
Liver_binned <- read.csv("C:/Research/Research_project/TCGAbinned/liverOutput.csv",sep = ",")
Liver_binned <- subset(Liver_binned, select = -c(1,2,5))

#Replace NAs with diploid
Liver_na <- mutate(Liver_binned, across(everything(), ~replace_na(.x, 2)))
Liver_pcf <- mutate_all(Liver_na, function(x) as.numeric(as.character(x)))


#Run pcf 
LiverFreq <- pcf(Liver_pcf, normalize = FALSE)

#Plot frequency

plotFreq(LiverFreq,thres.gain = c(2,5),thres.loss = 2, main = "Frequency Plot of TCGA : Phenotype-Liver, n=362")

#WGD correction
wgd_Liver <- read.csv("C:/Research/Research_project/tcga/wgdTCGALiver.csv",sep = ",",header = TRUE)
wgd_Liver$Sample <- gsub(wgd_Liver$Sample,pattern = "-",replacement = ".")
wgdLiversubset <- Liver_pcf[ ,wgd_Liver$Sample]
wgdcorrecteddf <- Liver_pcf
for (i in wgd_Liver$Sample){
  as <- wgd_Liver[wgd_Liver$Sample == i,]
  mean_ploidy <- as[[2]]
  if (mean_ploidy > 3 & mean_ploidy < 3.5){
    wgdcorrecteddf[i] <- wgdcorrecteddf[i]-1
  }
  else if (mean_ploidy >= 3.5){
    wgdcorrecteddf[i] <- wgdcorrecteddf[i]-2
  }
  
}

#Rerun pcf on wgd corrected df
wgdLiver_pcf <- pcf(wgdcorrecteddf,normalize = FALSE)

plotFreq(wgdLiver_pcf,thres.gain = 2,thres.loss = 2,main="Frequency Plot of TCGA : Liver, n=362")

####################################
#############TCGA_Lung##############
####################################


#Read binned TCGA data (Lung)
Lung_binned <- read.csv("C:/Research/Research_project/TCGAbinned/Lungoutput.csv",sep = ",")
Lung_binned <- subset(Lung_binned, select = -c(1,2,5))

#Replace NAs with diploid
library(dplyr)
library(tidyr)
Lung_na <- mutate(Lung_binned, across(everything(), ~replace_na(.x, 2)))
Lung_pcf <- mutate_all(Lung_na, function(x) as.numeric(as.character(x)))
#Lung_pcf$chromosome <- as.character(Lung_pcf$chromosome)

#Run pcf 
library(spatstat)
LungFreq <- pcf(Lung_pcf, normalize = FALSE)

#Plot frequency 
plotFreq(LungFreq, thres.gain = c(2,5),thres.loss = 2, main = "Frequency Plot of TCGA : Phenotype-Lung, n=950")

#WGD correction
wgd_Lung <- read.csv("C:/Research/Research_project/tcga/wgdTCGALung.csv",sep = ",",header = TRUE)
wgd_Lung$Sample <- gsub(wgd_Lung$Sample,pattern = "-",replacement = ".")
wgdLungsubset <- Lung_pcf[ ,wgd_Lung$Sample]

wgdcorrecteddf <- Lung_pcf

for (i in wgd_Lung$Sample){
  as <- wgd_Lung[wgd_Lung$Sample == i,]
  mean_ploidy <- as[[2]]
  if (mean_ploidy > 3 & mean_ploidy < 3.5){
    wgdcorrecteddf[i] <- wgdcorrecteddf[i]-1
  }
  else if (mean_ploidy >= 3.5){
    wgdcorrecteddf[i] <- wgdcorrecteddf[i]-2
  }
  
}

#Rerun pcf on wgd corrected df
wgdLung_pcf <- pcf(wgdcorrecteddf,normalize = FALSE)

plotFreq(wgdLung_pcf,thres.gain = 2,thres.loss = 2,main="Frequency Plot of TCGA : Lung, n=950")
