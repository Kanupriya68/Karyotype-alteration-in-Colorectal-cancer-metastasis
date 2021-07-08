#Read TCGA segments data (11084 unique patients data)
TCGA_segs <- read.csv("C:/Research/Research_project/TCGA_data/TCGA_mastercalls.abs_segtabs.fixed.txt",
                      sep = "\t")

#Select columns 
WGDdf <- subset(TCGA_segs, select = c(1,2,6,9))

#Calculate weighted mean copy number
WGD <- WGDdf %>% group_by(Sample) %>% #group by patient sample IDs
  summarise(meanploidy = weighted.mean(Modal_Total_CN, Length)) 
#Round weighted mean upto two digits
WGD$meanploidy <- round(WGD$meanploidy, digits = 2)

#Convert mean ploidy copy number values into category
WGD$WGD <- cut(WGD$meanploidy, breaks = c(0,3,Inf), labels = c("No","Yes"), right = FALSE)

#Clean sample ID names to match with biopsy_locations data frame 
WGD$sample <- substr(WGD$Sample, 1, 12)
WGD <- WGD %>% relocate(sample, .after = Sample)
WGD <- subset(WGD, select = -c(1))
colnames(WGD)[1] <- "Sample"

#import df for site annotation
anno <- read.csv("C:/Research/Research_project/tcga/site.csv",header = TRUE,sep = ",")

#subset sample, biopsy location for annotation
anno <- subset(unique(anno, select = c(1,2)))

#match sample IDs and merge both data frames
newdf <- left_join(WGD,anno,by="Sample")

#remove unmatched rows
newdf <- na.omit(newdf)

wgd <- newdf %>% arrange(Biopsy.Location) 
wgd <- wgd[!grepl("No",wgd$WGD),]  

wgdColon <- wgd %>% filter(Biopsy.Location == 'Colon')
#write.csv(wgdColon,"C:/Research/Research_project/tcga/wgdTCGAColon.csv",row.names = FALSE)
wgdLiver <- wgd %>% filter(Biopsy.Location == 'Liver')
#write.csv(wgdLiver,"C:/Research/Research_project/tcga/wgdTCGALiver.csv", row.names = FALSE)
wgdLung <- wgd %>% filter(Biopsy.Location == 'Lung') 
#write.csv(wgdLung,"C:/Research/Research_project/tcga/wgdTCGALung.csv", row.names = FALSE)
  
#Read binned TCGA data (Liver)
Liver_binned <- read.csv("C:/Research/Research_project/TCGAbinned/liverOutput.csv",sep = ",")
Liver_binned <- subset(Liver_binned, select = -c(1,2,5))

#Replace NAs with diploid
Liver_na <- mutate(Liver_binned, across(everything(), ~replace_na(.x, 2)))
Liver_pcf <- mutate_all(Liver_na, function(x) as.numeric(as.character(x)))

#Filter whole genome doubled samples for Liver
wgdLiverSamples <- wgdLiver$Sample
wgdLiverSamples <- gsub(wgdLiverSamples,pattern = "-",replacement = ".")
wgdLiversubset <- Liver_pcf[,wgdLiverSamples]
chromStart <- Liver_pcf[,c(1,2)]
wgdLiver_pcf <- cbind(chromStart,wgdLiversubset)

#Run pcf 
wgdLiverFreq <- pcf(wgdLiver_pcf, normalize = FALSE)

#Plot frequency
plotFreq(wgdLiverFreq,thres.gain = c(2.5,5),thres.loss = 2, main = "WGD Frequency Plot of TCGA : Phenotype-Liver, n=113")

#Read binned TCGA data (Lung)
Lung_binned <- read.csv("C:/Research/Research_project/TCGAbinned/Lungoutput.csv",sep = ",")
Lung_binned <- subset(Lung_binned, select = -c(1,2,5))

#Replace NAs with diploid
Lung_na <- mutate(Lung_binned, across(everything(), ~replace_na(.x, 2)))
Lung_pcf <- mutate_all(Lung_na, function(x) as.numeric(as.character(x)))

#Filter whole genome doubled samples for Liver
wgdLungSamples <- wgdLung$Sample
wgdLungSamples <- gsub(wgdLungSamples,pattern = "-",replacement = ".")
wgdLungsubset <- Lung_pcf[,wgdLungSamples]
chromStart <- Lung_pcf[,c(1,2)]
wgdLung_pcf <- cbind(chromStart,wgdLungsubset)

#Run pcf 
wgdLungFreq <- pcf(wgdLung_pcf, normalize = FALSE)

#Plot frequency
plotFreq(wgdLungFreq,thres.gain = c(2.5,5),thres.loss = 2, main = "WGD Frequency Plot of TCGA : Phenotype-Lung, n=376") 

#Read binned TCGA data (Colon)
Colon_binned <- read.csv("C:/Research/Research_project/TCGAbinned/Colonoutput.csv",sep = ",")
Colon_binned <- subset(Colon_binned, select = -c(1,2,5))

#Replace NAs with diploid
Colon_na <- mutate(Colon_binned, across(everything(), ~replace_na(.x, 2)))
Colon_pcf <- mutate_all(Colon_na, function(x) as.numeric(as.character(x)))

#Filter whole genome doubled samples for Liver
wgdColonSamples <- wgdColon$Sample
wgdColonSamples <- gsub(wgdColonSamples,pattern = "-",replacement = ".")
wgdColonsubset <- Colon_pcf[,wgdColonSamples]
chromStart <- Colon_pcf[,c(1,2)]
wgdColon_pcf <- cbind(chromStart,wgdColonsubset)

#Run pcf 
wgdColonFreq <- pcf(wgdColon_pcf, normalize = FALSE)

#Plot frequency
plotFreq(wgdColonFreq,thres.gain = c(2.5,5),thres.loss = 2, main = "WGD Frequency Plot of TCGA : Phenotype-Colon, n=95") 

#############################Hartwig WGD################################

Hartwig_WGD_Samples_Liver_Lung <- read.csv("C:/Research/Research_project/Hartwig_WGD_Samples_Liver_Lung.csv",header = TRUE,sep = ";")

wgdHartwig_Liver <- Hartwig_WGD_Samples_Liver_Lung$LiverWGD_SampleIDs  

Hartwig_data <- read.csv("C:/Research/Research_project/TCGA_data/dataCN.csv",header = TRUE, sep = ",")  
base <- subset(Hartwig_data, select = -c(1,2,5,6))

base[] <- sapply(base, as.numeric)
base <- mutate(base, across(everything(), ~replace_na(.x, 2)))

wgdLiversubset <- base[,wgdHartwig_Liver]
chromStart <- base[,c(1,2)]
wgdLiver_pcf <- cbind(chromStart,wgdLiversubset)


wgdLiverFreq <- pcf(wgdLiver_pcf, normalize = FALSE)

#Plot frequency
plotFreq(wgdLiverFreq,thres.gain = c(2.5,5),thres.loss = 2, main = "WGD Frequency Plot of Hartwig : Phenotype-Liver, n=226") 

###Hartwig Lung WGD####

wgdHartwig_Lung <- c("CPCT02010474T","CPCT02040102T","CPCT02040188T","CPCT02040226T","CPCT02040264T","CPCT02050184T","CPCT02060259T","CPCT02070232T","CPCT02080084T","CPCT02160038T","CPCT02190007T","CPCT02210064T","CPCT02220009T","CPCT02330036T","CPCT02330043T","CPCT02380015T","CPCT02380016T","DRUP01050025T","DRUP01050031T","DRUP01330008T")

Hartwig_data <- read.csv("C:/Research/Research_project/TCGA_data/dataCN.csv",header = TRUE, sep = ",")  
base <- subset(Hartwig_data, select = -c(1,2,5,6))

base[] <- sapply(base, as.numeric)
base <- mutate(base, across(everything(), ~replace_na(.x, 2)))

wgdLungsubset <- base[,wgdHartwig_Lung]
chromStart <- base[,c(1,2)]
wgdLung_pcf <- cbind(chromStart,wgdLungsubset)


wgdLungFreq <- pcf(wgdLung_pcf, normalize = FALSE)

#Plot frequency
plotFreq(wgdLungFreq,thres.gain = c(2.5,5),thres.loss = 2, main = "WGD Frequency Plot of Hartwig : Phenotype-Lung, n=20") 
