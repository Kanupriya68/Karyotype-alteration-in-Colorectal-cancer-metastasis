TCGA_cndata <- read.csv(file = "C:/Research/Research_project/TCGA_data/TCGA_mastercalls.abs_segtabs.fixed.txt",
                        sep = "\t",header = TRUE)
TCGA_data <- subset(TCGA_cndata,select = c(1:4))
TCGA_data$copynumber <- TCGA_cndata$Modal_Total_CN



#Create TCGA data bins and append to make a final dataframe 
zero_fifty <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-0-50.csv",
         sep = ",",header = TRUE, na.strings=c("","NA"))

#Read all TCGA binned files one by one
fiftyone_100 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-51-100.csv",
                       sep = ",",header = TRUE,na.strings=c("","NA"))
fiftyone_100 <- subset(fiftyone_100,select = -c(1:6))

hundredone_150 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-101-150.csv",
                         sep = ",",header = TRUE,na.strings=c("","NA"))
hundredone_150 <- subset(hundredone_150,select = -c(1:6))

one51_200 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-151-200.csv",
                           sep = ",",header = TRUE,na.strings=c("","NA"))
one51_200 <- subset(one51_200,select = -c(1:6))

two01_250 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-201-250.csv",
                      sep = ",",header = TRUE,na.strings=c("","NA"))
two01_250 <- subset(two01_250,select = -c(1:6))

two51_300 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-251-300.csv",
                      sep = ",",header = TRUE,na.strings=c("","NA"))
two51_300 <- subset(two51_300,select = -c(1:6))

three01_350 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-301-350.csv",
                      sep = ",",header = TRUE,na.strings=c("","NA"))
three01_350 <- subset(three01_350,select = -c(1:6))

three51_400 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-351-400.csv",
                        sep = ",",header = TRUE,na.strings=c("","NA"))
three51_400 <- subset(three51_400,select = -c(1:6))

four01_450 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-401-450.csv",
                        sep = ",",header = TRUE,na.strings=c("","NA"))
four01_450 <- subset(four01_450,select = -c(1:6))

four51_500 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-451-500.csv",
                        sep = ",",header = TRUE,na.strings=c("","NA"))
four51_500 <- subset(four51_500,select = -c(1:6))

five01_550 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-501-550.csv",
                        sep = ",",header = TRUE,na.strings=c("","NA"))
five01_550 <- subset(five01_550,select = -c(1:6))

five51_600 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-551-600.csv",
                        sep = ",",header = TRUE,na.strings=c("","NA"))
five51_600 <- subset(five51_600,select = -c(1:6))

six01_650 <- read.csv(file = "C:/Research/Research_project/TCGAbinned/Tcga-601-650.csv",
                        sep = ",",header = TRUE,na.strings=c("","NA"))
six01_650 <- subset(six01_650,select = -c(1:6))

#Merge all binned TCGA files to create a dataframe
base <- cbind(zero_fifty,fiftyone_100)
base1 <- cbind(base,hundredone_150)
base2 <- cbind(base1,one51_200)
base3 <- cbind(base2,two01_250)
base4 <- cbind(base3,two51_300)
base5 <- cbind(base4,three01_350)
base6 <- cbind(base5,three51_400)
base7 <- cbind(base6,four01_450)
base8 <- cbind(base7,four51_500)
base9 <- cbind(base8,five01_550)
base10 <- cbind(base9,five51_600)
base11 <- cbind(base10,six01_650)

#Remove X and Y chromosome bins
TCGA <- base11[-c(2898:3113),]
TCGA_pcf <- subset(TCGA, select = -c(1,2,5,6))
#TCGA_pcf <- imputeMissing(TCGA_pcf,method = "pcf")
#Replace NAs with diploid value : "2"
TCGA_na <- mutate(TCGA_pcf, across(everything(), ~replace_na(.x, 2)))
TCGA_pcf <- mutate_all(TCGA_na, function(x) as.numeric(as.character(x)))


library(spatstat)

TCGA_freq <- pcf(TCGA_pcf)




#Import Hartwig copynumber data 
Hartwig_data <- read.csv("C:/Research/Research_project/TCGA_data/dataCN.csv",header = TRUE, sep = ",")

#Create dataframe for pcf on Hartwig copynumber data
#Remove bn.binnumber and arm column
Hartwig_pcf <- subset(Hartwig_data, select = -c(1,2,5,6))
#Remove X Y chromosomes
Hartwig_pcf <- Hartwig_pcf[-c(2898:3113),]


library(dplyr)
Hartwig_pcf <- mutate_all(Hartwig_pcf, function(x) as.numeric(as.character(x)))

library(spatstat)
 Hartwig_f <- pcf(Hartwig_pcf)

#Plot frequency graph of Colcc metastatic copynumber data (Hartwig)
plotFreq(segments=Hartwig_f,thres.gain=2.5,thres.loss=1.5,main = "Copynumber Loss & Gain : Hartwig")

#Plot frequency graph of primary tumour copynumber data (TCGA)
plotFreq(segments=TCGA_freq,thres.gain=2.5,thres.loss=1.5,main = "Copynumber Loss & Gain : TCGA")

write.csv(Hartwig_freq, file = "C:/Research/Research_project/test/Hartwig_freq.csv",row.names = FALSE)
