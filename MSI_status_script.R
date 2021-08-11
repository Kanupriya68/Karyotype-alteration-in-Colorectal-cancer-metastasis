#Read MSS and MSI status sample wise (metastatic)
metadata <- read.csv(file = "C:/Research/Research_project/Colcc_metadata.csv",header = TRUE, sep = ";")

#Extract sampleid and MSS status
status <- subset(metadata, select = c(1,9,42))

#Correct status with missing information
status$msIndelsPerMb[68] <- 'MSS'
status$msIndelsPerMb[81] <- 'MSS'
status$msIndelsPerMb[145] <- 'MSS'
status$msIndelsPerMb[268] <- 'MSS'
status$msIndelsPerMb[484] <- 'MSI'
status$msIndelsPerMb[489] <- 'MSS'
status$msIndelsPerMb[597] <- 'MSS'
status$msIndelsPerMb[623] <- 'MSS'
status$msIndelsPerMb[629] <- 'MSS'

#Filter Liver and lung samples 
library(dplyr)
metastatic_status <- status %>% filter(biopsyDate == "Liver"|biopsyDate == "Lung")
metastatic_status <- subset(metastatic_status,select=c(1,3))
colnames(metastatic_status)[2] <- "MSI.Status"

#Read MSS and MSI status sample wise (primary)
metadata_primary <- read.csv("C:/Research/Research_project/test/TCGA_clinical_MSS_MSI_POLE.csv",header = TRUE,sep = ";")

#Extract sampleid and MSS status
status_pri <- subset(metadata_primary, select = c(1,14))
colnames(status_pri)[1] <- "sampleId"

#Combine all rows from metastatic and primary MSI status dataframes
MSI_status <- rbind(metastatic_status,status_pri)

#Correct names of MSI-H and MSI-L to MSI
MSI_status["MSI.Status"][MSI_status["MSI.Status"] == "MSI-H"] <- "MSI"
MSI_status["MSI.Status"][MSI_status["MSI.Status"] == "MSI-L"] <- "MSI"

#Read cluster numbers for each primary and metastatic sample
cluster_status <- read.csv("C:/Research/Research_project/test/kmean_clusters.csv",header = TRUE, sep = ",")
#Correct cluster order
cluster_status[cluster_status == "cluster2"] <- "clust4"
cluster_status[cluster_status == "cluster1"] <- "clust5"
cluster_status[cluster_status == "cluster3"] <- "clust1"
cluster_status[cluster_status == "cluster4"] <- "clust2"
cluster_status[cluster_status == "cluster5"] <- "clust3"
colnames(cluster_status)[1] <- "sampleId"

#Correct sample ID strings
MSI_status$sampleId <- gsub("-",".",MSI_status$sampleId)


#Join MSI status and cluster status info
cluster_MSI <- left_join(cluster_status,MSI_status,by="sampleId")


#Calculate percentage of MSI and MSS status
test <- cluster_MSI %>% group_by(Cluster) %>% filter(MSI.Status == "MSS" | MSI.Status == "MSI") 
clust1 <- cluster_MSI %>% filter(Cluster == "clust1")
clust1 <- na.omit(clust1)
clust1_rows <- nrow(clust1)
clust1_MSI <- clust1 %>% filter(MSI.Status == "MSI")
clust1_MSS <- clust1 %>% filter(MSI.Status == "MSS")
#Cluster 1 percentages
clust1_MSI_perc <- (nrow(clust1_MSI)/clust1_rows)*100
clust1_MSS_perc <- (nrow(clust1_MSS)/clust1_rows)*100

clust2 <- cluster_MSI %>% filter(Cluster == "clust2")
clust2 <- na.omit(clust2)
clust2_rows <- nrow(clust2)
clust2_MSI <- clust2 %>% filter(MSI.Status == "MSI")
clust2_MSS <- clust2 %>% filter(MSI.Status == "MSS")
#Cluster 2 percentages
clust2_MSI_perc <- (nrow(clust2_MSI)/clust2_rows)*100
clust2_MSS_perc <- (nrow(clust2_MSS)/clust2_rows)*100


clust3 <- cluster_MSI %>% filter(Cluster == "clust3")
clust3 <- na.omit(clust3)
clust3_rows <- nrow(clust3)
clust3_MSI <- clust3 %>% filter(MSI.Status == "MSI")
clust3_MSS <- clust3 %>% filter(MSI.Status == "MSS")
#Cluster 3 percentages
clust3_MSI_perc <- (nrow(clust3_MSI)/clust3_rows)*100
clust3_MSS_perc <- (nrow(clust3_MSS)/clust3_rows)*100

clust4 <- cluster_MSI %>% filter(Cluster == "clust4")
clust4 <- na.omit(clust4)
clust4_rows <- nrow(clust4)
clust4_MSI <- clust4 %>% filter(MSI.Status == "MSI")
clust4_MSS <- clust4 %>% filter(MSI.Status == "MSS")
#Cluster 4 percentages
clust4_MSI_perc <- (nrow(clust4_MSI)/clust4_rows)*100
clust4_MSS_perc <- (nrow(clust4_MSS)/clust4_rows)*100

clust5 <- cluster_MSI %>% filter(Cluster == "clust5")
clust5 <- na.omit(clust5)
clust5_rows <- nrow(clust5)
clust5_MSI <- clust5 %>% filter(MSI.Status == "MSI")
clust5_MSS <- clust5 %>% filter(MSI.Status == "MSS")
#Cluster 4 percentages
clust5_MSI_perc <- (nrow(clust5_MSI)/clust5_rows)*100
clust5_MSS_perc <- (nrow(clust5_MSS)/clust5_rows)*100
