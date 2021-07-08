


##load TCGA data for site Liver
TCGAliver <-  read.csv("C:/Research/Research_project/TCGAbinned/liverOutput.csv", header = TRUE,sep = ",")

#Subset first few column and just retain sample columns
TCGAliver <- subset(TCGAliver,select = -c(1:5))

##load TCGA data for site Lung
TCGAlung <-  read.csv("C:/Research/Research_project/TCGAbinned/Lungoutput.csv", header = TRUE,sep = ",")

TCGAlung <- subset(TCGAlung,select = -c(1:5))

##load TCGA data for site Colon
TCGAcolon <-  read.csv("C:/Research/Research_project/TCGAbinned/Colonoutput.csv", header = TRUE,sep = ",")

TCGAcolon <- subset(TCGAcolon,select = -c(1:5))

#Import 1Mb Copy number data from Hartwig Liver and Lung
HartwigLiver <- read.csv("C:/Research/Research_project/tcga/HartwigLiver.csv",header = TRUE,sep = ",")
HartwigLung <- read.csv("C:/Research/Research_project/tcga/HartwigLung.csv",header = TRUE,sep = ",")
#Subset Hartwiglung
HartwigLung <- subset(HartwigLung,select = -c(1:3))

#Combine both Liver and lung dataframes
Hartwig <- cbind(HartwigLiver,HartwigLung)

#remove X and Y chromosome data
Hartwig <- Hartwig[-c(2898:3113),]

Hartwig <- subset(Hartwig,select = -c(1:3))

#Create final dataframe (P and q arms for all Primary and Metastasis data)
a <- cbind(Hartwig,TCGAliver)
b <- cbind(a,TCGAlung)
Alldata <- cbind(b,TCGAcolon)

write.csv(Alldata,file = "C:/Research/Research_project/test/Alldata.csv",row.names = FALSE)

library(tidyverse)
hcdf <- tibble::rowid_to_column(Alldata, "bin")

#Convert integer column into numeric
hcdf<- as.data.frame(lapply(hcdf,as.numeric))
df <- subset(hcdf,select=-c(1))

#Hierarchial clustering for non WGD corrected data
library("cluster")
hc_df <- df %>% 
  scale() %>% # Scale the data
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2") # Compute hierarchical clustering

#Visualise clusters
# Visualize using factoextra
# Cut in 5 groups and color by groups

install.packages(c("factoextra", "dendextend"))

par(cex=0.3, mar=c(5, 8, 4, 1))
library("factoextra")
fviz_dend(hc_df, k = 5, # Cut in five groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#a020f0", "#E7B800", "#FC4E07","#bada55"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE,# Add rectangle around groups
main = "Primary and Metastasis Cluster Dendogram (with ploidy)")

#Check if clusters found are truly clusters with pvclust
library(pvclust)

#Parallel P-value based clustering of groups
parpv_result <- parPvclust(cl=NULL, hcdf, method.hclust = "ward.D2",
           method.dist = "euclidean", nboot = 10,
           iseed = NULL)

#Visualise clusters found 

par(cex=1,font=3)
plot(parpv_result, print.pv=TRUE, print.num=TRUE, float=0.01,
     col.pv=c(si=4, au=2, bp=3, edge=8), cex.pv=0.8,edge.width=0.1,
     edge.lty=0.1,labels=FALSE)

#Put rectangles around groups highly supported by data
pvrect(parpv_result, alpha=0.95)

## print the result of multiscale bootstrap resampling
print(parpv_result, digits=3)


#Elucidate significant clusters 
sig_clusters <- pvpick(parpv_result, alpha=0.95,pv="au")

#Print a cluster with high P-value
sig_clusters$clusters[[2]]

#Print its edge number
sig_clusters$edges[2]

#Count sample fraction in each cluster

sampleDF <- data.frame(cluster = numeric(), edge = numeric(), sampleId = character())
names(sampleDF)<- c("cluster","edge","sampleId")


for (i in 1:269){
  print(i)
  print(sig_clusters$edges[[i]])
  print(sig_clusters$clusters[[i]])
  for(j in 1:length(sig_clusters$clusters[[i]])){
    newDframe <- data.frame(c(i),c(sig_clusters$edges[[i]]),c(sig_clusters$clusters[[i]][j]))
    names(newDframe)<- c("cluster","edge","sampleId")
    sampleDF <- rbind(sampleDF,newDframe) 
  }
}

#Extract sampleIDs in an array
sampleId <- colnames(hcdf)

#Create sample type array
SampleSource <- c(rep("Hartwig Liver",327),rep("Hartwig Lung",27),
                  rep("TCGA Liver",362),rep("TCGA Lung",950),
                  rep("TCGA Colon",290))

#Combine Sample IDs and sample source array and create a data frame
ID_source <- cbind(sampleId,SampleSource)

#Merge sampleDF with ID_Source to get data source annotation
fraction_df <- left_join(sampleDF,ID_source,by="sampleId",copy=TRUE)

#Calculate sample fraction in each cluster

siteFraction <- data.frame(cluster=numeric(),HartwigLiver=numeric(),
                           HartwigLung=numeric(),TCGALiver=numeric(),
                           TCGALung=numeric(),TCGAcolon=numeric())

for (i in unique(fraction_df$cluster)){
  df <- fraction_df[fraction_df$cluster == i,]
  
  
  HLi <- df[df$SampleSource == "Hartwig Liver",]
  countHLi <- length(HLi$SampleSource)
  HLu <- df[df$SampleSource == "Hartwig Lung",]
  countHLu <- length(HLu$SampleSource)
  TC <- df[df$SampleSource == "TCGA Colon",]
  countTC <- length(TC$SampleSource)
  TLi <- df[df$SampleSource == "TCGA Liver",]
  countTLi <- length(TLi$SampleSource)
  TLu <- df[df$SampleSource == "TCGA Lung",]
  countTLu <- length(TLu$SampleSource)
  
  totalcount <- length(df$SampleSource)
  
  fractionHLi <- countHLi/totalcount
  fractionHLu <- countHLu/totalcount
  fractionTC <- countTC/totalcount
  fractionTLi <- countTLi/totalcount
  fractionTLu <- countTLu/totalcount
  
  
  tempdf <- data.frame(i,fractionHLi,fractionHLu,fractionTLi,fractionTLu,fractionTC)
  
  names(tempdf) <- c("cluster","HartwigLiver",
                     "HartwigLung","TCGALiver",
                     "TCGALung","TCGAcolon")
  siteFraction <- rbind(siteFraction,tempdf)
}

#Run K-mean clustering on significant clusters of interest
# K-Means Clustering with 5 clusters
#Some cleaning
library(tidyr)
hcdf <- replace_na(hcdf,as.list(colMeans(hcdf,na.rm=T)))
hcdf <- scale(hcdf) #standardize variables
hcdf[is.nan(hcdf)] <- 0

###Partitioning###
# Determine number of clusters
#wss <- (nrow(hcdf)-1)*sum(apply(hcdf,2,var))
#for (i in 2:15) wss[i] <- sum(kmeans(hcdf,
                                     #centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",
     #ylab="Within groups sum of squares")

# Model Based Clustering
#library(mclust)
#model1 <- Mclust(hcdf)
#plot(model1) # plot results
#summary(model1) # display the best model

# K-Means Cluster Analysis
fit <- kmeans(hcdf, 5) 
# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(hcdf, fit$cluster[[34]])


