
acrocentric <- p_and_q

load("p_and_q_arms_hg19.RData")

#remove acrocentric arms 13,14,15,21,22
load("p_and_q_arms_hg19_no_acrocentric.RData")

allCNdf <- read.csv("C:/Research/Research_project/test/Alldata.csv",header = TRUE)

chr_start_stop <- subset(allCNdf,select = c(1,2,3))
allCNdf <- subset(allCNdf,select = -c(1,2,3))




meanNewData <- data.frame(sampleID=character(),meanploidy=numeric())
for (i in colnames(allCNdf)){
  average <- mean(allCNdf[[i]],na.rm=TRUE)
  df_newmean <- data.frame(i,average)
  names(df_newmean) <- c("sampleID","meanploidy")
  meanNewData <- rbind(meanNewData,df_newmean)
}


#Ploidy correction (remove whole genome doubling)

for (i in meanNewData$sampleID){
  as <- meanNewData[meanNewData$sampleID == i,]
  mean_ploidy <- as[[2]]
  if (mean_ploidy > 3 & mean_ploidy < 3.5){
    allCNdf[i] <- allCNdf[i]-1
  }
  else if (mean_ploidy >= 3.5 & mean_ploidy < 4.5){
    allCNdf[i] <- allCNdf[i]-2
  }
  else if (mean_ploidy >= 4.5 & mean_ploidy < 5.5){
    allCNdf[i] <- allCNdf[i]-3
  }
  
  else if (mean_ploidy >= 5.5 & mean_ploidy < 6.5){
    allCNdf[i] <- allCNdf[i]-4
  }
  
  else if (mean_ploidy >= 6.5 & mean_ploidy < 7.5){
    allCNdf[i] <- allCNdf[i]-5
  }
  
  else if (mean_ploidy >= 7.5 & mean_ploidy < 8.5){
    allCNdf[i] <- allCNdf[i]-6
  }
  
}


allCNdf <- cbind(chr_start_stop,allCNdf)

#write.csv(allCNdf,file = "C:/Research/Research_project/test/allCNWGDCorrected.csv")



allCNdf$arm = "p"
library(dplyr)
allCNdf <- allCNdf %>% relocate(arm, .before = CPCT02010414T)

#Calculate pand q arms for df
for (i in 1:nrow(acrocentric)) {
  allCNdf$arm[(allCNdf$chromosome == acrocentric$Chromosome[i] & allCNdf$start >= acrocentric$Start[i] & allCNdf$stop <= acrocentric$End[i])] <- acrocentric$cyto[i]
  allCNdf$arm[(allCNdf$chromosome == acrocentric$Chromosome[i] & allCNdf$start >= acrocentric$Start[i] & allCNdf$stop >= acrocentric$End[i])] <- "q"
  
} 

#Subset start,stop column 
allCNdf <- subset(allCNdf,select = -c(2,3))

#Remove p arms of acrocentric chromosomes
allCNdf1 <- allCNdf[-c(2094:2111,2210:2227,2318:2336,2797:2810,2846:2860),]

#Paste chromosome and arm string together in a new column
allCNdf1$autosomes <- apply(allCNdf1[,c("chromosome","arm")],1,paste,collapse="")
allCNdf1 <- allCNdf1 %>% relocate(autosomes,.after=arm)




allCNdf1 <- subset(allCNdf1,select=-c(1,2))

#write.csv(allCNdf1,file = "C:/Research/Research_project/test/allCNdf1.csv", row.names = FALSE)
dfmain <- subset(allCNdf1,select = -c(1))

for(i in 1:ncol(dfmain)){
  dfmain[is.na(dfmain[,i]), i] <- mean(dfmain[,i], na.rm = TRUE)
}

dfmain$autosomes <- allCNdf1$autosomes
allCNdf1 <- dfmain %>% relocate(autosomes, .before = CPCT02010414T)

#Calculate mean of every chromosome p and q arms
PQarmdf <- allCNdf1 %>% group_by(autosomes) %>% summarise(across(everything(), mean))

auto <- PQarmCnDf$autosomes
PQarmCnDf <- subset(PQarmCnDf,select=-c(1))

write.csv(PQarmdf,file = "C:/Research/Research_project/tcga/PQarmCNdfFinal.csv",row.names = FALSE)

#Convert integer column into numeric
PQarmCnDf<- as.data.frame(lapply(PQarmCnDf,as.numeric))
row.names(PQarmCnDf) <- auto

#transpose the data frame
tPQarms <- as.data.frame(t(PQarmCnDf))


#Visualise clusters
install.packages(c("factoextra", "dendextend"))

par(cex=0.3, mar=c(5, 8, 4, 1))
library("factoextra")
fviz_dend(PQclust, k = 5, # Cut in three groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#a020f0", "#FC4E07","#FE4C78","#DAFF33"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE,# Add rectangle around groups
          main = "Primary and Metastatic Cluster : Samples = 1956",
          ylab = "Euclidean Distance", xlab = "Patient samples")

dev.off()

#Cluster Validation
install.packages("pvclust")
library("pvclust")

#data cleaning before pvclust
#Replace Nans with mean column value
dfgh <- data.frame(arm = character())

#Replace na with mean of corresponding columns
for(i in 1:ncol(PQarmCnDf)){
  PQarmCnDf[is.na(PQarmCnDf[,i]), i] <- mean(PQarmCnDf[,i], na.rm = TRUE)
}

NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
PQarmCnDf <- replace(PQarmCnDf, TRUE, lapply(PQarmCnDf, NA2mean))

PQarmCnDf <- as.data.frame(PQarmCnDf)

#Replace Nan with zero
PQarmCnDf <- rapply( PQarmCnDf, f=function(x) ifelse(is.nan(x),0,x), how="replace" )

#Cluster samples based on p-value
#Calculate p-value with pvclust and perform clustering
pvpq <- pvclust(PQarmCnDf, method.hclust="ward.D2",
                method.dist="euclidean",
                nboot=100, parallel=FALSE,
                store=FALSE, weight=FALSE, iseed=NULL, quiet=FALSE)

#Visualise clusters found 

par(cex=1,font=3)
plot(pvpq, print.pv=TRUE, print.num=FALSE, float=0.01,
     col.pv=c(si=4, au=2, bp=3, edge=8), cex.pv=0.7,edge.width=0.1,
     edge.lty=0.1,labels=FALSE,cex = 0.8)


#Put rectangles around groups highly supported by data
pvrect(pvpq, alpha=0.95)

## print the result of multiscale bootstrap resampling
print(pvpq, digits=3)

#Elucidate significant clusters 
sig_clusters <- pvpick(pvpq, alpha=0.95,pv="au")

#Print a cluster with high P-value
sig_clusters$clusters[[2]]

#Print its edge number
sig_clusters$edges[2]

#Count sample fraction in each cluster

sampleDF <- data.frame(cluster = numeric(), edge = numeric(), sampleId = character())
names(sampleDF)<- c("cluster","edge","sampleId")


for (i in 1:189){
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
sampleId <- colnames(PQarmCnDf)

#Create sample type array
SampleSource <- c(rep("Hartwig Liver",327),rep("Hartwig Lung",27),
                  rep("TCGA Liver",362),rep("TCGA Lung",950),
                  rep("TCGA Colon",290))

#Combine Sample IDs and sample source array and create a data frame
ID_source <- cbind(sampleId,SampleSource)

#Merge sampleDF with ID_Source to get data source annotation
fraction_df <- left_join(sampleDF,ID_source,by="sampleId",copy=TRUE)

interstingClust <- fraction_df[c(663:666,1211:1217),3]

#Extract samples that are clustered
kmeandf <- subset(PQarmCnDf,select = c(interstingClust))
mydata <- kmeandf

#mydata <- na.omit(mydata)
#mydata <- scale(mydata)
#mydata <- as.data.frame(mydata)
#mydata <- rapply(mydata, f=function(x) ifelse(is.nan(x),0,x), how="replace" )

# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#K-means cluster analysis
set.seed(123)
fit <- kmeans(mydata,2)

# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)

# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

smydata <- scale(mydata)

# Ward Hierarchical Clustering
d <- dist(smydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D")

fviz_dend(fit, k = 2, # Cut in two groups
          cex = 0.5, # label size
          k_colors = c("#0033cc","#ff0066"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE,# Add rectangle around groups
          main = "Whole Genome corrected Karyotypes driving Clustering",
          ylab = "Euclidean Distance", xlab = "Karyotypes")



# K-Means Clustering with 2 clusters
mydata_ <- subset(mydata,select = -c(12))
tmydata <- t(mydata_)
set.seed(124)
tmydata <- scale(tmydata)
fit1 <- kmeans(tmydata, 2)

# append cluster assignment
newdf1 <- data.frame(tmydata, fit1$cluster)
newdf1 <- newdf1 %>% relocate(fit1.cluster,.before="X.1p")


# get cluster means
aggregate(newdf1,by=list(fit1$cluster),FUN=mean)
fit1$size
fit1$cluster


# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster)

clusplot(newdf1, fit1$cluster, color=TRUE, shade=TRUE,
         labels=2, lines=0,main = "2D representation of Clusters")
dev.off()
# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(newdf1, fit1$cluster)

#Better cluster plot
library(ggpubr)
fviz_cluster(fit1, newdf1,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

# Dimension reduction using PCA
res.pca <- prcomp(newdf1,  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(fit1$cluster)
# Add Species groups from the original data sett
ind.coord$Samples <- c("Hartwig Liver","TCGA Liver","TCGA Lung","TCGA Colon",
                       "Hartwig Liver","TCGA Liver",rep("TCGA Lung",5))
# Data inspection
head(ind.coord)

# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

cluster1plot <- ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = c("blue","red"), ellipse = TRUE, ellipse.type = "convex",
  shape = "Samples", size = 2,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )"),
) +
  stat_mean(aes(color = cluster), size = 4)

#Add title
ggpar(cluster1plot,main = "2D representation of Clusters with Whole Genome Doubling (n=459)")
















