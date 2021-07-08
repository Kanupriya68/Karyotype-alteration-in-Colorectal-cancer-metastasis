
acrocentric <- p_and_q

load("p_and_q_arms_hg19.RData")

#remove acrocentric arms 13,14,15,21,22
load("p_and_q_arms_hg19_no_acrocentric.RData")

allCNdf <- read.csv("C:/Research/Research_project/test/Alldata.csv",header = TRUE)


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

autosomeL <- subset(allCNdf1,select=c(1))

allCNdf2 <- subset(allCNdf1,select=-c(1))




for(i in 1:ncol(allCNdf2)){
  allCNdf2[is.na(allCNdf2[,i]), i] <- mean(allCNdf2[,i], na.rm = TRUE)
}

allCNdf3 <- cbind(autosomeL,allCNdf2)

#Calculate mean of every chromosome p and q arms
PQarmCnDf <- allCNdf3 %>% group_by(autosomes) %>% summarise(across(everything(), mean))



write.csv(PQarmCnDf,file = "C:/Research/Research_project/test/PQArmWgdCorrectAllCN.csv",row.names = FALSE)


hm3 <- t(PQarmCnDf)

names(hm4) <- hm4[1,]
hm4 <- hm4[-1,]

hm5 <- hm4

hm5$sampleid <- sampleId

hm5 <- hm5[40,1:39]




write.csv(PQarmCnDf,file = "C:/Research/Research_project/test/PQarmmeanall.csv",row.names = FALSE)

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


for (i in 1:187){
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

interstingClust <- fraction_df[542:991,3]

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
fit <- kmeans(mydata,4)

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
          k_colors = c("#333399","#ff3300"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE,# Add rectangle around groups
          main = "Karyotypes driving Clustering",
          ylab = "Euclidean Distance", xlab = "Karyotypes")



# K-Means Clustering with 4 clusters
mydata_ <- subset(mydata,select = -c(451:453))
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
ind.coord$Samples <- c(rep("Hartwig Liver",7),
                       rep("TCGA Liver",175),rep("TCGA Lung",187),
                       rep("TCGA Colon",81))
# Data inspection
head(ind.coord)

# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = c("blue","orange"), ellipse = TRUE, ellipse.type = "convex",
  shape = "Samples", size = 2,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )"),
) +
  stat_mean(aes(color = cluster), size = 4)

#Add title
ggpar(clusterplot,main = "2D representation of Clusters with Whole Genome Doubling")















