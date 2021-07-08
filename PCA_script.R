#Read copy number data both from TCGA and Hartwig
TCGAliver <- read.csv("C:/Research/Research_project/TCGAbinned/liverOutput.csv",header=TRUE)
TCGAlung <- read.csv("C:/Research/Research_project/TCGAbinned/LungOutput.csv",header=TRUE)
TCGAcolon <- read.csv("C:/Research/Research_project/TCGAbinned/ColonOutput.csv",header=TRUE)
Hartwigliver <- read.csv("C:/Research/Research_project/tcga/HartwigLiver.csv",header=TRUE)
Hartwiglung <- read.csv("C:/Research/Research_project/tcga/HartwigLung.csv",header=TRUE)

#Save bn array
a <- subset(TCGAliver,select=c(1))

#subset sample ID columns
TCGAliver <- subset(TCGAliver,select=-c(1:5))
TCGAlung <- subset(TCGAlung, select=-c(1:5))
TCGAcolon <- subset(TCGAcolon, select=-c(1:5))
Hartwiglung <- subset(Hartwiglung, select=-c(1:3))

#join all the data frame into one
test <- cbind(Hartwigliver,Hartwiglung)
test <- test[-c(2898:3113),]
test1 <- cbind(test,TCGAliver)
test2 <- cbind(test1,TCGAlung)
all_samples <- cbind(test2,TCGAcolon)
all_samples <- subset(all_samples, select=-c(1:3))

#Add bin column

PCAdata <- cbind(all_samples,a)
PCAdata <- PCAdata %>% relocate(bn,.before = "CPCT02010414T")
PCAdata <- lapply(PCAdata, as.numeric)
PCAdata <- as.data.frame(PCAdata)
#Remove zero column to unit variance
#PCAdata <- PCAdata[ , which(apply(PCAdata, 2, var) != 0)]


#transpose dataframe PCAdata
PCAdata$bn <- as.character(PCAdata$bn)
tPCA <- t(PCAdata)
tPCA <- as.data.frame(tPCA)
names(tPCA) <- as.matrix(tPCA[1, ])
tPCA <- tPCA[-1, ]
#tPCA[] <- lapply(tPCA, function(x) type.convert(as.character(x)))
tPCA$sampleid <- row.names(tPCA)
tPCA <- tPCA %>% relocate(sampleid,.before="1")
biopsySite <- c(rep("HartwigLiver",327),rep("HartwigLung",27),rep("TCGALiver",362),rep("TCGALung",950),rep("TCGAColon",290))

#add biopsysite annotation to sampleids
finalPCA <- cbind(tPCA,biopsySite)

#position column
finalPCA <- finalPCA %>% relocate(biopsySite,.after="sampleid")

#Remove sample id and biopsysite column
PCAdf <- subset(finalPCA, select = -c(1,2))

#Convert columns into numeric class
PCAdf[] <- lapply(PCAdf, function(x) type.convert(as.numeric(x)))

#Replace NAs of TCGA data with column mean
for(i in 1:ncol(PCAdf)){
  PCAdf[is.na(PCAdf[,i]), i] <- mean(PCAdf[,i], na.rm = TRUE)
}


#Perform PCA
PCAresult <- prcomp(PCAdf,scale. = TRUE,center = TRUE)

install.packages("factoextra","dendextend")
install.packages("ggplot2")
library("factoextra")
#Plot scree plot
fviz_eig(PCAresult)

#Graph
fviz_pca_ind(PCAresult,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)




pca.var <- PCAresult$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per

##Retrieve percentage variance and cumulative variance
summ = summary(PCAresult)
summ_var = summ$importance[2,] * 100
barplot(summ_var[1:15])

#Create two-dimensional PCA score
Xscores = PCAresult$x
par(mfrow = c(1,1))
my_cols <- c(5,7,8,10,12)
names(my_cols) <- c("Hartwig Liver","Hartwig Lung","TCGA Liver","TCGA Lung","TCGA Colon")
plot(Xscores[,1], Xscores[,2],pch=19,col=my_cols,xlab = "PC1",ylab = "PC2",main = "PCA plot")
legend("topright",legend=c("Hartwig Liver","Hartwig Lung","TCGA Liver","TCGA Lung","TCGA Colon"), fill=c(5,7,8,10,12), border = NA, bty="n",x.intersp = 1,cex=0.7, title.adj = "adj")

#PCAdataframe <- as.data.frame(PC1 = Xscores[,1],PC2 = Xscores[,2],biopsySite = biopsySite)

WGD <- read.csv("C:/Research/Research_project/tcga/WGDsamplesAllData.csv",header = TRUE,sep = ";")

PCAdf["sampleid"] <- rownames(PCAdf)
test <- left_join(PCAdf,WGD,by="sampleid")
#Reoder
test1 <- test[,c(2898,2899,1:2897)]
#Replace na with No string
library(tidyr)
test1$WGD <- test1$WGD %>% replace_na("No") 
WholeGenomeDoubling <- test1$WGD

x <- cbind(PCAdf,Xscores[,1:2])
x$biopsysite <- biopsySite
x$WGD <- WholeGenomeDoubling
PCAdataframe <- x %>% relocate(c(biopsysite,PC1,PC2),.before="1")


#Plot venn diagram
ggplot(x,aes(PC1,PC2,fill=biopsySite,colour=biopsySite))+
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = " "))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = " "))+
  ggtitle("PCA graph : Primary and Metastasis samples:1MB Bins")+
  stat_ellipse(geom = "polygon",col="black",alpha = 0.3)+
  #geom_point(pch = 19, colour = "black",size=0.5)+
  geom_point(aes(shape=WGD))+
  geom_point(pch = 21,size=0.2)+
  
  

  
  


#Visualize PCA results
plot(PCAresult)


#Perform clustering
#Distance
df <- scale(PCAdf)
dtd <- dist(df, method = "euclidean")
head(dtd)

#Linkage
cluster <- hclust(dtd, method = "ward.D")

library(dendextend)
color_branches(cluster)
par(cex=1,font=1)
plot(cluster,labels = NULL, hang = -1,xlab = NULL,sub = NULL)

library(pvclust)
binpv <- pvclust(PCAdf,nboot=100,parallel = FALSE)
plot(binpv)

HCA <- as.dendogram(cluster)
#####This looks ewww#########

##Create new version of PCAdata
##Select 100 bins with highest standard deviation 
## Add column with standard deviation for each bin
PCAdata <- subset(PCAdata,select=-c(1))
x_sd <- transform(PCAdata,SD=apply(PCAdata,1,sd,na.rm= TRUE))

# Sort rows by SD
sort_X <- x_sd[with(x_sd,order(-SD)),]

# Extract the top 100 bins
#bin_X <- head(sort_X,100)

# Remove the SD column
bin_X <- sort_X[1:(length(sort_X)-1)]

# Converting table into matrix
mx_bin <- as.matrix(bin_X)

# Create Dendogram of Genes
# Calculate Euclidean distance

dtd_bins <- dist(mx_bin, method = "euclidean", diag = FALSE , upper = FALSE , p = 2)


#Create bin cluster 
cluster_bins <- hclust(dtd_bins, method = "complete", members = NULL)
library(dendextend)
cluster_bins <- color_branches(cluster_bins)
plot(cluster_bins,hang = -1, main = "Bin Dendogram")

##Create heatmap 
par(mfrow = c(1,1))
rc <- c(rep("green",327),rep("red",27),rep("blue",362),rep("yellow",950), rep("grey",290),ncol(1956))
#mx_bin <- mapply(mx_bin, FUN=as.numeric)
#mx_bin <- matrix(data=mx_bin, ncol=1956, nrow=100)
#Replace missing values with median 
for(i in 1:ncol(mx_bin)){
  mx_bin[is.na(mx_bin[,i]), i] <- mean(mx_bin[,i], na.rm = TRUE)
}


Heatmap(mx_bin,width = unit(8, "cm"), height = unit(8, "cm"))





