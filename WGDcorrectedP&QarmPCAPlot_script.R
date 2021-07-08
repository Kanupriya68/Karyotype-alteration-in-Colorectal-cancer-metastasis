#Read copynumber csv
PQcorrected <- read.csv("C:/Research/Research_project/test/PQArmWgdCorrectAllCN.csv",header = TRUE,sep = ",")
#Transpose dataframe
PQc <- t(PQcorrected)
PQc <- as.data.frame(PQc)
#Name first row as column names
colm <- PQc[1,]
colnames(PQc) <- colm
#Remove first row
PQc <- PQc[-1,]

#Import list of whole genome doubled samples
wgd <- read.csv("C:/Research/Research_project/tcga/WGDsamplesAllData.csv",header = TRUE,sep = ";")

PQc["sampleid"] <- rownames(PQc)
test <- left_join(PQc,wgd,by="sampleid")
#Reoder
test1 <- test[,c(40:41,1:39)]

#Replace na with No string
library(tidyr)

test1$WGD <- test1$WGD %>% replace_na("No") 
Wgdarray <- test1$WGD

PQc <- subset(PQc,select = -c(40))

#Convert columns into numeric class
PQc[] <- lapply(PQc, function(x) type.convert(as.numeric(x)))

#Replace NAs of TCGA data with column mean
for(i in 1:ncol(PQc)){
  PQc[is.na(PQc[,i]), i] <- mean(PQc[,i], na.rm = TRUE)
}

#Remove zero 
PQc <- PQc[ , which(apply(PQc, 2, var) != 0)]


#Perform PCA
PQcPCA <- prcomp(PQc,scale. = TRUE,center = TRUE)

install.packages("factoextra","dendextend")
install.packages("ggplot2")
library("factoextra")

#Plot scree plot
fviz_eig(PQcPCA)

#Graph
fviz_pca_ind(PQcPCA,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

pca.var <- PQcPCA$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per

##Retrieve percentage variance and cumulative variance
summ = summary(PQcPCA)
summ_var = summ$importance[2,] * 100
barplot(summ_var[1:15])

#Create two-dimensional PCA score
Xscores = PQcPCA$x
par(mfrow = c(1,1))
my_cols <- c(5,7,8,10,12)
names(my_cols) <- c("Hartwig Liver","Hartwig Lung","TCGA Liver","TCGA Lung","TCGA Colon")
plot(Xscores[,1], Xscores[,2],pch=19,col=my_cols,xlab = "PC1",ylab = "PC2",main = "PCA plot")
legend("bottomleft",legend=c("Hartwig Liver","Hartwig Lung","TCGA Liver","TCGA Lung","TCGA Colon"), fill=c(5,7,8,10,12), border = NA, bty="n",x.intersp = 1,cex=0.7, title.adj = "adj")

#Create dataframe with X scores
x <- cbind(PQ1,Xscores[,1:2])
x$biopsysite <- c(rep("Hartwig Liver",327),rep("Hartwig Lung",27),rep("TCGA Liver",362),rep("TCGA Lung",950),rep("TCGA Colon",290))
x$WGD <- Wgdarray
dfmain <- x[,c(40:43,1:39)]

#salpie1 <- subset(dfmain,select = -c(1,2))
#write.csv(salpie1,"C:/Research/Research_project/test/dfwithoutWGD.csv",row.names = FALSE)

#Plot venn diagram
ggplot(dfmain,aes(PC1,PC2,fill=biopsysite,colour=biopsysite))+
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = " "))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = " "))+
  ggtitle("PCA graph : Primary and Metastatic samples:P & Q arms(WGD corrected)")+
  stat_ellipse(geom = "polygon",col="black",alpha = 0.2)+
  geom_point(pch = 21, colour = "black",size=1.5)


