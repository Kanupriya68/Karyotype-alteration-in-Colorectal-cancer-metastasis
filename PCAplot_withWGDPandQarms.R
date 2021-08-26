#Read copynumber csv
PQ <- read.csv("C:/Research/Research_project/tcga/PQarmCNdf.csv",header = TRUE,sep = ",")
#Transpose dataframe
PQ1 <- t(PQ)
PQ1 <- as.data.frame(PQ1)
#Name first row as column names
colm <- PQ1[1,]
colnames(PQ1) <- colm
#Remove first row
PQ1 <- PQ1[-1,]

#Import list of whole genome doubled samples
wgd <- read.csv("C:/Research/Research_project/tcga/WGDsamplesAllData.csv",header = TRUE,sep = ";")
wgd$sampleid <- gsub('-',".",wgd$sampleid)
PQ1["sampleid"] <- rownames(PQ1)

library(dplyr)
test <- left_join(PQ1,wgd,by="sampleid")
#Reoder
test1 <- test[,c(40,41,1:39)]
#Replace na with No string
library(tidyr)
test1$WGD <- test1$WGD %>% replace_na("No") 
samplesource <- c(rep("Hartwig Liver",327),rep("Hartwig Lung",27),rep("TCGA Liver",362),
                  rep("TCGA Lung",950),rep("TCGA Colon",292))
test1$SampleSource <- samplesource
test2 <- test1 %>% relocate(SampleSource,.after=sampleid)
test3 <- subset(test2,select = c(1,2))
write.csv(test3,file = "C:/Research/Research_project/test/samplesource.csv",row.names = FALSE)
Wgdarray <- test1$WGD
PQ1 <- subset(PQ1,select = -c(40))

#Convert columns into numeric class
PQ1[] <- lapply(PQ1, function(x) type.convert(as.numeric(x)))

#Remove zero 
#PQ1 <- PQ1[ , which(apply(PQ1, 2, var) != 0)]


#Perform PCA
PQpca <- prcomp(PQ1,scale. = TRUE,center = TRUE)

install.packages("factoextra","dendextend")
install.packages("ggplot2")
library("factoextra")

#Plot scree plot
fviz_eig(PQpca)

#Graph
fviz_pca_ind(PQpca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

pca.var <- PQpca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per

##Retrieve percentage variance and cumulative variance
summ = summary(PQpca)
summ_var = summ$importance[2,] * 100
barplot(summ_var[1:15])

#Create two-dimensional PCA score
Xscores = PQpca$x
par(mfrow = c(1,1))
my_cols <- c(5,7,8,10,12)
names(my_cols) <- c("Hartwig Liver","Hartwig Lung","TCGA Liver","TCGA Lung","TCGA Colon")
plot(Xscores[,1], Xscores[,2],pch=19,col=my_cols,xlab = "PC1",ylab = "PC2",main = "PCA plot")
legend("topright",legend=c("Hartwig Liver","Hartwig Lung","TCGA Liver","TCGA Lung","TCGA Colon"), fill=c(5,7,8,10,12), border = NA, bty="n",x.intersp = 1,cex=0.7, title.adj = "adj")

#Create dataframe with X scores
x <- cbind(PQ1,Xscores[,1:2])
x$biopsysite <- c(rep("Hartwig Liver",327),rep("Hartwig Lung",27),rep("TCGA Liver",362),rep("TCGA Lung",950),rep("TCGA Colon",290))
Wgdarray <- as.data.frame(Wgdarray)
Wgdarray <- Wgdarray[-c(1958,1956),]
x$WGD <- Wgdarray
dfmain <- x[,c(40:43,1:39)]

salpie <- subset(dfmain,select = -c(1,2))
write.csv(salpie,"C:/Research/Research_project/test/dfwithWGD.csv",row.names = FALSE)

#Plot venn diagram
ggplot(dfmain,aes(PC1,PC2,fill=biopsysite,colour=biopsysite))+
  xlab(paste("PC1 - ",pca.var.per[1],"%",sep = " "))+
  ylab(paste("PC2 - ",pca.var.per[2],"%",sep = " "))+
  ggtitle("PCA graph : Primary and Metastatic samples:P & Q arms(WGD)")+
  stat_ellipse(geom = "polygon",col="black",alpha = 0.2)+
  geom_point(pch = 19, colour = "black",size=0.5)+
  geom_point(aes(shape=WGD))

ggplot(dfmain, aes(fill=WGD, y=value, x=biopsysite)) + 
  geom_bar(position="stack", stat="identity")  

specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

