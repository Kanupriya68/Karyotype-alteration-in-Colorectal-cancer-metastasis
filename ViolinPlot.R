

wgdCorrectDf <- read.csv("C:/Research/Research_project/test/allCNWGDCorrected.csv")

newdata <- wgdCorrectDf[c(-1,-2,-3,-4)]


dfg <- data.frame(id=numeric(),abbrFrac=numeric())


for (i in 1:ncol(newdata)) {
  samp = newdata[,i]
  abbrFraction <- (length(samp[(samp<=1.5)]) + length(samp[(samp>=2.5)]))/length(samp)
  de <- data.frame(i,abbrFraction)
  names(de)<-c("id","abbrFrac")
  dfg <- rbind(dfg, de)
  
}

#Violin Plot for WGDcorrected
library(tidyverse)
library(plotly)
library(IRdisplay)
dfg %>% ggplot(aes(x="Whole Genome Corrected Violin Plot",y=abbrFrac))+
  geom_violin() + geom_boxplot(width=0.1)+
  theme(axis.title.x =  element_blank())


wgdDf <- read.csv("C:/Research/Research_project/test/Alldata.csv")

ndata <- wgdDf[c(-1,-2,-3)]


dfge <- data.frame(id=numeric(),abbrFrac=numeric())


for (i in 1:ncol(ndata)) {
  samp = ndata[,i]
  abbrFraction <- (length(samp[(samp<=1.5)]) + length(samp[(samp>=2.5)]))/length(samp)
  de1 <- data.frame(i,abbrFraction)
  names(de1)<-c("id","abbrFrac")
  dfge <- rbind(dfge, de1)
  
}

#Violin Plot for WGD data

dfge %>% ggplot(aes(x="Whole Genome doubled Violin Plot",y=abbrFrac))+
  geom_violin(draw_quantiles = TRUE) + geom_boxplot(width=0.1)+
  theme(axis.title.x =  element_blank() )+
  scale_fill_manual(values = "#F5FCC2")

p <- dfge %>% plot_ly(x="",y=~abbrFrac, type = "violin",
                      box=list(visible=TRUE,width=0.2))
htmlwidgets::saveWidget(p,"p.html")
display_html('<iframe src="p.html" width=600 height=600 frameborder="0"></iframe>')



biopsySiteDF <- read.csv("C:/Research/Research_project/test/samplesource.csv",header = TRUE)

biopsySiteDF$SampleSource[biopsySiteDF$SampleSource == 'Hartwig Liver'] <- 'Metastatic Liver'
biopsySiteDF$SampleSource[biopsySiteDF$SampleSource == 'Hartwig Lung'] <- 'Metastatic Lung'

biopsySiteDF$SampleSource[biopsySiteDF$SampleSource == 'TCGA Liver'] <- 'Primary Liver'
biopsySiteDF$SampleSource[biopsySiteDF$SampleSource == 'TCGA Lung'] <- 'Primary Lung'
biopsySiteDF$SampleSource[biopsySiteDF$SampleSource == 'TCGA Colon'] <- 'Primary Colon'


biopsySiteDF$SampleSource <- as.factor(biopsySiteDF$SampleSource)

wgdCorrectDf <- read.csv("C:/Research/Research_project/test/allCNWGDCorrected.csv")

newdata <- wgdCorrectDf[c(-1,-2,-3,-4)]

colors <- c("#FFFFFF","#F5FCC2")


dfg <- data.frame(sampleid=numeric(),AbberationFraction=numeric())


for (i in colnames(newdata)) {
  samp = newdata[,i]
  samp <- na.omit(samp)
  abbrFraction <- (length(samp[(samp<=1.5)]) + length(samp[(samp>=2.5)]))/length(samp)
  de <- data.frame(i,abbrFraction)
  names(de)<-c("sampleid","AbberationFraction")
  dfg <- rbind(dfg, de)
  
}

wgdCorretedDF <- dfg

wgdCorretedDF <- left_join(wgdCorretedDF, biopsySiteDF, by = "sampleid")




wgdCorretedDF %>% ggplot(aes(x=SampleSource,AbberationFraction, fill=SampleSource))+
  geom_violin(show.legend = FALSE) + geom_boxplot(width=0.1)+
  theme(axis.title.x =  element_blank() )



wgdDf <- read.csv("C:/Research/Research_project/test/Alldata.csv")

ndata <- wgdDf[c(-1,-2,-3)]


dfge <- data.frame(sampleid=numeric(),AbberationFraction=numeric())


for (i in colnames(ndata)) {
  samp = ndata[,i]
  samp <- na.omit(samp)
  abbrFraction <- (length(samp[(samp<=1.5)]) + length(samp[(samp>=2.5)]))/length(samp)
  de1 <- data.frame(i,abbrFraction)
  names(de1)<-c("sampleid","AbberationFraction")
  dfge <- rbind(dfge, de1)
  
}

WgDoubledDF <- dfge

WgDoubledDF <- left_join(WgDoubledDF, biopsySiteDF, by = "sampleid")



WgDoubledDF %>% ggplot(aes(x=SampleSource,AbberationFraction, fill=SampleSource))+
  geom_violin(show.legend = FALSE) + geom_boxplot(width=0.1)+
  theme(axis.title.x =  element_blank() )



combined <- rbind(dfg,dfge)


#combined %>% ggplot(aes(x=Whole_Genome_Doubling,AbberationFraction, fill=Whole_Genome_Doubling))+
  #geom_violin(show.legend = FALSE) + geom_boxplot(width=0.1)+
  #theme(axis.title.x =  element_blank() )


anova1 <- aov(WgDoubledDF$AbberationFraction~WgDoubledDF$SampleSource)
summary(anova1)


TukeyHSD(anova1)
plot(TukeyHSD(anova1),las=1)


anova2 <- aov(wgdCorretedDF$AbberationFraction~wgdCorretedDF$SampleSource)


summary(anova2)


TukeyHSD(anova2)
plot(TukeyHSD(anova2),las=1)




