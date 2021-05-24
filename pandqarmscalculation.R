pandq <- p_and_q
acrocentric <- p_and_q


#cyto <- read.table("/Users/nowins01/Documents/PCa/hg19cytoBandUCSC.txt")

load("p_and_q_arms_hg19.RData")

#remove acrocentric arms 13,14,15,21,22
load("p_and_q_arms_hg19_no_acrocentric.RData")


##load CN data
dataCN <-  read.csv("C:/Research/Research_project/test/CNdata.csv", header = TRUE,sep = ",")

dataCN$arm = "p"
library(dplyr)
dataCN <- dataCN %>% relocate(arm, .after = stop)


for (i in 1:nrow(pandq)) {
  dataCN$arm[(dataCN$chromosome == pandq$Chromosome[i] & dataCN$start >= pandq$Start[i] & dataCN$stop <= pandq$End[i])] <- pandq$cyto[i]
  dataCN$arm[(dataCN$chromosome == pandq$Chromosome[i] & dataCN$start >= pandq$Start[i] & dataCN$stop >= pandq$End[i])] <- "q"
    
} 

######################################################
#Aberration Fraction of Hartwig data for p and q arms#
######################################################

#Create categorical variables with defined range of CN aberrations
datacn1 = dataCN[,c(7:639)] > 1.5 & dataCN[,c(7:639)] < 2.5
stick <- subset(dataCN, select = c(1:6))
df <- cbind(stick,datacn1)
df <- df[-c(2898:3113), ]


#Create base dataframe with chromosome arm id
firstdf <- data.frame(chromosome=c(1:22,1:22),arm=c(rep("p",22),rep("q",22)))
firstdf$id<- paste(firstdf$chromosome,",",firstdf$arm)
merged_df <- subset(firstdf,select = c(id))


#Loop to calculate and append further sample fractions

j=6
for (i in colnames(df[7:639])){
  j=j+1
  dfy <- subset(df, select = c("chromosome","arm",i))
  
  
  print(dfy)
  
  if(length(unique(dfy[,3]))<2){
    print(paste("All TRUE in sample: ",i ))
    skipList <- c(skipList, i)
    
    truedf <- data.frame(id = merged_df[1], i = 1)
    names(truedf)[names(truedf) == "i"] <- i
    merged_df <- left_join(merged_df,truedf, by = "id")
    
    
    
  }
  else{
    
    test <- dfy%>%
      group_by(chromosome,arm,df[j])%>%
      summarise(count=n())%>%
      spread(key = 3, value = count, fill = 0)%>%
      magrittr::set_colnames(c("chromosome","arm","total_not","total_yes"))
    print("CheckA")
    print(class(test))
    test <- as.data.frame(test)
    #if (ncol(test) < 4){
      #test$tru <- 0
      
    #}
    test$fraction <- round(test[4]/(test[4]+test[3]),digits = 4)
    names(test)[names(test) == "fraction"] <- i
    print(i)
    test$id <- paste(test$chromosome,",",test$arm)
    print("check1")
    colnames(test)
    print(test)
    #test <- subset(test, select = c("id",i))
    
    
    test <- as.data.frame(test[,c("id",i)])
    print("check2")
    
    print(i)
    print("print")
    test <- data.frame(id = test[1][[1]], i = test[2][[1]])
    names(test)[names(test) == "total_yes"] <- i
    
    
    #if (j<8){
    merged_df <- left_join(merged_df,test, by = "id")
      
    #}
#    else if(i=="CPCT02020265T"){
#      print("inside else if")
#      print(i)  }
    #else{
      #merged_df <- left_join(merged_df,test, by = "id")
      
    #}
    
    
  }
}
View(merged_df)

library(stringr)
pqdf <- str_split_fixed(merged_df$id, ",",2)
pqdf <- as.data.frame(pqdf)
merged_df$Chromosome <- pqdf$V1
merged_df$Arms <- pqdf$V2
merged_df <- merged_df %>% relocate(Chromosome, .after = id)%>%
  relocate(Arms, .after = Chromosome)
merged_df <- subset(merged_df,select = -c(1))

#Remove X and Y chromosome rows
#merged_df <- merged_df[-c(45,46), ]

#Arrange chromosomes as per order
merged_df <- merged_df[order(as.integer(merged_df$Chromosome),decreasing = FALSE), ]

#Remove .TRUE from patient sample column names
#for ( col in 3:ncol(merged_df)){
#colnames(merged_df)[col] <-  sub("$'TRUE'*", "", colnames(merged_df)[col])
#}


write.table(merged_df, file = "C:/Research/pqarms/pqarmsHartwig.csv",sep = ",",row.names = FALSE)



######################################################
####################Histogram codes###################
######################################################

pandqArmsCN <- subset(dataCN, select = -c(1:6))

#write.csv(pandqArmsCN,file = "C:/Research/Research_project/test/pandqarmsCN.csv",row.names = FALSE)

p_arms_per_chrom_loss = list()
p_arms_per_chrom_gain = list()
q_arms_per_chrom_loss = list()
q_arms_per_chrom_gain = list()

for (i in 1:22) {
  if (i != 13 | i != 14 | i != 15 | i != 21 | i != 22) {
    p_arms = dataCN[(dataCN$chromosome == i & dataCN$arm == "p"),]
    
    #how much of p arm is lossed
    p_arms_per_chrom_loss_ = colSums(p_arms[,c(7:639)] < 2)/nrow(p_arms)
    #p_arms_per_chrom_loss_ = p_arms_per_chrom_loss_[p_arms_per_chrom_loss_ != 0]
    p_arms_per_chrom_loss = c(p_arms_per_chrom_loss, p_arms_per_chrom_loss_)
    
    #how much of p arm is gained
    p_arms_per_chrom_gain_ = colSums(p_arms[,c(7:639)] > 3 & p_arms[,c(7:639)] < 5 )/nrow(p_arms)
    p_arms_per_chrom_gain_ = p_arms_per_chrom_gain_[p_arms_per_chrom_gain_ != 0]
    p_arms_per_chrom_gain = c(p_arms_per_chrom_gain, p_arms_per_chrom_gain_)
  }
  
  q_arms = dataCN[(dataCN$chromosome == i & dataCN$arm == "q"),]
  
  #how much of q arm is lossed
  q_arms_per_chrom_loss_ = colSums(q_arms[,c(7:639)] < 2)/nrow(q_arms)
  #q_arms_per_chrom_loss_ = q_arms_per_chrom_loss_[q_arms_per_chrom_loss_ != 0]
  q_arms_per_chrom_loss = c(q_arms_per_chrom_loss, q_arms_per_chrom_loss_)
  
  #how much of q arm is gained
  q_arms_per_chrom_gain_ = colSums(q_arms[,c(7:639)] > 3 & q_arms[,c(7:639)] < 5)/nrow(q_arms)
  q_arms_per_chrom_gain_ = q_arms_per_chrom_gain_[q_arms_per_chrom_gain_ != 0]
  q_arms_per_chrom_gain = c(q_arms_per_chrom_gain, q_arms_per_chrom_gain_)
}


par(mfrow=c(2,2))

hist(as.numeric(p_arms_per_chrom_loss), main="amount of p arm lost", breaks=100, col="blue")

hist(as.numeric(q_arms_per_chrom_loss), main="amount of q arm lost", breaks=100, col="blue")

hist(as.numeric(p_arms_per_chrom_gain), main="amount of p arm gained", breaks=100, col="red")

hist(as.numeric(q_arms_per_chrom_gain), main="amount of q arm gained", breaks=100, col="red")
