#Create categorical variables with defined range of CN aberrations
datacn1 = dataCN[,c(7:639)] > 3 & dataCN[,c(7:639)] < 5
stick <- subset(dataCN, select = c(1:6))
df <- cbind(stick,datacn1)

skipList = list()

j=6
for (i in colnames(df[7:639])){
  j=j+1
  dfy <- subset(df, select = c("chromosome","arm",i))
  
  
  print(dfy)
  
  if(length(unique(dfy[,3]))<2){
    print(paste("All False in sample: ",i ))
    skipList <- c(skipList, i)
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
  if (ncol(test) < 4){
    test$tru <- 0
    
  }
  test$fraction <- round(test[4]/(test[4]+test[3]),digits = 2)
  names(test)[names(test) == "fraction"] <- i
  print(i)
  test$id <- paste(test$chromosome,",",test$arm)
  print("check1")
  colnames(test)
  print(test)
  #test <- subset(test, select = c("id",i))
  
    
  test <- as.data.frame(test[,c("id",i)])
  print("check2")
  print(test)
  if (j<8){
    merged_df <- test
  }
  else if(i=="CPCT02020265T"){
    print("inside else if")
    print(i)  }
  else{
    merged_df <- left_join(merged_df,test, by = "id")
  
  }
  
  
  }
}
View(merged_df)

pqdf <- str_split_fixed(merged_df$id, ",",2)
pqdf <- as.data.frame(pqdf)
merged_df$Chromosome <- pqdf$V1
merged_df$Arms <- pqdf$V2
merged_df <- merged_df %>% relocate(Chromosome, .after = id)%>%
  relocate(Arms, .after = Chromosome)
merged_df <- subset(merged_df,select = -c(1))

#Remove X and Y chromosome rows
merged_df <- merged_df[-c(45,46), ]

#Arrange chromosomes as per order
merged_df <- merged_df[order(as.integer(merged_df$Chromosome),decreasing = FALSE), ]

#Remove .TRUE from patient sample column names
#for ( col in 3:ncol(merged_df)){
  #colnames(merged_df)[col] <-  sub("$'TRUE'*", "", colnames(merged_df)[col])
#}


write.table(merged_df, file = "C:/Research/pqarms/pqarmsWithFraction.csv",sep = ",",row.names = FALSE)





