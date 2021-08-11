#Read dataframe containing copy number calls with P and Q arm annotation
data=read.table("C:/Research/Research_project/tcga/PQarmCNdf.csv", header=T, sep=",")

#Round off copy number values from 2nd column
data[,2:ncol(data)] = round(data[,2:ncol(data)])
#Name all the columns
names = names(data)[2:ncol(data)]
#Convert into a dataframe
names = as.data.frame(names)
#Calculate mean of each column and filter greater than 3 values in new column
names$WGD = colMeans(data[,2:ncol(data)]) > 3
#Read csv where whole genome doubling is removed 
ann = read.table("C:/Research/Research_project/test/dfwithoutWGD.csv", header=T, sep=",")
#Use above csv information to annotate each biopsy location as per sample ID
names$ann = ann$biopsysite
#Subset
sub = names[,c(3,2)]
#Convert into table
tab = table(sub)
#Convert into dataframe
tab = as.data.frame(tab)
#Render using ggplot
install.packages("ggplot")
ggplot(tab, aes(fill=WGD, y=Freq, x=ann)) + 
  geom_bar(position="fill", stat="identity") +
  theme_cowplot(12) +
  xlab("")
