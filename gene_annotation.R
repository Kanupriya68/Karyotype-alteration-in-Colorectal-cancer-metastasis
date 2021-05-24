load("C:/Users/kanup/Downloads/gene_annotations.RData")

TCGA_data <- read.csv("C:/Research/Research_project/TCGA_data/TCGA_mastercalls.abs_segtabs.fixed.txt",header = TRUE,sep = "")

#Read gene annotation data and save it in a variable
results <- read.csv("C:/Research/Research_project/TCGA_data/results.csv",sep = ",",header = TRUE)

#Read all p-values of 6 combinations of contingency tables
we <- read.csv("C:/Research/Research_project/fishertestdf/allPvalues.csv",sep = ";",header = TRUE)
#Change column name 'chromosome' to 'chr'
colnames(we)[2] <- "chr"

#Change gene annotation variable name to gene_annotation
gene_annotation <- results
#Convert 'chr' column into character class
we$chr <- as.character(we$chr)
#Replace chr with chromosome name in p-value data frame
we$chr <- gsub("chr", "", we$chr)
#Create empty column named 'gene' with NAs in p-value data frame
we$gene <- "NA"
#Create for loop to annotate each bin row with gene names
for (i in 1:nrow(we)) {
  sub_we_chr <- subset(gene_annotation, we$chr[i] == gene_annotation$chromosome_name)
  subb <- subset(sub_we_chr, we$start[i] <= sub_we_chr$start_position & we$stop[i] >= sub_we_chr$end_position)
  list <- paste(c(unique(subb$hgnc_symbol)), collapse=", ")
  if (nrow(subb)<2) {
    list <- "no genes"
  }
  we$gene[i] <- list
}

#Write gene annotated data frame into a csv file
write.csv(we,file = "C:/Research/Research_project/fishertestdf/gene_annotation.csv")

#Filter gene_annotation dataframe as per significant bins
significantBins <- we %>% filter(gainloss_fisherPvalue < 0.05|cnLOHloss_fisherPvalue <0.05|
                      cnLOHgain_fisherPvalue < 0.05|cnLOHamp_fisherPvalue < 0.05|
                      amploss_fisherPvalue < 0.05 | ampgain_fisherPvalue < 0.05)
  
#Subset significantBins dataframe to have only bin locations and genes
geneSet <- subset(significantBins, select = c(1,2,3,4,17))

#Create csv file with all significant bins of all 6 combinations of contingency tables
write.csv(significantBins,file = "C:/Research/Research_project/fishertestdf/significant_bins.csv")
#Create csv file of geneSet dataframe
write.csv(geneSet,file = "C:/Research/Research_project/fishertestdf/geneSet.csv")

#Split geneSet daraframe gene name wise in each row
splitGenes <- separate_rows(geneSet,gene,convert = TRUE)
#Remove bins with 'no' and 'genes' as it has no genes
splitGenes = filter(splitGenes, gene != "no" & gene != "genes" & gene != "")
colnames(splitGenes)[1] <- "bn"

#Read copynumber dataframe with chromosome column
chr_data <- read.csv(file = "C:/Research/Research_project/test/dataCN.csv", header = TRUE, sep =",")

#Merge geneSet data and copynumber data
gseaDf <- left_join(splitGenes,chr_data,by = "bn" )
#Remove duplicate columns
gseaDf <- subset(gseaDf,select = -c(1:4))

#Extract above in csv
write.csv(gseaDf,file = "C:/Research/Research_project/test/gseaDf.csv",row.names = FALSE)

gseaDf <-  subset(gseaDf, select = -c(2:5))
#Clean column names
for ( col in 2:ncol(gseaDf)){
  colnames(gseaDf)[col] <-  sub("*X.", "", colnames(gseaDf)[col])
}

write.csv(df, file = "C:/Research/Research_project/binnedCnDataframe/biopsydf.csv", row.names = FALSE)

#Import metadata csv into a dataframe
biopsydf <- read.csv("C:/Research/Research_project/binnedCnDataframe/biopsydf.csv", sep = ",",header = TRUE)
#Select biopsy location list and sampleid only
biopsydf <- subset(biopsydf,select = c(1,2))

#transpose copy number data frame
sort <- t(gseaDf)
#Convert matrix into dataframe again
sort <- as.data.frame(sort)
#Make first row as column names (gene names)
colnames(sort) <- sort[1,]
#Delete first row 
sort=sort[-1,]
#create sampleid column and insert sampleid names into it
sort$sampleid <- row.names(sort)
#move newly added sampleid column from last column position to fisrt position
sort <- sort[,c(ncol(sort),1:(ncol(sort)-1))]
#join sort data with biopsy location dataframe by 'sampleid' common column
sorteddf <- left_join(sort,biopsydf,by="sampleid")
#Sort column positions again and tidy-up
sorteddf <- sorteddf[,c(ncol(sorteddf),1:(ncol(sorteddf)-1))]

#Sort only Lung Liver sample
sorteddf <- filter(sorteddf, biopsySite == "Liver" | biopsySite == "Lung")
sorteddf <-  arrange(sorteddf, biopsySite)
#Create phenotype array for GSEA 
biopsySiteAnno <- subset(sorteddf,select = c(1))
biopsySiteAnno <- t(biopsySiteAnno)
biopsySiteAnno <- as.data.frame(biopsySiteAnno)
firstRow <- c("354 2 1",rep("",353))
secondRow <- c("# Liver Lung",rep("",353))
biopsySiteAnno <- rbind(firstRow,secondRow,biopsySiteAnno)

write.csv(biopsySiteAnno, file = "C:/Research/Research_project/test/biopsyanno.csv",row.names = FALSE)

write.csv(sorteddf, file = "C:/Research/Research_project/test/significantBinswithGenes.csv",row.names = FALSE)

gseaDataset <- subset(sorteddf, select = -c(1))
gseaDataset <- t(gseaDataset)
colnames(gseaDataset) <- gseaDataset[1,]
gseaDataset=gseaDataset[-1,]

write.csv(gseaDataset, file = "C:/Research/Research_project/test/GseaDataSet.csv")

