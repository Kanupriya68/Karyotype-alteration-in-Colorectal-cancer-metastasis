#Create dataframe from mean table

# read file path
cnv_paths <-
  list.files(path = "C://Research/Research_project/meanCN",
             pattern = "*.tsv",
             full.names = TRUE)

# read file content
cnv_contents <-
  cnv_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = ",",
         encoding = "UTF-8")

# read file name
cnv_filenames <- cnv_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
cnv_lists <- mapply(c,cnv_contents, cnv_filenames, SIMPLIFY = FALSE)

all_result <- rbindlist(cnv_lists, fill = T)


dfhj <- subset(all_result, select = c(2,3,4))


# change column name
names(dfhj)[3] <- "File.Path"


#Pivot the dataframe for each column to have file name
wide_data <- pivot_wider(dfhj, names_from = File.Path, values_from = meanCopyNumber)
#View(wide_data)

#Remove extra row with NAs
final_CNdataframe <- wide_data[-c(25),]

#Rename columns as per sample name
colnames(final_CNdataframe) <- gsub('.purple.cnv.somatic.tsv','',colnames(final_CNdataframe))

#Round-off dataframe
rounded_df <- subset(final_CNdataframe,select=-c(1))

rounded_off <- round(rounded_df, digits = 2)
rounded_off$chromosome <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y") 
heatmap_df <-rounded_off[,c(634, 1:633)]


write.csv(final_CNdataframe, file = "C://Research/Research_project/heatmap_df.csv")


