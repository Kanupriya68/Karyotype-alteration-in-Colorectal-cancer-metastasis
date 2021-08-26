##Create Heatmap of copynumbers and 633 patient sample
pdf("Figure1_A.pdf", width=10, height=10)

##Load libraries
x<-c("openxlsx", "tidyverse", "survminer", "ComplexHeatmap", "RColorBrewer", "ggbeeswarm", "MASS", "grid", "gtable", "gsubfn", "proto", "circlize", "xlsx", "RColorBrewer", "ComplexHeatmap", "grid", "data.table", "dplyr", "survminer", "ggpubr", "magrittr", "survival", "openxlsx", "reshape2", "reshape", "copynumber", "BiocGenerics", "parallel", "ggplot2", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
lapply(x, require, character.only = TRUE)
library(methods)

##Read sampledata csv and convert into a dataframe
CN_csv <- read.csv(file = "C:/Research/Research_project/heatmap_df.csv", sep = ",", header = TRUE)

BiopsySite <- read.csv(file = "C:/Research/Research_project/biopsy_metadata.csv",sep = ";",header = TRUE)
colnames(BiopsySite)[1] <- "sampleid"
sites <- BiopsySite$biopsySite
sites <- as.data.frame(sites)
BiopsySites <- sites[-c(424:575),]

write.csv(heatmap_df,file = "C:/Research/Research_project/heatmap_df.csv",row.names = FALSE)

hmdf <- left_join(heatmap_df,BiopsySite,by="sampleid")
hmdf <- hmdf%>% relocate(biopsySite,.after = "sampleid")

#Remove NULL samples
hm <- hmdf[-c(424:575),-c(1:2)]
hm <- hm[,-c(2898:3112)]
CNmat <- as.matrix(hm)
#chr <- heatmap_df$chromosome
idlist <- hmdf$sampleid



#Heatmap(CNmat)
library("circlize")
colFun = colorRamp2(c(-1, 3 , 6), c("blue", "white", "red"))
colFun(seq(-3, 3))



CNs <- Heatmap(CNmat, name = "CN",cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,col = colFun, na_col = "grey",
        column_names_side = "bottom", row_names_side = "right", width = unit(15,"cm"),
        height = unit(15, "cm"), heatmap_legend_param = list(color_bar = "discrete",
        labels = c("homozygous loss","loss", "diploid", "gain","amplification"),border="black"),
        grid.lines(sum_chr_ends,c(-7,5),gp = gpar(lty = "dotted", lwd = 1)))
       

df = data.frame(BiopsySites)
df$BiopsySites <- gsub("dorsal from vena cava inferior/medial right Kidney","Kidney",df$BiopsySites)



type = c("#ea7e5d", "#5deac5", "#5d83ea","#66d9ff",
         "#660066","#d9d9d9","#006600","#ffb3ff",
         "#006666","#a6ff4d","#e600ac","#cccc00",
         "#ff4d4d","#ccccff","#ff0066","#ffff00",
         "#ff6600","#4c0080","#ff4000","#ff9900",
         "#333300","#ffffcc","#b3ff66","#ccccff",
         "#FF0000","#990033","#ffe6ee","#e0b3ff",
         "#1f004d","#ffcc99","#00ffff","#9cb200",
         "#b78d91","#ffcc15","#ff3315","#623315")

names(type) = c("Abdomen", "abdominal", "Abdominal metastatic",
                "Abdominal wall","Abdominal wall right","adnex",
                "Adrenal gland","belly button","Bone","Breast",
                "Colon","cutane","Cutaneous",
                "ilium","Intra abdominal mass","Kidney","Liver",
                "Lung","Lymph node","Mamma","neck right","Omental cake",
                "Omentum","peritoneale","Peritoneum","Primary","rectum",
                "sigmoid","skin","Soft tissue",
                "sternum","subcutane laesion","Subcutaneous",
                "Thoracic wall","upper abdomen","Vagina")

colors = list(type = type)

ha <- HeatmapAnnotation(df = df, which = "row", col=colors, width = unit(0.5, "cm"))

draw(ha + CNs, row_sub_title_side = "left", gap = unit(0.5, "cm"))


dev.off()





