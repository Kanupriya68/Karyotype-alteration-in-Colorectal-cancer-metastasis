##Create Heatmap of copynumbers and 633 patient sample
pdf("Figure1_A.pdf", width=10, height=10)

##Load libraries
x<-c("openxlsx", "tidyverse", "survminer", "ComplexHeatmap", "RColorBrewer", "ggbeeswarm", "MASS", "grid", "gtable", "gsubfn", "proto", "circlize", "xlsx", "RColorBrewer", "ComplexHeatmap", "grid", "data.table", "dplyr", "survminer", "ggpubr", "magrittr", "survival", "openxlsx", "reshape2", "reshape", "copynumber", "BiocGenerics", "parallel", "ggplot2", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
lapply(x, require, character.only = TRUE)
library(methods)

##Read sampledata csv and convert into a dataframe
CN_csv <- read.csv(file = "C:/Research/Research_project/heatmap_df.csv", sep = ",", header = TRUE)

CNmat <- as.matrix(heatmap_df[,-1])
chr <- heatmap_df$chromosome
rownames(CNmat) <- chr

#Heatmap(CNmat)
colFun = colorRamp2(c(-1, 3 , 6), c("blue", "white", "red"))
colFun(seq(-3, 3))
Heatmap(t(CNmat), name = "CN",cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = TRUE,col = colFun, na_col = "grey",
        column_names_side = "bottom", row_names_side = "right", width = unit(10,"cm"),
        height = unit(10, "cm"), heatmap_legend_param = list(at = c(-2,1,4,7,10), labels = c("superloss","loss", "diploid", "gain","amp")))
        

df = data.frame(type = c(rep("Adenomas", ncol(Ds)), rep("Ca-in-Ads", ncol(ca_in_ads)), rep("Carcinomas", ncol(CRCs))))

type = c("#ea7e5d", "#5deac5", "#5d83ea")
names(type) = c("Adenomas", "Ca-in-Ads", "Carcinomas")
colors = list(type = type)

ha <- HeatmapAnnotation(df = df, which = "row", col=colors, width = unit(1, "cm"))

draw(ha + CNs, row_hclust_side = "left", row_sub_title_side = "left", gap = unit(0.5, "cm"))

