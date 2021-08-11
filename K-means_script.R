library(dplyr)
library(data.table)
library(biomaRt)
library(ggplot2)
library(reshape)
library(reshape2)
library(openxlsx)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(ggbeeswarm)
x<-c("gsubfn", "proto", "circlize", "xlsx", "RColorBrewer", "ComplexHeatmap", "grid", "data.table", "dplyr", "survminer", "ggpubr", "magrittr", "survival", "openxlsx", "reshape2", "reshape", "copynumber", "BiocGenerics", "parallel", "ggplot2", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
lapply(x, require, character.only = TRUE)
library(tidyverse)
#kanu = read.table("/Users/nowins01/Downloads/PQArmWgdCorrectAllCN.csv", header=T, sep=",")
kanu = read.table("C:/Users/kanup/Downloads/PQarmCNdfFinal.csv", header=T, sep=",")
ann = read.table("C:/Users/kanup/Downloads/dfwithoutWGD.csv", header=T, sep=",")
kanu_location_colour = c("black", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
ann = ann[,1:3]
cnas = kanu[,2:ncol(kanu)]
cnas = round(cnas)
cohort_bin_mat = cnas
row.names(cohort_bin_mat) <- kanu$autosomes

save(cohort_bin_mat, file="KANU_noWGD.RData")
df = as.data.frame(ann$biopsysite)
names(df) = "location"
dff = as.data.frame(ann$WholeGenomeDoubling)
names(dff) = "GenomeDoubling"
location = c("black", "#E2D200", "#46ACC8", "#E58601", "#B40F20")
names(location) = unique(ann$biopsysite)
colors = list(location = location)
GenomeDoubling = c("black", "white")
names(GenomeDoubling) = unique(ann$WholeGenomeDoubling)
colorss = list(location = GenomeDoubling)
rownames(cohort_bin_mat) = kanu$autosomes
pdf("Copynumber_heatmap_Kanu_clustered_final.pdf",width=20, height=10)
CopyNumber <- Heatmap(t(cohort_bin_mat),
                      name = "CopyNumber",
                      cluster_rows = TRUE, 
                      cluster_columns = FALSE,
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      height = 300,
                      row_dend_width = unit(6, "cm"),
                      row_km = 5,
                      # row_split = as.character(df$type),
                      col = colorRamp2(c(1, 2, 3,4,5,6,7), c("blue", "white", "pink","red", "orange", "darkorange", "purple")),
                      heatmap_legend_param = list(at = c(1, 2, 3,4,5,6,7), labels = c("CN1", "CN2", "CN3", "CN4","CN5", "CN6", "CN7")))
ha <- HeatmapAnnotation(df = df, which = "row", col=colors, simple_anno_size = unit(2, "cm"), width = unit(4, "cm"))
haa <- HeatmapAnnotation(df = dff, which = "row", col=colorss, simple_anno_size = unit(2, "cm"), width = unit(4, "cm"))
draw(ha+haa + CopyNumber, gap = unit(0.5, "cm"))
dev.off()
CopyNumber = draw(CopyNumber)
kmean_cluster <- row_order(CopyNumber)
lapply(kmean_cluster, function(x) length(x))
#extract sample id
mat <- t(cohort_bin_mat)
for (i in 1:length(kmean_cluster)){
   if (i == 1) {
     clu <- t(t(row.names(mat[kmean_cluster[[i]],])))
     out <- cbind(clu, paste("cluster", i, sep=""))
     colnames(out) <- c("sampleID", "Cluster")
     } else {
       clu <- t(t(row.names(mat[kmean_cluster[[i]],])))
       clu <- cbind(clu, paste("cluster", i, sep=""))
       out <- rbind(out, clu)
       }
}

#export
write.csv(out, file= "C:/Research/Research_project/test/kmean_clusters.csv", row.names=FALSE)
