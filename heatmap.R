###################################################
#         Prj: SAGL1
#         Assignment: heatmap
#         Author: Bin Zhao
#         Date: 2022.1.20
#         Locate: HENU                            
###################################################
setwd("C:/Users/Administrator/Desktop")
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library("openxlsx")
a<- read.xlsx("C:/Users/Administrator/Desktop/Table SFX-heatmap.xlsx",3)
a1<-a[-1,-1]
a3<-a[-1,]
a1<-apply(a1,2,as.numeric)
rownames(a1) <-a3$Name
##pheatmap
pheatmap(a1,scale="row",cluster_cols = FALSE,cluster_rows =TRUE,
         clustering_distance_rows="correlation",clustering_distance_cols="none",
         clustering_method="average",show_rownames = T,show_colnames=T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         treeheight_row=10,treeheight_col=5,cellwidth=6,cellheight=6
         ,fontsize=6,fontsize_row =4,fontsize_col =4,
         filename = "./o.pdf"
)