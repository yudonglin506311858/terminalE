#合并小鼠和人类的单细胞转录组数据

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)

setwd("D:/analysis/forpublication/bulk_sc")

YDL<-readRDS("D:/analysis/forpublication/bulk_sc/terminalE.rds")
#saveRDS(YDL,"terminalE.rds")


DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5) + NoLegend()
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "tsne", pt.size = 1.5, label=TRUE)
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE)
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,group.by="celltype")

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,group.by="Phase")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident",group.by="Phase")
DimPlot(YDL, reduction = "tsne", pt.size = 1.5, label=TRUE,group.by="Phase")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident",group.by="Phase")


allmarker<-read.csv("allmarker_PURE.csv")
head(allmarker)
table(allmarker$cluster)


setwd("D:/analysis/forpublication/bulk_sc/int")
library(pheatmap)
cellInfo<- read.table("Cell.Info.txt", sep = "\t", header = T, row.names = 1)

celltype = subset(cellInfo,select = 'CellType')
head(celltype)
write.csv(celltype,"celltype.txt")
AUCmatrix<-read.table("AUCell.txt")
BINmatrix <- read.table("binary_mtx.txt", header = T)
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

colnames(AUCmatrix)<-celltype$CellType
library(gtools)
AUCmatrix<- AUCmatrix[,mixedorder(colnames(AUCmatrix))]
AUCmatrix[1:4,1:4]
dim(AUCmatrix)
table(colnames(AUCmatrix))


#计算每个cluster的score的平均值
E0<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E0")],1,mean) #方差
E1<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E1")],1,mean) #方差
E2<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E2")],1,mean) #方差
E3<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E3")],1,mean) #方差
E4<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E4")],1,mean) #方差
E5<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E5")],1,mean) #方差



data<-cbind(E0,E1,E2,E3,E4,E5)

data<-data[which(rowSums(data) > 0),]#用R去除全是0的行
pheatmap(data,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(data,fontsize = 7,scale = "row",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "single",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "average",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "complete",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "complete",show_rownames =T,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "complete",show_rownames =F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




rownames(data)<-gsub("\\_.*","",rownames(data))

library(Hmisc)
rownames(data)<-toupper(rownames(data))

write.csv(data,"mouse_regulon——E0-5.csv")



setwd("D:/analysis/forpublication/human_bm_ubc/all/SCENIC/int")
library(pheatmap)
cellInfo<- read.table("Cell.Info.txt", sep = "\t", header = T, row.names = 1)

celltype = subset(cellInfo,select = 'CellType')
head(celltype)
write.csv(celltype,"celltype.txt")
AUCmatrix<-read.table("AUCell.txt")
BINmatrix <- read.table("binary_mtx.txt", header = T)
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

colnames(AUCmatrix)<-celltype$CellType
library(gtools)
AUCmatrix<- AUCmatrix[,mixedorder(colnames(AUCmatrix))]
AUCmatrix[1:4,1:4]
dim(AUCmatrix)
table(colnames(AUCmatrix))


#计算每个cluster的score的平均值
E0<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("ProE/BasoE")],1,mean) #方差
E1<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("Early-PolyE")],1,mean) #方差
E2<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("Late-PolyE")],1,mean) #方差
E3<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("Early-OrthoE")],1,mean) #方差
E4<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("Late-OrthoE")],1,mean) #方差


data<-cbind(E0,E1,E2,E3,E4)
colnames(data)<-c("ProE/BasoE","Early-PolyE","Late-PolyE","Early-OrthoE","Late-OrthoE")
data<-data[which(rowSums(data) > 0),]#用R去除全是0的行
pheatmap(data,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(data,fontsize = 7,scale = "row",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "single",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "average",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "complete",show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "complete",show_rownames =T,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "complete",show_rownames =F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




rownames(data)<-gsub("\\_.*","",rownames(data))

library(Hmisc)
rownames(data)<-toupper(rownames(data))

write.csv(data,"human_regulon——E0-5.csv")


setwd("D:/analysis/forpublication/summary/human_mouse")
mouse<-read.csv("D:/analysis/forpublication/bulk_sc/int/mouse_regulon——E0-5.csv")
colnames(mouse)<-c("GENE","MOUSE_E0","MOUSE_E1","MOUSE_E2","MOUSE_E3","MOUSE_E4","MOUSE_E5")
mouse$GENE<-toupper(mouse$GENE)

human<-read.csv("D:/analysis/forpublication/human_bm_ubc/all/SCENIC/int/human_regulon——E0-5.csv")
colnames(human)<-c("GENE","HUMAN_ProE/BasoE","HUMAN_Early-PolyE","HUMAN_Late-PolyE","HUMAN_Early-OrthoE","HUMAN_Late-OrthoE")
all<-merge(mouse,human,by=c("GENE"))


mouse<-mouse[!duplicated(mouse$GENE),] #删掉所有列上都重复的


human<-human[!duplicated(human$GENE),] #删掉所有列上都重复的



merge_tf0<-mouse$GENE
merge_tf1<-human$GENE


library (ggvenn)
x<-list("mouse regulons"=merge_tf0,
        "human regulons"=merge_tf1)
ggvenn(x)



#合并不同的矩阵，列名弄成一样长，然后补上0值。

ID<-c(mouse$GENE,human$GENE)
name<-as.data.frame(ID)
name<-name[!duplicated(name$ID),] #删掉所有列上都重复的
name<-as.data.frame(name)
colnames(name)<-c("GENE")
head(name)



DATA_cluster0<-dplyr::bind_rows(name,mouse)
DATA_cluster0<-DATA_cluster0[order(DATA_cluster0[,2],decreasing=T),]#以第一列降序排列
DATA_cluster0<-DATA_cluster0[!duplicated(DATA_cluster0$GENE),] #删掉所有列上都重复的
DATA_cluster0[is.na(DATA_cluster0)] <- 0#把NA值全部替换为0
#colnames(DATA_cluster0)<-c("GENE","NONPRO_CLUSTER0","NONPRO_CLUSTER1","NONPRO_CLUSTER2","NONPRO_CLUSTER3","NONPRO_CLUSTER4","NONPRO_CLUSTER5","NONPRO_CLUSTER6")

head(DATA_cluster0)
#DATA_cluster0[c(2:7)]<-scale(DATA_cluster0[,c(2:7)],center = T,scale = T)

DATA_cluster1<-dplyr::bind_rows(name,human)
DATA_cluster1<-DATA_cluster1[order(DATA_cluster1[,2],decreasing=T),]#以第一列降序排列
DATA_cluster1<-DATA_cluster1[!duplicated(DATA_cluster1$GENE),] #删掉所有列上都重复的
DATA_cluster1[is.na(DATA_cluster1)] <- 0#把NA值全部替换为0
#colnames(DATA_cluster1)<-c("GENE","HUMAN_CLUSTER0","HUMAN_CLUSTER1","HUMAN_CLUSTER2","HUMAN_CLUSTER3","HUMAN_CLUSTER4","HUMAN_CLUSTER5","HUMAN_CLUSTER6","HUMAN_CLUSTER7","HUMAN_CLUSTER8")

head(DATA_cluster1)
#DATA_cluster1[c(2:6)]<-scale(DATA_cluster1[,c(2:6)],center = T,scale = T)


DATA<-merge(DATA_cluster0,DATA_cluster1,by=c("GENE"))

# DATA<-DATA[,c("ID","cluster0","cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7")]
# 
# dim(DATA)
# head(DATA)
# DATA<-DATA[!duplicated(DATA$ID),] #删掉所有列上都重复的
rownames(DATA)<-DATA[,1]
DATA<-DATA[,-1]
DATA<-as.data.frame(DATA)
# https://blog.csdn.net/c1z2w3456789/article/details/79467095

head(DATA)



write.csv(DATA,"regulon打分.csv")
DATA<-read.csv("regulon打分.csv",row.names = 1)

data<-DATA

# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pheatmap(data,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))
pheatmap(data,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,scale = "row",fontsize = 7,show_rownames = F,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,scale = "row",fontsize = 7,show_rownames = T,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




pheatmap(data,scale = "row",fontsize = 6,filename = "new2.pdf",
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()
pheatmap(data,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,scale = "row",fontsize = 7,filename = "all_regulon.pdf",width = 10,height = 40,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,scale = "row",fontsize = 7,filename = "all_regulon.jpg",width = 10,height = 40,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




mouse_specific<-setdiff(merge_tf0,merge_tf1)
human_specific<-setdiff(merge_tf1,merge_tf0)
merge_gene<-intersect(merge_tf0,merge_tf1)


write.csv(mouse_specific,"mouse_specific.csv")
write.csv(human_specific,"human_specific.csv")
write.csv(merge_gene,"mouse_human_common.csv")
write.csv(newdata,"mouse_human_打分.csv")




#https://blog.csdn.net/c1z2w3456789/article/details/79467095

merge(mouse_specific,human_specific)


newdata<-data[c(merge_gene),]
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))


pheatmap(newdata,scale = "row",fontsize = 7,show_rownames = T,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))


pheatmap(newdata,scale = "row",
         fontsize = 7,show_rownames = T,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))



pheatmap(newdata,scale = "row",fontsize = 7,filename = "共有_regulon.pdf",width = 10,height = 12,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))


