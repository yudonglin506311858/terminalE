YDL<- readRDS("SINGLET_terminalE_reti.rds")
YDL
set.seed(123)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
#先运行几个降维算法之后再读取文件
###prepare normalized expression for Cellrouter input
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2"))
write.csv(mydata,file = "projection.csv")
write.table(as.matrix(YDL@assays$RNA@counts),"YDL.normalized_expression.txt",sep="\t")
write.csv(colnames(YDL@assays$RNA@counts),"YDL.cell_names.csv")
write.csv(rownames(YDL@assays$RNA@counts),"YDL.gene_names.csv")
write.table(as.matrix(YDL@meta.data),"YDL.meta.data.txt")




#run this script in linux system and oracle java
source('/data/yudonglin/nopro/Cellrouter/cellrouter/CellRouter_Class.R')
libdir <- '/data/yudonglin/nopro/Cellrouter/cellrouter/CellRouter/'
set.seed(123)
library(dplyr)
library(plotrix)
matrix=read.table("projection.csv",sep=",",header=T,row.names=1)
colnames(matrix) <- c('UMAP1','UMAP2')
#rownames(matrix) = colnames (matrix)
rownames(matrix)=gsub("-1","",rownames(matrix))
ndata <- read.table('YDL.normalized_expression.txt',sep="\t",header=T,row.names=1)
genes <-as.vector(rownames(ndata))
map <- data.frame(id=rownames(ndata),symbol=genes,stringsAsFactors = FALSE)
ndata <- averageIds(ndata,map,'symbol')
#Remove genes with zero variance across all cells
var <- apply(ndata,1,var)
var <- var[which(var > 0)]
ndata <- ndata[names(var),]
### selecting genes to use as regulated along developmental trajectories.
#pca <- prcomp(t(ndata),scale=TRUE,center=TRUE)
#loadings <- pca$rotation
#num_pc <- 5
#quantile <- 0.975
#genes2use <- unique(as.vector(unlist(apply(loadings[,1:num_pc],2,function(x){names(x[which(abs(x) >= quantile(x,quantile))])}))))
genes2use=rownames(ndata)
ggrn <- buildGRN('Mm',ndata,genes2use,2,'results/GRN.R') #original 5
rownames(matrix) = colnames(ndata)
#下次可以直接load这个GRN文件
#ggrn <- get(load('results/GRN.R'))
### Subpopulation identification and gene signatures with CellRouter
cellrouter <- CellRouter(expdata=ndata,annotations=colnames(ndata))
cellrouter@rdimension <- matrix
pdf("kNN_network.pdf")
cellrouter <- findsubpopulations(cellrouter,90,'jaccard','results/kNN_network.gml')
dev.off()
df=read.table("YDL.meta.data.txt",header=T,row.names=1,sep="\t")
df$sample_id=rownames(df)
df=merge(cellrouter@sampTab,df,by="sample_id",all=T)
write.table(df,"YDL.meta.data.withSP.txt",sep="\t")
lengths(cellrouter@graph$subpopulation)
cellrouter <- diffexpr(cellrouter,column='population',pvalue = 0.05)
markers <- findmarkers(cellrouter)
write.table(markers,"results/YDL.markers.txt",sep="\t")
plotReducedDimension(cellrouter,5,5,filename='results/YDL.tSNE.pdf')
table(cellrouter@sampTab$population)
write.table(cellrouter@sampTab,"results/YDL.cellrouter_sampTab.txt",sep="\t")
######## Trajectory Detection using CellRouter ###
pdf("kNN_network_trajectory.pdf")
cellrouter <- createKNN(cellrouter,90,'jaccard','results/paths/kNN_network_trajectory.gml') #10 before this 90
dev.off()
filename <- "results/paths/cell_edge_weighted_network.txt"
write.table(cellrouter@graph$edges,file=filename,sep='\t',row.names=FALSE,col.names = FALSE,quote=FALSE) #input network
saveRDS(cellrouter,"cellrouter.RDS")
dev.off()


