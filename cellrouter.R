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


##select starting subpopulation,all other subpopulations are targets
sources <- c('SP_2') #from SP_9 to SP_16
targets <- setdiff(as.vector(cellrouter@sampTab$population),sources)
methods <- c("euclidean","maximum","manhattan","canberra","binary",'graph') #graph for distances in KNN
cellrouter <- findpaths(cellrouter,libdir,paste(getwd(),'results/paths',sep='/'),method="graph")
ranks <- c('path_cost','path_flow','rank','length')
cellrouter <- processtrajectories(cellrouter,genes2use,path.rank=ranks[3],num.cells = 3,neighs = 1)
names <- unique(names(cellrouter@pathsinfo$distr))
clusters.show <- names
cellrouter <- correlationpseudotime(cellrouter,type='spearman')
cellrouter <- topgenes(cellrouter,0.85,0.15)
cellrouter <- smoothdynamics(cellrouter,names)
cellrouter <- clusterGenesPseudotime(cellrouter,10)
save(cellrouter,file='results/CellRouter_StemID_Processed.R')
saveRDS(cellrouter,"cellrouter_1.RDS")
cellrouter@signatures$SP_1$subpopulation="SP_1"
cellrouter@signatures$SP_2$subpopulation="SP_2"
cellrouter@signatures$SP_3$subpopulation="SP_3"
cellrouter@signatures$SP_4$subpopulation="SP_4"
cellrouter@signatures$SP_5$subpopulation="SP_5"
cellrouter@signatures$SP_6$subpopulation="SP_6"
cellrouter@signatures$SP_7$subpopulation="SP_7"
cellrouter@signatures$SP_8$subpopulation="SP_8"
cellrouter@signatures$SP_9$subpopulation="SP_9"
cellrouter@signatures$SP_10$subpopulation="SP_10"
cellrouter@signatures$SP_11$subpopulation="SP_11"
cellrouter@signatures$SP_12$subpopulation="SP_12"
cellrouter@signatures$SP_13$subpopulation="SP_13"
cellrouter@signatures$SP_14$subpopulation="SP_14"
cellrouter@signatures$SP_15$subpopulation="SP_15"
cellrouter@signatures$SP_16$subpopulation="SP_16"
cellrouter@signatures$SP_17$subpopulation="SP_17"
data=rbind(cellrouter@signatures$SP_1,cellrouter@signatures$SP_2,cellrouter@signatures$SP_3,cellrouter@signatures$SP_4,cellrouter@signatures$SP_5,cellrouter@signatures$SP_6,cellrouter@signatures$SP_7,cellrouter@signatures$SP_8,cellrouter@signatures$SP_9,cellrouter@signatures$SP_10,cellrouter@signatures$SP_11,cellrouter@signatures$SP_12,cellrouter@signatures$SP_13,cellrouter@signatures$SP_14,cellrouter@signatures$SP_15,cellrouter@signatures$SP_16,cellrouter@signatures$SP_17)
data$gene_names=rownames(data)
write.table(data,"results/ydl.dif_genes.txt",sep="\t")
## GRN score for selected transitions
tfs <- find_tfs(species = 'Mm')
save(tfs,file="results/tfs.R")
dev.off()
p <- c('SP_2.SP_1')
pdf("SP_2.SP_1.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_2')
pdf("SP_2.SP_2.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_3')
pdf("SP_2.SP_3.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_4')
pdf("SP_2.SP_4.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_5')
pdf("SP_2.SP_5.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_6')
pdf("SP_2.SP_6.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_7')
pdf("SP_2.SP_7.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_8')
pdf("SP_2.SP_8.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_9')
pdf("SP_2.SP_9.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_10')
pdf("SP_2.SP_10.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_11')
pdf("SP_2.SP_11.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_12')
pdf("SP_2.SP_12.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_13')
pdf("SP_2.SP_13.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_14')
pdf("SP_2.SP_14.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_15')
pdf("SP_2.SP_15.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_16')
pdf("SP_2.SP_16.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_17')
pdf("SP_2.SP_17.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()





savehistory("cellrouter_code.txt")





#只展示前14个--单独
p <- c('SP_2.SP_1')
pdf("SP_2.SP_1-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_2')
pdf("SP_2.SP_2-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_3')
pdf("SP_2.SP_3-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_4')
pdf("SP_2.SP_4-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_5')
pdf("SP_2.SP_5-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_6')
pdf("SP_2.SP_6-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_7')
pdf("SP_2.SP_7-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_8')
pdf("SP_2.SP_8-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_9')
pdf("SP_2.SP_9-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_10')
pdf("SP_2.SP_10-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_11')
pdf("SP_2.SP_11-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_12')
pdf("SP_2.SP_12-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_13')
pdf("SP_2.SP_13_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()

p <- c('SP_2.SP_14')
pdf("SP_2.SP_14_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


p <- c('SP_2.SP_15')
pdf("SP_2.SP_15_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


p <- c('SP_2.SP_16')
pdf("SP_2.SP_16_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


p <- c('SP_2.SP_17')
pdf("SP_2.SP_17_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


##只展示前14个--合并
pdf("SP_2.top.pdf",width =5,height=10)
p <- c('SP_2.SP_1')
#pdf("SP_2.SP_1-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_2')
#pdf("SP_2.SP_2-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_3')
#pdf("SP_2.SP_3-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_4')
#pdf("SP_2.SP_4-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_5')
#pdf("SP_2.SP_5-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_6')
#pdf("SP_2.SP_6-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_7')
#pdf("SP_2.SP_7-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_8')
#pdf("SP_2.SP_8-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_9')
#pdf("SP_2.SP_9-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_10')
#pdf("SP_2.SP_10-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_11')
#pdf("SP_2.SP_11-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_12')
#pdf("SP_2.SP_12-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_13')
#pdf("SP_2.SP_13_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()

p <- c('SP_2.SP_14')
#pdf("SP_2.SP_14_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
p <- c('SP_2.SP_15')
#pdf("SP_2.SP_15_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
p <- c('SP_2.SP_16')
#pdf("SP_2.SP_16_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
p <- c('SP_2.SP_17')
#pdf("SP_2.SP_17_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()




#展示前40个
pdf("SP_2.all.pdf",width =8,height=16)
p <- c('SP_2.SP_1')
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_2')
#pdf("SP_2.SP_2-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_3')
#pdf("SP_2.SP_3-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_4')
#pdf("SP_2.SP_4-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_5')
#pdf("SP_2.SP_5-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_6')
#pdf("SP_2.SP_6-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_7')
#pdf("SP_2.SP_7-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_8')
#pdf("SP_2.SP_8-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_9')
#pdf("SP_2.SP_9-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_10')
#pdf("SP_2.SP_10-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_11')
#pdf("SP_2.SP_11-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_12')
#pdf("SP_2.SP_12-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_13')
#pdf("SP_2.SP_13_3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_14')
#pdf("SP_2.SP_14_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_15')
#pdf("SP_2.SP_15_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_16')
#pdf("SP_2.SP_16_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_17')
#pdf("SP_2.SP_17_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))

dev.off()



transitions <- c('SP_2.SP_13','SP_2.SP_14')
pdf("grnscores_all_13_14.pdf",width =13,height=5)
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=20, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
dev.off()

transitions <- c('SP_2.SP_11','SP_2.SP_12')
pdf("grnscores_all_11_12.pdf",width =13,height=5)
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=20, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
dev.off()



pdf("grnscores_total.pdf",width =13,height=5)
transitions <- c('SP_2.SP_1','SP_2.SP_3')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_4','SP_2.SP_5')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_6','SP_2.SP_7')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_8','SP_2.SP_9')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_10','SP_2.SP_11')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_12','SP_2.SP_13')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_15','SP_2.SP_14')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_16','SP_2.SP_17')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')

dev.off()




##select starting subpopulation,all other subpopulations are targets
sources <- c('SP_15') #from SP_9 to SP_16
targets <- setdiff(as.vector(cellrouter@sampTab$population),sources)
methods <- c("euclidean","maximum","manhattan","canberra","binary",'graph') #graph for distances in KNN
cellrouter <- findpaths(cellrouter,libdir,paste(getwd(),'results/paths',sep='/'),method="graph")
ranks <- c('path_cost','path_flow','rank','length')
cellrouter <- processtrajectories(cellrouter,genes2use,path.rank=ranks[3],num.cells = 3,neighs = 1)
names <- unique(names(cellrouter@pathsinfo$distr))
clusters.show <- names
cellrouter <- correlationpseudotime(cellrouter,type='spearman')
cellrouter <- topgenes(cellrouter,0.85,0.15)
cellrouter <- smoothdynamics(cellrouter,names)
cellrouter <- clusterGenesPseudotime(cellrouter,10)
save(cellrouter,file='results/CellRouter_StemID_Processed.R')
saveRDS(cellrouter,"cellrouter_1.RDS")
cellrouter@signatures$SP_1$subpopulation="SP_1"
cellrouter@signatures$SP_2$subpopulation="SP_2"
cellrouter@signatures$SP_3$subpopulation="SP_3"
cellrouter@signatures$SP_4$subpopulation="SP_4"
cellrouter@signatures$SP_5$subpopulation="SP_5"
cellrouter@signatures$SP_6$subpopulation="SP_6"
cellrouter@signatures$SP_7$subpopulation="SP_7"
cellrouter@signatures$SP_8$subpopulation="SP_8"
cellrouter@signatures$SP_9$subpopulation="SP_9"
cellrouter@signatures$SP_10$subpopulation="SP_10"
cellrouter@signatures$SP_11$subpopulation="SP_11"
cellrouter@signatures$SP_12$subpopulation="SP_12"
cellrouter@signatures$SP_13$subpopulation="SP_13"
cellrouter@signatures$SP_14$subpopulation="SP_14"
cellrouter@signatures$SP_15$subpopulation="SP_15"
cellrouter@signatures$SP_16$subpopulation="SP_16"
cellrouter@signatures$SP_17$subpopulation="SP_17"
data=rbind(cellrouter@signatures$SP_1,cellrouter@signatures$SP_2,cellrouter@signatures$SP_3,cellrouter@signatures$SP_4,cellrouter@signatures$SP_5,cellrouter@signatures$SP_6,cellrouter@signatures$SP_7,cellrouter@signatures$SP_8,cellrouter@signatures$SP_9,cellrouter@signatures$SP_10,cellrouter@signatures$SP_11,cellrouter@signatures$SP_12,cellrouter@signatures$SP_13,cellrouter@signatures$SP_14,cellrouter@signatures$SP_15,cellrouter@signatures$SP_16,cellrouter@signatures$SP_17)
data$gene_names=rownames(data)
write.table(data,"results/ydl.dif_genes.txt",sep="\t")
## GRN score for selected transitions
tfs <- find_tfs(species = 'Mm')
save(tfs,file="results/tfs.R")
dev.off()
p <- c('SP_2.SP_1')
pdf("SP_2.SP_1.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_2')
pdf("SP_2.SP_2.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_3')
pdf("SP_2.SP_3.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_4')
pdf("SP_2.SP_4.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_5')
pdf("SP_2.SP_5.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_6')
pdf("SP_2.SP_6.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_7')
pdf("SP_2.SP_7.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_8')
pdf("SP_2.SP_8.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_9')
pdf("SP_2.SP_9.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_10')
pdf("SP_2.SP_10.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_11')
pdf("SP_2.SP_11.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_12')
pdf("SP_2.SP_12.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_13')
pdf("SP_2.SP_13.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_14')
pdf("SP_2.SP_14.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_15')
pdf("SP_2.SP_15.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_16')
pdf("SP_2.SP_16.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()
p <- c('SP_2.SP_17')
pdf("SP_2.SP_17.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()





savehistory("cellrouter_code.txt")





#只展示前14个--单独
p <- c('SP_2.SP_1')
pdf("SP_2.SP_1-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_2')
pdf("SP_2.SP_2-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_3')
pdf("SP_2.SP_3-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_4')
pdf("SP_2.SP_4-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_5')
pdf("SP_2.SP_5-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_6')
pdf("SP_2.SP_6-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_7')
pdf("SP_2.SP_7-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_8')
pdf("SP_2.SP_8-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_9')
pdf("SP_2.SP_9-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_10')
pdf("SP_2.SP_10-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_11')
pdf("SP_2.SP_11-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_12')
pdf("SP_2.SP_12-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()
p <- c('SP_2.SP_13')
pdf("SP_2.SP_13_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()

p <- c('SP_2.SP_14')
pdf("SP_2.SP_14_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


p <- c('SP_2.SP_15')
pdf("SP_2.SP_15_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


p <- c('SP_2.SP_16')
pdf("SP_2.SP_16_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


p <- c('SP_2.SP_17')
pdf("SP_2.SP_17_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()


##只展示前14个--合并
pdf("SP_2.top.pdf",width =5,height=10)
p <- c('SP_2.SP_1')
#pdf("SP_2.SP_1-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_2')
#pdf("SP_2.SP_2-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_3')
#pdf("SP_2.SP_3-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_4')
#pdf("SP_2.SP_4-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_5')
#pdf("SP_2.SP_5-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_6')
#pdf("SP_2.SP_6-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_7')
#pdf("SP_2.SP_7-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_8')
#pdf("SP_2.SP_8-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_9')
#pdf("SP_2.SP_9-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_10')
#pdf("SP_2.SP_10-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_11')
#pdf("SP_2.SP_11-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_12')
#pdf("SP_2.SP_12-1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_13')
#pdf("SP_2.SP_13_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()

p <- c('SP_2.SP_14')
#pdf("SP_2.SP_14_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
p <- c('SP_2.SP_15')
#pdf("SP_2.SP_15_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
p <- c('SP_2.SP_16')
#pdf("SP_2.SP_16_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
p <- c('SP_2.SP_17')
#pdf("SP_2.SP_17_1.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
dev.off()




#展示前40个
pdf("SP_2.all.pdf",width =8,height=16)
p <- c('SP_2.SP_1')
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_2')
#pdf("SP_2.SP_2-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_3')
#pdf("SP_2.SP_3-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_4')
#pdf("SP_2.SP_4-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_5')
#pdf("SP_2.SP_5-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_6')
#pdf("SP_2.SP_6-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_7')
#pdf("SP_2.SP_7-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_8')
#pdf("SP_2.SP_8-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_9')
#pdf("SP_2.SP_9-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_10')
#pdf("SP_2.SP_10-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_11')
#pdf("SP_2.SP_11-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_12')
#pdf("SP_2.SP_12-3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_13')
#pdf("SP_2.SP_13_3.pdf",width =5,height=10)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_14')
#pdf("SP_2.SP_14_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_15')
#pdf("SP_2.SP_15_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_16')
#pdf("SP_2.SP_16_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))
#dev.off()
p <- c('SP_2.SP_17')
#pdf("SP_2.SP_17_3.pdf",width =6,height=12)
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.7,q.down=0.3,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoin
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
plottr(cellrouter, p, x[[p]]$scores, cluster=FALSE, 2, 4.5, 5.5, paste('results/', p, 'Regulators_transition.pdf',sep=''))

dev.off()



transitions <- c('SP_2.SP_13','SP_2.SP_14')
pdf("grnscores_all_13_14.pdf",width =13,height=5)
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=20, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
dev.off()

transitions <- c('SP_2.SP_11','SP_2.SP_12')
pdf("grnscores_all_11_12.pdf",width =13,height=5)
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=20, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
dev.off()



pdf("grnscores_total.pdf",width =13,height=5)
transitions <- c('SP_2.SP_1','SP_2.SP_3')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_4','SP_2.SP_5')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_6','SP_2.SP_7')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_8','SP_2.SP_9')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_10','SP_2.SP_11')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_12','SP_2.SP_13')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_15','SP_2.SP_14')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')
transitions <- c('SP_2.SP_16','SP_2.SP_17')
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=40, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=40, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_down')

dev.off()



