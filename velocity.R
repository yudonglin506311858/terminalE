

conda activate velocity
export OPENBLAS_NUM_THREADS=1
#导入包
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(velocyto.R)
library(pagoda2)
YDL<- readRDS("SINGLET_terminalE.rds")

head(colnames(YDL))
write.csv(YDL$orig.ident,"data.csv")


#pro
x1 <-velocyto.R::read.loom.matrices(file = "/data/yudonglin/software/cellranger-4.0.0/pro/velocyto/pro.loom", engine = "hdf5r")
splice<-x1$spliced
unsplice<-x1$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste(substring(colnames(splice),5,20),"_1",sep="")
colnames(unsplice) <- paste(substring(colnames(unsplice),5,20),"_1",sep="")
head(colnames(splice))
head(colnames(unsplice))
x1$spliced<-splice
x1$unspliced<-unsplice

#baso1
x2 <-velocyto.R::read.loom.matrices(file = "/data/yudonglin/software/cellranger-4.0.0/baso/velocyto/baso.loom", engine = "hdf5r")
splice<-x2$spliced
unsplice<-x2$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste(substring(colnames(splice),6,21),"_1",sep="")
colnames(unsplice) <- paste(substring(colnames(unsplice),6,21),"_1",sep="")
head(colnames(splice))
head(colnames(unsplice))
x2$spliced<-splice
x2$unspliced<-unsplice



#poly
x4 <-velocyto.R::read.loom.matrices(file = "/data/yudonglin/software/cellranger-4.0.0/poly/velocyto/poly.loom", engine = "hdf5r")
splice<-x4$spliced
unsplice<-x4$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste(substring(colnames(splice),6,21),"_2",sep="")
colnames(unsplice) <- paste(substring(colnames(unsplice),6,21),"_2",sep="")
head(colnames(splice))
head(colnames(unsplice))
x4$spliced<-splice
x4$unspliced<-unsplice

#ortho
x5<-velocyto.R::read.loom.matrices(file = "/data/yudonglin/software/cellranger-4.0.0/ortho/velocyto/ortho.loom", engine = "hdf5r")
splice<-x5$spliced
unsplice<-x5$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste(substring(colnames(splice),7,22),"",sep="")
colnames(unsplice) <- paste(substring(colnames(unsplice),7,22),"",sep="")
head(colnames(splice))
head(colnames(unsplice))
x5$spliced<-splice
x5$unspliced<-unsplice



#reti
x6<-velocyto.R::read.loom.matrices(file = "/data/yudonglin/software/cellranger-4.0.0/reti/velocyto/reti.loom", engine = "hdf5r")
splice<-x6$spliced
unsplice<-x6$unspliced
head(colnames(splice))
head(colnames(unsplice))
colnames(splice) <- paste(substring(colnames(splice),6,21),"-1_2",sep="")
colnames(unsplice) <- paste(substring(colnames(unsplice),6,21),"-1_2",sep="")
head(colnames(splice))
head(colnames(unsplice))
x6$spliced<-splice
x6$unspliced<-unsplice


#整合reti
spliced <- cbind(x1[["spliced"]], x2[["spliced"]],x3[["spliced"]], x4[["spliced"]], x5[["spliced"]], x6[["spliced"]])
unspliced <- cbind(x1[["unspliced"]], x2[["unspliced"]],x3[["unspliced"]], x4[["unspliced"]], x5[["unspliced"]], x6[["unspliced"]])
#未整合reti
spliced <- cbind(x1[["spliced"]], x2[["spliced"]], x4[["spliced"]], x5[["spliced"]])
unspliced <- cbind(x1[["unspliced"]], x2[["unspliced"]], x4[["unspliced"]], x5[["unspliced"]])



emat <- spliced
nmat <- unspliced
emat <- emat[,colSums(is.na(emat))<nrow(emat)]
nmat <- nmat[,colSums(is.na(nmat))<nrow(nmat)]
seurat.object<-YDL
emb <- seurat.object@reductions$umap@cell.embeddings
# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$umap@cell.embeddings)))

# I'm not sure what this parameter does to be honest. 0.02 default
# perform gamma fit on a top/bottom quantiles of expression magnitudes
fit.quantile <- 0.02
# Main velocity estimation

rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,n.cores=48)
# This section gets the colors out of the seurat umap object so that my seurat and velocyto plots use the same color scheme.
library("Seurat")
library("ggplot2")

gg <- UMAPPlot(seurat.object)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)

p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=10,show.grid.flow=T,
                                     min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.5,
                                     do.par=F,cell.border.alpha = 0.1,
                                     main="RNA Velocity")
p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=5,show.grid.flow=T,
                                     min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.5,
                                     do.par=F,cell.border.alpha = 0.1,
                                     main="RNA Velocity")



p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=5,show.grid.flow=T,
                                     min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                     do.par=F,cell.border.alpha = 0.1,
                                     n.cores=48,main="RNA Velocity")
DimPlot(YDL, reduction = 'umap', label=T)
p2<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                                   cell.colors = ac(colors, alpha = 0.5),
                                   cex = 0.8, arrow.scale = 20, show.grid.flow = T,
                                   min.grid.cell.mass = 0.5, grid.n = 40,
                                   arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")
p2<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                                   cell.colors = ac(colors, alpha = 0.5),
                                   cex = 0.8, arrow.scale =10, show.grid.flow = T,
                                   min.grid.cell.mass = 0.5, grid.n = 40,
                                   arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")

dev.off()










emat <- spliced
nmat <- unspliced
seurat.object<-YDL
emb <- seurat.object@reductions$tsne@cell.embeddings
# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$tsne@cell.embeddings)))


# I'm not sure what this parameter does to be honest. 0.02 default
# perform gamma fit on a top/bottom quantiles of expression magnitudes
fit.quantile <- 0.02
# Main velocity estimation

rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,n.cores=48)
# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
library("Seurat")
library("ggplot2")

gg <- TSNEPlot(seurat.object)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)
p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=2,show.grid.flow=T,
                                     min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                     do.par=F,cell.border.alpha = 0.1,
                                     n.cores=48,main="RNA Velocity")
DimPlot(YDL, reduction = 'tsne', label=T)
p2<-show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = 'sqrt',
                                   cell.colors = ac(colors, alpha = 0.5),
                                   cex = 0.8, arrow.scale = 20, show.grid.flow = T,
                                   min.grid.cell.mass = 0.5, grid.n = 40,
                                   arrow.lwd = 1,do.par = T, cell.border.alpha = 0.1,n.cores=48,main="RNA Velocity")

dev.off()

