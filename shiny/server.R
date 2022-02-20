#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
options(encoding = "UTF-8")
library(rsconnect)
Sys.setenv(LANGUAGE = "en") #
gc()
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(shiny)
library('stringr') 
library(ggplot2) 
library(plyr) 
library(rpart) 
library(RColorBrewer) 
library(reshape2)
#options(rsconnect.max.bundle.size=5145728000)
options(shiny.maxRequestSize=7000*1024^2)
library(shinydashboard)
#seurat_obj<-readRDS("D:/analysis/forpublication/bulk_sc/terminalE.rds")
seurat_obj<-readRDS("terminalE.rds")

# seurat_obj<-subset(seurat_obj, downsample = 500)
# saveRDS(seurat_obj,"terminalE_less.rds")
# YDL.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# # install.packages("magrittr") # package installations are only needed the first time you use it
# # install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
# YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# #marker
# 
# write.csv(YDL.markers,file="allmarker_YDL.csv")
#YDL.markers<-read.csv("D:/analysis/forpublication/bulk_sc/allmarker_YDL.csv",row.names = 1)
YDL.markers<-read.csv("allmarker_YDL.csv",row.names = 1)
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_logFC)
# Define server logic required to draw a histogram


merge_tf<-read.csv("E0_data.csv")
label<-merge_tf$label

AUCmatrix <- readRDS("3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(seurat_obj, AUCmatrix)
scRNAauc@assays$integrated <- NULL
#saveRDS(scRNAauc,'scRNAauc.rds')

BINmatrix <- readRDS("4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(seurat_obj, BINmatrix)
scRNAbin@assays$integrated <- NULL
#saveRDS(scRNAbin, 'scRNAbin.rds')

shinyServer(function(input, output) {
  output$contents <-  DT::renderDataTable({
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)
    }else{
      seurat_obj<- readRDS(inFile$datapath)
      seurat_obj
      req(seurat_obj)
      
    }
  })
  
  
  output$distPlot <- renderPlot({
    DimPlot(seurat_obj, label = TRUE, pt.size = 1.5,reduction   = input$comment) 
    
  })
  
  
  output$VlnPlot<-renderPlot({
    VlnPlot(seurat_obj, features = input$gene, pt.size = 0, ncol = 1)
  })
  
  
  output$RidgePlot<-renderPlot({
    RidgePlot(seurat_obj, features = input$gene)
  })
  
  
  output$FeaturePlot<-renderPlot({
    FeaturePlot(seurat_obj, features = input$gene, reduction= input$comment, pt.size = 1.5, ncol = 1,cols = c("gray", "red"))
  })
  
  
  output$DotPlot<-renderPlot({
    DotPlot(seurat_obj, features = input$gene)#+ RotatedAxis()
  })
  
  
  output$DoHeatmap<-renderPlot({
    DoHeatmap(seurat_obj, features = top10$gene)+ NoLegend()
  })
  
  output$Regulon<-renderPlot({
    p1 = FeaturePlot(scRNAbin, features=input$Regulon, label=T, reduction = input$comment,cols = c("gray", "red"))
    p2 = FeaturePlot(scRNAauc, features=input$Regulon, label=T, reduction =input$comment,cols = c("gray", "red"))
    feature<-str_sub(input$Regulon,1,str_locate(input$Regulon,"_")[1]-1)
    p3<-FeaturePlot(object = seurat_obj, reduction =input$comment,features = feature,cols = c("gray", "red"))#actin
    plotc = p1|p2|p3
    print(plotc)
  })
  
  
  output$marker<-renderDataTable(
    
    df <-data.frame( YDL.markers),
    options = list(
      pageLength = 10
    ))
  
  
  
})




#reference:https://www.jianshu.com/p/de33438a72a2?utm_campaign=haruki
