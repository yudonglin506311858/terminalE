#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
options(encoding = "UTF-8")
Sys.setenv(LANGUAGE = "en") #
gc()
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(shiny)
library('stringr') 
library(rsconnect)
#library(shinydashboard)
#options(rsconnect.max.bundle.size=5145728000)

#seurat_obj<-readRDS("D:/analysis/forpublication/bulk_sc/terminalE.rds")
seurat_obj<-readRDS("terminalE.rds")
merge_tf<-read.csv("E0_data.csv")
label<-merge_tf$label
# Seurat可视化SCENIC结果
# 
# 把SCENIC结果中最重要的regulonAUC矩阵导入Seurat，这样得到的可视化结果更容易与我们之前的分析联系起来。

##导入原始regulonAUC矩阵
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

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(seurat_obj, BINmatrix)
scRNAbin@assays$integrated <- NULL
#saveRDS(scRNAbin, 'scRNAbin.rds')


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("single cell landscape of mouse terminal erythropoiesis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    
    
    sidebarPanel(
      
      
      conditionalPanel(
        condition = "input.smoother == ture",
        selectInput("comment","Reduction",list("umap","tsne","pca"))
      ),
      
      
      
      conditionalPanel(
        condition = "input == ture",
        selectInput("gene",'Gene',as.list(c(rownames(seurat_obj)) ))
        
      ),
      
      
      
      conditionalPanel(
        condition = "input == ture",
        selectInput("Regulon",'Regulon',as.list(c(label) ))
        
      )
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      #h3("DimPlot"),
      
      tabsetPanel(
        tableOutput('contents'),
        tabPanel("Reduction",plotOutput("distPlot")),
        tabPanel("VlnPlot",plotOutput("VlnPlot")),
        tabPanel("RidgePlot",plotOutput("RidgePlot")),
        tabPanel("FeaturePlot",plotOutput("FeaturePlot")),
        tabPanel("DotPlot",plotOutput("DotPlot")),
        tabPanel("DoHeatmap",plotOutput("DoHeatmap")),
        tabPanel("Regulon",plotOutput("Regulon")),
        tabPanel('marker',dataTableOutput('marker'))
        
      )
      
    )
  )
  
))


#reference:https://www.jianshu.com/p/de33438a72a2?utm_campaign=haruki
