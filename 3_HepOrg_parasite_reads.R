##Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(DESeq2)
library(RColorBrewer)
library(enrichR)
library(clustree)
library(DoMultiBarHeatmap)
library(rlang)
library(magrittr)


##Set working directory
setwd("~/HepOrgMalaria/parasitereads/")


##Load in rawdata and metadata
rawdata <- read.csv("HepOrgMalaria_parasitereads_rawdata.csv",sep=",")
names   <- make.unique(rawdata[,1])
rownames(rawdata) <- names
rawdata <- rawdata[,-1]

rawdata[is.na(rawdata)] <- 0

rawdata <- rawdata[grep("ERCC-",rownames(rawdata), invert=TRUE),]
rawdata <- rawdata[grep("MT-",rownames(rawdata), invert=TRUE),]


metadata <- read.csv("HepOrgMalaria_parasitereads_metadata.csv", sep=",")
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]


#Creat Seurat object
pbmc <- CreateSeuratObject(counts = rawdata, 
                           meta.data = metadata
                          )

pbmc <- subset(x = pbmc, 
               subset = nCount_RNA > 50 & nCount_RNA < Inf
              )

##SCTransform
pbmc <- SCTransform(pbmc, vars.to.regress = c("nCount_RNA","nFeature_RNA","plate"), do.scale = T, verbose = T)

# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc  <- RunPCA(object = pbmc , verbose = T)
pbmc  <- RunTSNE(object = pbmc , verbose = T, check_duplicates = F)
pbmc  <- RunUMAP(object = pbmc , dims = 1:10, verbose = T)
pbmc  <- FindNeighbors(object = pbmc , dims = 1:10, verbose = T)
pbmc  <- FindClusters(object = pbmc , verbose = T)
pbmc  <- FindClusters(object = pbmc , resolution = 0.4, verbose = T)

# Hexadecimal color specification 
cl.cols <- 6
clustercols <- colorRampPalette(brewer.pal(6, "Set2"))(cl.cols)

