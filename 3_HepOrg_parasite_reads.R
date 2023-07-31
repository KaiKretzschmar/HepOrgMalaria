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

