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

##SCTransform, dimensional reduction and clustering of cleaned data
pbmc <- SCTransform(pbmc, vars.to.regress = c("nCount_RNA","nFeature_RNA","plate"), do.scale = T, verbose = T)
pbmc  <- RunPCA(object = pbmc , verbose = T)
pbmc  <- RunTSNE(object = pbmc , verbose = T, check_duplicates = F)
pbmc  <- RunUMAP(object = pbmc , dims = 1:10, verbose = T)
pbmc  <- FindNeighbors(object = pbmc , dims = 1:10, verbose = T)
pbmc  <- FindClusters(object = pbmc , verbose = T)
pbmc  <- FindClusters(object = pbmc , resolution = 0.4, verbose = T)

#Define color codes
cl.cols <- 6
clustercols <- colorRampPalette(brewer.pal(6, "Set2"))(cl.cols)

#Make cluster plots
pdf("Figure_4A_initial_clustering_tsne.pdf")
DimPlot(pbmc, 
        reduction = "tsne", 
        pt.size = 3, 
        cols = clustercols
       )
dev.off()

##Assessing different experimenal categories
pdf("Figure_4B_stages_tsne.pdf", width = 5, height = 4.5)
DimPlot(pbmc, 
        group.by = c("stage"), 
        pt.size = 3, 
        reduction = "tsne", 
        cols = c("#5ab4ac","#d8b365")
       )
dev.off()

#Safe initial cluster IDs
Idents(pbmc)
pbmc$initialclusters <- Idents(pbmc)
Idents(pbmc)

Idents(pbmc) <- pbmc@meta.data$stage
Idents(pbmc)

##Differentially expressed (DE) genes comparing liver stage and blood stage
stage.markers <- FindAllMarkers(object = pbmc, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                               thresh.use = 0.25
                                   )
write.csv(stage.markers,"Supplementary_Data_4.csv")

#Plot heatmap for stage marker genes
stagecols <- c("#5ab4ac","#d8b365")
colorlist <- list(stage=stagecols)

stage.markers %>%
        group_by(stage) %>%
        top_n(n = 50, wt = avg_log2FC) -> top50

pdf("Figure_4C_heatmap_stage_marker_genes.pdf", width = 10, height = 5)
DoMultiBarHeatmap(pbmc, 
                  features = top50$gene, 
                  group.by="stage", 
                  cols.use=colorlist
                 ) 
+ scale_fill_gradientn(colors = viridis(10))

#Plot selected DE genes up in liver stage
pdf("Figure_4D_liver_stage_genes.pdf")
my_comparisonp <- list( c("liver","blood"))
csp <- VlnPlot(pbmc, log = T, features = "XM001351086.1", cols = c("#5ab4ac","#d8b365"), group.by = "stage", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisonp) + theme(legend.position = "none") + ggtitle('CSP (PF3D7_0304600)') 
lisp <- VlnPlot(pbmc, log = T, features = "XM001348316.1", cols = c("#5ab4ac","#d8b365"), group.by = "stage", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisonp) + theme(legend.position = "none") + ggtitle('LISP1 (PF3D7_1418100)') 
slarp <- VlnPlot(pbmc, log = T, features = "XM001348111.1", cols = c("#5ab4ac","#d8b365"), group.by = "stage", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisonp) + theme(legend.position = "none") + ggtitle('SLARP (PF3D7_1147000)') 
CombinePlots(list(csp,lisp,slarp),ncol = 3)
dev.off()


