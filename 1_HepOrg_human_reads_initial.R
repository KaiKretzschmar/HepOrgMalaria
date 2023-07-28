##Load packages
library(Seurat)
library(ggplot2)
library(sctransform)
library(DESeq2)
library(RColorBrewer)
library(enrichR)
library(clustree)


##Set working directory
setwd("~/HepOrgMalaria/humanreads/InitialAnalysis/")


##Load in rawdata and metadata
rawdata <- read.csv("HepOrgMalaria_humanreads_rawdata.csv",sep=",")
names   <- make.unique(rawdata[,1])
rownames(rawdata) <- names
rawdata <- rawdata[,-1]

rawdata[is.na(rawdata)] <- 0

rawdata <- rawdata[grep("ERCC-",rownames(rawdata), invert=TRUE),]
rawdata <- rawdata[grep("MT-",rownames(rawdata), invert=TRUE),]


metadata <- read.csv("HepOrgMalaria_humanreads_metadata.csv", sep=",")
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]


##Create Seurat object
initial <- CreateSeuratObject(counts = rawdata, 
                              meta.data = metadata
                             )

initial <- subset(initial, 
                  subset = nCount_RNA > 2000 & nCount_RNA < Inf
                 )


##SCTransform
initial  <- SCTransform(initial, 
                        vars.to.regress = c("nCount_RNA","nFeature_RNA","plate"), 
                        do.scale = T, 
                        verbose = T
                       )

##Analysis of necrosis
#See van den Brink et al. (2017), Nature Methods - PMID: 28960196
Necrosisgenes <- c("FOSB","FOS","JUN","JUNB","ATF3","EGR1","HSPA1A","HSPA1B","HSPB1","IER3","IER2","DUSP1")

#Computes an enrichment score for necrosis genes
initial  <- AddModuleScore(
  object = initial,
  features = Necrosisgenes,
  ctrl = 100,
  name = "Necrosis_Scoring"
)

#Filtering out of cells with high necrosis gene score
cleaned  <- subset(initial, subset = Necrosis_Scoring1 < 1)


##SCTransform, dimensional reduction and clustering of cleaned data
cleaned  <- SCTransform(cleaned, 
                        vars.to.regress = c("nCount_RNA","nFeature_RNA","plate"), 
                        do.scale = T, 
                        verbose = T
                       )

#These are now standard steps in the Seurat workflow for visualization and clustering
cleaned  <- RunPCA(object = cleaned , verbose = T)
cleaned  <- RunTSNE(object = cleaned , verbose = T)
cleaned  <- RunUMAP(object = cleaned , dims = 1:10, verbose = T)
cleaned  <- FindNeighbors(object = cleaned , dims = 1:10, verbose = T)
cleaned  <- FindClusters(cleaned, verbose = T, resolution = 0.4, algorithm = 1)

#Safe Seurat clusters
cleaned$initialclusters <- Idents(cleaned)

#Define colours
cl.cols <- 5
clustercols <- colorRampPalette(brewer.pal(5, "Set1"))(cl.cols)
conditioncolors <- c("#5ab4ac","#d8b365")
daycolors <- c("#1f78b4","#33a02c")
gatingcolors <- c("#7fbf7b","#af8dc3")
infectionstatuscolors <- c("grey","red")
trafficlightinfectioncolors <- c("#A6D854","#FFD92F","#E41A1C")
cellcyclecolors <- c("#b3cde3","#ccebc5","#fbb4ae")


#Make cluster plots 
pdf("Figure_2A_tsne.pdf")
DimPlot(cleaned, 
        reduction = "tsne", 
        pt.size = 1, 
        cols = clustercols
       ) 
+ ggtitle('Cell Clusters') 
dev.off()

pdf("Figure_S3B_umap.pdf")
DimPlot(cleaned, 
        reduction = "umap", 
        pt.size = 1, 
        cols = clustercols
       ) 
+ ggtitle('Cell Clusters') 
dev.off()

pdf("Figure_S3C_pca.pdf")
DimPlot(cleaned, 
        reduction = "pca", 
        pt.size = 1, 
        cols = clustercols
       ) 
+ ggtitle('Cell Clusters') 
dev.off()


##Make SNN graph - source: https://romanhaa.github.io/projects/scrnaseq_workflow/#snn-graph
SCT_snn <- cleaned@graphs$SCT_snn %>%
  as.matrix() %>%
  ggnetwork() %>%
  left_join(cleaned@meta.data %>% mutate(vertex.names = rownames(.)), by = 'vertex.names')

pdf("Figure_S3A_SNN_graph.pdf")
ggplot(SCT_snn, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = 'grey50', alpha = 0.05) +
  geom_nodes(aes(color = initialclusters), size = 1) +
  scale_color_manual(
    name = 'Cluster', values = clustercols,
    guide = guide_legend(ncol = 1, override.aes = list(size = 2))
  ) +
  theme_blank() +
  theme(legend.position = 'none') +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(cleaned@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5)
dev.off()


##Make cluster tree
pdf("Figure_S3A_cluster_tree.pdf")
cleaned <- BuildClusterTree(
  cleaned ,
  dims = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE
)


##Cell cycle analysis
#Read in list of cell cycle markers from Tirosh et al. (2016), Nature - PMID: 27806376
cc.genes <- readLines(con = "CellCycleGenes.txt")

#Cell cycle scoring
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
cleaned <- CellCycleScoring(cleaned, 
                            assay = 'SCT',
                            s.features = s.genes, 
                            g2m.features = g2m.genes
                           )

#Make tsne plot displaying cell cycle phases
pdf("Figure_S3F_cell_cycle.pdf")
DimPlot(cleaned, 
        group.by = c("Phase"), 
        reduction = "tsne", 
        cols = cellcyclecolors, 
        pt.size = 1
       ) 
+ ggtitle('Cell Cycle Phases') 
dev.off()

##Define hepatocyte genes
hepatocytemarker = c("ALB","AFP","RBP4","FABP1","SERPINA1","ASGR2","ASGR1","APOA2","APOC3")

#Make tsne plots for hepatocyte marker gene expression
pdf("Figure_S4A_hepatocyte_marker_expression.pdf")
FeaturePlot(cleaned, 
            features = hepatocytemarker, 
            cols = viridis(10), 
            reduction = "tsne", 
            pt.size = 1, 
            ncol = 3
           )
dev.off()

#Enrichment scoring for hepatocyte marker genes
cleaned <- AddModuleScore(
  object = cleaned,
  features = hepatocytemarker,
  ctrl = 100,
  name = "Hepatocyte_Scoring"
           )

pdf("Figure_S4B_hepatocyte_marker_scoring.pdf")
FeaturePlot(cleaned, 
            features = "Hepatocyte_Scoring1", 
            cols = inferno(10), 
            reduction = "tsne", 
            pt.size = 1
           )
dev.off()


##Define cholangiocyte genes
cholangiocytemarker = c("KRT19","KRT8","KRT18","EPCAM","KRT7")

#Make tsne plots for cholangiocyte marker gene expression
pdf("Figure_S5A_cholangiocyte_marker_expression.pdf")
FeaturePlot(cleaned, 
            features = cholangiocytemarker, 
            cols = viridis(10), 
            reduction = "tsne", 
            pt.size = 1, 
            ncol = 3
           )
dev.off()

#Enrichment scoring for cholangiocyte marker genes
cleaned <- AddModuleScore(
  object = cleaned,
  features = cholangiocytemarker,
  ctrl = 100,
  name = "Cholangiocyte_Scoring"
)

pdf("Figure_S5B_cholangiocyte_marker_scoring.pdf")
FeaturePlot(cleaned, features = "Cholangiocyte_Scoring1", cols = inferno(10), reduction = "tsne", pt.size = 1)
dev.off()


##Differentially expressed (DE) genes per cluster - only positive markers are reported per Seurat cluster
cleaned.markers <- FindAllMarkers(object = cleaned, only.pos = TRUE, min.pct = 0.25, 
                                       thresh.use = 0.25)
write.csv(cleaned.markers,"Supplemental_Data_1_Cluster_markers_Wilcox.csv")


##Make violin plots for lineage marker gene expression per Seurat cluster
pdf("Figure_2B_lineage_marker.pdf", width = 3, height = 3)
VlnPlot(obj = cleaned, 
        features = c("ALB","MKI67","KRT7","EPCAM"), 
        cols = clustercols, 
        fill.by = "ident", 
        flip = TRUE, 
        stack = TRUE
       ) 
+ NoLegend()
dev.off()



