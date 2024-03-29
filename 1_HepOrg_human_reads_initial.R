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
infectionstatuscolors <- c("grey","black")
parasitetranscriptscolors <- c("#A6D854","#FFD92F","#E41A1C")
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
dev.off()

##Assessing different experimenal categories
pdf("Figure_S3D_left_experimental_conditions_tsne.pdf")
DimPlot(cleaned, 
        group.by = "condition", 
        reduction = "tsne", 
        cols = conditioncolors, 
        pt.size = 2
       ) 
dev.off()

pdf("Figure_S3D_centre_clusters_split_by_experimental_conditions_tsne.pdf")
DimPlot(cleaned, 
        group.by = "initialclusters", 
        split.by = "condition", 
        reduction = "tsne", 
        ncol = 2, 
        cols = clustercols, 
        pt.size = 2
       ) 
dev.off()

ggplot(cleaned@meta.data, 
       aes(x=condition, 
           fill=initialclusters)
      ) 
+ geom_bar(position = "fill") 
+ theme_classic() 
+ scale_fill_manual(values=clustercols)
ggsave("Figure_S3D_right_condition_vs_clusters.pdf")


pdf("Figure_S3E_collection_day_tsne.pdf")
DimPlot(cleaned, 
        group.by = "day", 
        reduction = "tsne", 
        cols = daycolors, 
        pt.size = 2
       ) 
dev.off()

pdf("Figure_S3G_left_gating_tsne.pdf")
DimPlot(cleaned, 
        group.by = "gating", 
        reduction = "tsne", 
        cols = gatingcolors
        pt.size = 2
       )
dev.off()

pdf("Figure_S3G_centre_infectionstatus_tsne.pdf")
DimPlot(cleaned, 
        group.by = "infectionstatus", 
        reduction = "tsne", 
        cols = infectionstatuscolors, 
        pt.size = 2
       ) 
dev.off()

ggplot(condition$sample@meta.data, 
       aes(x=gating, 
           fill=infectionstatus)
      ) 
+ geom_bar(position="fill") 
+ theme_classic() 
+ scale_fill_manual(values=infectionstatuscolors) 
ggsave("Figure_S3G_right_gating-vs-infectionstatus.pdf")


##Differentially expressed (DE) genes per cluster - only positive markers are reported per Seurat cluster
cleaned.markers <- FindAllMarkers(object = cleaned, only.pos = TRUE, min.pct = 0.25, 
                                       thresh.use = 0.25)
write.csv(cleaned.markers,"Supplementary_Data_1_Cluster_markers_Wilcox.csv")

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

pdf("Figure_S4B_hepatocyte_marker_scoring_tsne.pdf")
FeaturePlot(cleaned, 
            features = "Hepatocyte_Scoring1", 
            cols = inferno(10), 
            reduction = "tsne", 
            pt.size = 1
           )
dev.off()

pdf("Figure_2C_hepatocyte_marker_scoring_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "Hepatocyte_Scoring1", 
        cols = clustercols, 
        group.by = "initialclusters", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, 
               fill="white"
              ) 
+ NoLegend() 
+ labs(title = "Hepatocyte Marker", 
       x = "Clusters", 
       y = "Score"
      ) 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     ref.group = "3"
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

pdf("Figure_S5B_cholangiocyte_marker_scoring_tsne.pdf")
FeaturePlot(cleaned, features = "Cholangiocyte_Scoring1", cols = inferno(10), reduction = "tsne", pt.size = 1)
dev.off()


pdf("Figure_S5C_cholangiocyte_marker_scoring_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "Cholangiocyte_Scoring1", 
        cols = clustercols, 
        group.by = "initialclusters", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, 
               fill="white"
              ) 
+ NoLegend() 
+ labs(title = "Cholangiocyte Marker", 
       x = "Clusters", 
       y = "Score"
      ) 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     ref.group = "4"
                    )
dev.off()


##Infection analysis
pdf("Figure_2D_parasite_transcripts_tsne.pdf")
DimPlot(cleaned, group.by = "parasitetranscripts", reduction = "tsne", cols = parasitetranscriptscolors, pt.size = 1)
dev.off()

ggplot(cleaned@meta.data, aes(x=parasitetranscripts, fill=initialclusters)) 
+ geom_bar(position = "fill") 
+ theme_classic() 
+ scale_fill_manual(values =  clustercols)
ggsave("Figure_2E_parasite_transcripts-vs-clusters.pdf")

ggplot(cleaned@meta.data, aes(x=initialclusters, fill=parasitetranscripts)) 
+ geom_bar(position = "fill") 
+ theme_classic() 
+ scale_fill_manual(values=parasitetranscriptscolors)
ggsave("Figure_2F_clusters-vs-parasite_transcripts.pdf")

ggplot(cleaned@meta.data, aes(x=day, fill=parasitetranscripts)) 
+ geom_bar(position= "fill") 
+ theme_classic() 
+ scale_fill_manual(values=parasitetranscriptscolors) 
ggsave("Figure_2G_days-vs-parasite_transcripts.pdf")


##Analysis of entry factors expression
#SCARB1
pdf("Figure_2H_SCARB1_tsne.pdf")
FeaturePlot(cleaned, 
            features = "SCARB1", 
            cols = viridis(10), 
            reduction = "tsne", 
            pt.size = 2
           )
dev.off()

pdf("Figure_2I_SCARB1_cluster_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "SCARB1", 
        cols = clustercols, 
        group.by = "initialclusters", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     ref.group = "3"
                    ) 
dev.off()


statsparasitetranscripts <- list( c("none", "low"),c("none", "high"),c("high", "low"))

pdf("Figure_2J_SCARB1_parasite_transcripts_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "SCARB1", 
        cols = c("#A6D854","#FFD92F","#FF4F51"), 
        group.by = "parasitetranscripts", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     comparisons = statsparasitetranscripts
                    ) 
dev.off()

pdf("Figure_S6D_SCARB1_split_by_condition_tsne.pdf")
FeaturePlot(cleaned, 
            features = "SCARB1", 
            group.by = "condition", 
            cols = viridis(10), 
            reduction = "tsne", 
            pt.size = 2
           ) 
dev.off()

pdf("Figure_S6D_SCARB1_vs_condition_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "SCARB1", 
        cols = c("#5ab4ac","#d8b365"), 
        group.by = "condition", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2,fill="white") 
+ stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "ctrl")
dev.off()

#CD81
pdf("Figure_S6A_left_CD81_tsne.pdf")
FeaturePlot(cleaned, 
            features = "CD81", 
            cols = viridis(10), 
            reduction = "tsne", 
            pt.size = 2
           )
dev.off()

pdf("Figure_S6A_right_CD81_cluster_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "CD81", 
        cols = clustercols, 
        group.by = "initialclusters", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     ref.group = "4"
                    ) 
dev.off()


statsparasitetranscripts <- list( c("none", "low"),c("none", "high"),c("high", "low"))

pdf("Figure_S6A_centre_CD81_parasite_transcripts_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "CD81", 
        cols = c("#A6D854","#FFD92F","#FF4F51"), 
        group.by = "parasitetranscripts", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     comparisons = statsparasitetranscripts
                    ) 
dev.off()


#EPHA2
pdf("Figure_S6B_left_EPHA2_tsne.pdf")
FeaturePlot(cleaned, 
            features = "EPHA2", 
            cols = viridis(10), 
            reduction = "tsne", 
            pt.size = 2
           )
dev.off()

pdf("Figure_S6B_right_CD81_cluster_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "EPHA2", 
        cols = clustercols, 
        group.by = "initialclusters", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     ref.group = "4"
                    ) 
dev.off()


statsparasitetranscripts <- list( c("none", "low"),c("none", "high"),c("high", "low"))

pdf("Figure_S6B_centre_EPHA2_parasite_transcripts_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "EPHA2", 
        cols = c("#A6D854","#FFD92F","#FF4F51"), 
        group.by = "parasitetranscripts", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     comparisons = statsparasitetranscripts
                    ) 
dev.off()

#Correlation plots for SCARB1 expression vs. hepatocyte, cholangiocyte and progenitor markers in sample plates only

condition <- SplitObject(cleaned, split.by = "condition")
pdf("Figure S6C_SCARB1_vs_marker_genes_in_sample_plates_only.pdf", width = 8, height = 8)
alb <- FeatureScatter(condition$sample, 
                      feature1 = "ALB", 
                      feature2 = "SCARB1", 
                      group.by = "parasitetranscripts", 
                      cols = c("#A6D854","#FFD92F","#FF4F51"), 
                      pt.size = 2, 
                      slot = "counts"
                     ) 
+ geom_smooth(method="lm", color = "#000000")
afb <- FeatureScatter(condition$sample, 
                      feature1 = "AFP", 
                      feature2 = "SCARB1", 
                      group.by = "parasitetranscripts", 
                      cols = c("#A6D854","#FFD92F","#FF4F51"), 
                      pt.size = 2, 
                      slot = "counts"
                     )  
+ geom_smooth(method="lm", color = "#000000")
rbp4 <- FeatureScatter(condition$sample, 
                       feature1 = "RBP4", 
                       feature2 = "SCARB1", 
                       group.by = "parasitetranscripts", 
                       cols = c("#A6D854","#FFD92F","#FF4F51"), 
                       pt.size = 2, 
                       slot = "counts"
                      ) 
+ geom_smooth(method="lm", color = "#000000")
epcam <- FeatureScatter(condition$sample, 
                        feature1 = "EPCAM", 
                        feature2 = "SCARB1", 
                        group.by = "parasitetranscripts", 
                        cols = c("#A6D854","#FFD92F","#FF4F51"), 
                        pt.size = 2, 
                        slot = "counts"
                       )  
+ geom_smooth(method="lm", color = "#000000")
krt19 <- FeatureScatter(condition$sample, 
                        feature1 = "KRT19", 
                        feature2 = "SCARB1", 
                        group.by = "parasitetranscripts", 
                        cols = c("#A6D854","#FFD92F","#FF4F51"), 
                        pt.size = 2, 
                        slot = "counts"
                       )  
+ geom_smooth(method="lm", color = "#000000")
krt8 <- FeatureScatter(condition$sample, 
                       feature1 = "KRT8", 
                       feature2 = "SCARB1", 
                       group.by = "parasitetranscripts", 
                       cols = c("#A6D854","#FFD92F","#FF4F51"), 
                       pt.size = 2, 
                       slot = "counts"
                      ) 
+ geom_smooth(method="lm", color = "#000000")
mki67 <- FeatureScatter(condition$sample, 
                        feature1 = "MKI67", 
                        feature2 = "SCARB1", 
                        group.by = "parasitetranscripts", 
                        cols = c("#A6D854","#FFD92F","#FF4F51"), 
                        pt.size = 2, 
                        slot = "counts"
                       )  
+ geom_smooth(method="lm", color = "#000000")
CombinePlots(list(alb,afb,rbp4,epcam,krt19,krt8,mki67), ncol = 3)
dev.off()


##GLUL-GS expression
pdf("Figure_S11_left_GLUL_tsne.pdf")
FeaturePlot(cleaned, 
            features = "GLUL", 
            cols = viridis(10), 
            reduction = "tsne", 
            pt.size = 2
           )
dev.off()

pdf("Figure_S11_right_GLUL_cluster_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "GLUL", 
        cols = clustercols, 
        group.by = "initialclusters", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     ref.group = "3"
                    ) 
dev.off()


statsparasitetranscripts <- list( c("none", "low"),c("none", "high"),c("high", "low"))
pdf("Figure_S11_centre_GLUL_parasite_transcripts_violin_plot.pdf")
VlnPlot(cleaned, 
        log = T, 
        features = "GLUL", 
        cols = c("#A6D854","#FFD92F","#FF4F51"), 
        group.by = "parasitetranscripts", 
        pt.size = 0
       ) 
+ geom_boxplot(width=0.2, fill="white") 
+ stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     comparisons = statsparasitetranscripts
                    ) 
dev.off()






