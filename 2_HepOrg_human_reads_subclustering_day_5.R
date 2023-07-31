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
setwd("~/HepOrgMalaria/humanreads/Day5Subclustering/")

##Initial clustering done as described in "1_HepOrg_human_reads_initial.R"


## Split dataset into separate datasets per batch
cleaned.list <- SplitObject(cleaned, split.by = "day")

dim(cleaned.list$day5)

## Batch Day 5
cleaned.list$day5 <- SCTransform(cleaned.list$day5, vars.to.regress = c("nCount_RNA","nFeature_RNA","plate"), do.scale = T, verbose = T)

# These are now standard steps in the Seurat workflow for visualization and clustering
cleaned.list$day5  <- RunPCA(object = cleaned.list$day5 , verbose = T)
cleaned.list$day5  <- RunTSNE(object = cleaned.list$day5 , verbose = T)
cleaned.list$day5  <- RunUMAP(object = cleaned.list$day5 , dims = 1:10, verbose = T)
cleaned.list$day5  <- FindNeighbors(object = cleaned.list$day5 , dims = 1:10, verbose = T)
cleaned.list$day5 <- FindClusters(cleaned.list$day5, verbose = T, resolution = 0.4, algorithm = 1)

#Safe clusters
Idents(cleaned.list$day5)
cleaned.list$day5$initialclusters <- Idents(cleaned.list$day5)
Idents(cleaned.list$day5)

##Define colours
infectionstatuscolors <- c("grey","black")

pdf("Figure_3A_infectionstatus_tsne.pdf",)
DimPlot(cleaned.list$day5, group.by = c("infectionstatus"), reduction = "tsne", cols = c("grey","black"), pt.size = 1)
dev.off()

##Differentially expressed (DE) genes comparing infected and non-infected HepOrg cells at day 5
DE.genes <- FindAllMarkers(object = cleaned.list$day5, only.pos = TRUE, min.pct = 0.25, 
                                  thresh.use = 0.25)
write.csv(DE.genes,"Supplementary_Data_2.csv")


##TOP DE genes
DE.genes %>%
        group_by(cluster) %>%
        top_n(n = 50, wt = avg_log2FC) -> top50

pdf("Figure_3B_TOP_DE_genes_heatmap.pdf")
colorlist <- list(infectionstatus=infectionstatuscolors)
names(colorlist[["infectionstatus"]]) <- c("no","yes")
DoMultiBarHeatmap(cleaned.list$day5, features = top50$gene, group.by="infectionstatus", cols.use=colorlist) + scale_fill_gradientn(colors = viridis(10))
dev.off()


pdf("Figure_3C_selected_DE_genes.pdf")
my_comparison1 <- list( c("no", "yes"))
cpt1a <- VlnPlot(cleaned.list$day5, log = T, features = "CPT1A", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
fasn <- VlnPlot(cleaned.list$day5, log = T, features = "FASN", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
apoa1 <- VlnPlot(cleaned.list$day5, log = T, features = "APOA1", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
hmgcr <- VlnPlot(cleaned.list$day5, log = T, features = "HMGCR", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
g6pc <- VlnPlot(cleaned.list$day5, log = T, features = "G6PC", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 

pcsk9 <- VlnPlot(cleaned.list$day5, log = T, features = "PCSK9", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
apob <- VlnPlot(cleaned.list$day5, log = T, features = "APOB", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
ppara <- VlnPlot(cleaned.list$day5, log = T, features = "PPARA", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
lss <- VlnPlot(cleaned.list$day5, log = T, features = "LSS", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
mttp <- VlnPlot(cleaned.list$day5, log = T, features = "MTTP", cols = infectionstatuscolors, group.by = "infectionstatus", pt.size = 0) + geom_boxplot(width=0.2,fill="white") + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparison1) + theme(legend.position = "none") 
CombinePlots(list(cpt1a,fasn,apoa1,hmgcr,g6pc,pcsk9,apob,ppara,lss,mttp),ncol = 5)
dev.off()


##Genes upregulated in infected HepOrg cells
infectionmarkers <- FindMarkers(cleaned.list$day5, 
                                ident.1 = "yes", 
                                ident.2 = "no", 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                thresh.use = 0.25
                               )

##Enrichr analysis
dbs <- listEnrichrDbs()
dbs <- "WikiPathway_2021_Human"

infectionmarker.genes <- rownames(infectionmarkers)
enriched.infected <- enrichr(infectionmarker.genes, dbs)
enriched.infected[["WikiPathway_2021_Human"]]
write.csv(enriched.infected[["WikiPathway_2021_Human"]],"Supplementary_Data_3.csv")


