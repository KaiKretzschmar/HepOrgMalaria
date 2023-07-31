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

pdf("Figure_3A_infectionstatus_tsne.pdf",)
DimPlot(cleaned.list$day5, group.by = c("infectionstatus"), reduction = "tsne", cols = c("grey","black"), pt.size = 1)
dev.off()

##Genes upregulated in infected HepOrg cells
infectionmarkers <- FindMarkers(cleaned.list$day5, ident.1 = "yes", ident.2 = "no", only.pos = TRUE, min.pct = 0.25, 
                           thresh.use = 0.25)

write.csv(infectionmarkers,"infection.marker.genes.csv")

##EnrichR analysis
dbs <- listEnrichrDbs()
dbs <- "WikiPathway_2021_Human"

infectionmarker.genes <- rownames(infectionmarkers)
enriched.infected <- enrichr(infectionmarker.genes, dbs)
enriched.infected[["WikiPathway_2021_Human"]]
write.csv(enriched.infected[["WikiPathway_2021_Human"]],"infected_WikiPathway_Human.csv")


