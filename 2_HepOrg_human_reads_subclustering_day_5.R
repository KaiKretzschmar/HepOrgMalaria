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
