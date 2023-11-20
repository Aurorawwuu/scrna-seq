# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("Seurat")
# BiocManager::install("SeuratObject")
# install.packages("devtools")
# remove.packages("Matrix")
# devtools::install_version("Matrix",version = "1.6-1.1")

library(tidyverse)
library(Seurat)
library(Matrix)
library(dplyr)
library(clustree)

### load data
dir1 <- "/home/aurora/scrna-seq/output/"
dir2 <- "/outs/filtered_feature_bc_matrix/"
files <- dir(dir1)

meta <- read.csv(file = "/home/aurora/scrna-seq/meta.txt", sep = ",")
rownames(meta) <- meta$Run
meta$Group <- substr(meta$source_name, start = 1, stop = 2)

### create seurat object list
seurat_list <- lapply(files, function(file){
  print(file)
  seurat <- Read10X(paste0(dir1,file,dir2))
  
  seurat_obj <- CreateSeuratObject(counts = seurat, project = "project")
  
  seurat_obj@meta.data$SampleID <- file
  seurat_obj@meta.data$GroupID <- meta$source_name[rownames(meta) == file]
  seurat_obj$Group <- meta$Group[rownames(meta) == file]
  return(seurat_obj)
})

### merge seurat obj in the same group
BM <- merge(x=seurat_list[[1]], y=seurat_list[2:5],
            add.cell.ids = c("BM150", "BM152","BM156","BM157","BM158"))
GM <- merge(x=seurat_list[[6]], y=seurat_list[7:13],
            add.cell.ids = c("GM136","GM143","GM144","GM147","GM148",
                             "GM183","GM184a","GM238"))
oral <- merge(BM,GM)
oral <- JoinLayers(oral)
oral[["RNA"]] <-split(oral[["RNA"]], f = oral$Group)

### Quality Control
oral[["percent.mt"]] <- PercentageFeatureSet(oral, pattern = "^MT-")
oral <- subset(oral, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

### Normalization
oral <- NormalizeData(oral)
oral <- FindVariableFeatures(oral, nfeatures = 4000)
save(oral, file = "/home/aurora/scrna-seq/Rdata/preanalysis.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/preanalysis.RData")

### Cell Cycle Scoring
oral <- JoinLayers(oral)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
oral <- CellCycleScoring(oral, 
                         s.features = s.genes,
                         g2m.features = g2m.genes,
                         set.ident = TRUE)
oral[["RNA"]] <-split(oral[["RNA"]], f = oral$Group)

### Pre-Integration analysis
oral <- ScaleData(oral, 
                   vars.to.regress= c("S.Score", "G2M.Score"),
                   features = rownames(oral))
save(oral, file = "/home/aurora/scrna-seq/Rdata/scaledOral.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/scaledOral.RData")
oral <- RunPCA(oral, npcs = 50)

### Pre-integration UMAP
oral_preumap <- RunUMAP(oral, reduction = "pca", dims = 1:50)
# DimPlot(oral_preumap, reduction = "umap",
#         group.by = "Group")
save(oral, file = "/home/aurora/scrna-seq/Rdata/preUMAP.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/preUMAP.RData")

### Integration
oral <- IntegrateLayers(oral, method = RPCAIntegration)
save(oral, file = "/home/aurora/scrna-seq/Rdata/integratedOral.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/integratedOral.RData")

oral <- FindNeighbors(oral, reduction = "pca", dims = 1:50)
oral <- FindClusters(oral, resolution = seq(from = 0.1, to = 0.3,  by = 0.1))
oral <- RunUMAP(oral, reduction = "pca", dims = 1:50)
# clustree(oral)
save(oral, file = "/home/aurora/scrna-seq/Rdata/analysised.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/analysised.RData")

# ### Cell Annotation
# markers <- FindAllMarkers(oral,
#                           only.pos = TRUE,
#                           min.pct = 0.25,
#                           logfc.threshold = 0.75)
