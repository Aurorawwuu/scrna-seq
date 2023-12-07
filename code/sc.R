# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("Seurat")
# BiocManager::install("SeuratObject")
# install.packages("devtools")
# remove.packages("Matrix")
# devtools::install_version("Matrix",version = "1.6-1.1")
# devtools::install_github('immunogenomics/presto')
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# devtools::install_version("dbplyr", version = "2.3.4")
# BiocManager::install("SingleCellExperiment")

library(tidyverse)
library(Seurat)
library(Matrix)
library(dplyr)
library(clustree)
library(presto)
library(SingleR)
library(celldex)
library(dbplyr)
library(SingleCellExperiment)
library(paletteer)
library(ggpubr)
library(reshape2)


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

### preqc vs postqc for visualization
preqc_BM <- BM
preqc_GM <- GM
preqc_BM[["percent.mt"]] <- PercentageFeatureSet(preqc_BM, pattern = "^MT-")
# save(preqc_BM, file = "/home/aurora/scrna-seq/Rdata/preqcBM.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/preqcBM.RData")
preqc_GM[["percent.mt"]] <- PercentageFeatureSet(preqc_GM, pattern = "^MT-")
# save(preqc_GM, file= "/home/aurora/scrna-seq/Rdata/preqcGM.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/preqcGM.RData")
postqc_BM <- subset(preqc_BM, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
# save(postqc_BM, file = "/home/aurora/scrna-seq/Rdata/postqcBM.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/postqcBM.RData")
postqc_GM <- subset(preqc_GM, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
# save(postqc_GM, file = "/home/aurora/scrna-seq/Rdata/postqcGM.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/postqcGM.RData")

oral <- merge(BM,GM)
oral <- JoinLayers(oral)
oral[["RNA"]] <-split(oral[["RNA"]], f = oral$Group)


### Quality Control
oral[["percent.mt"]] <- PercentageFeatureSet(oral, pattern = "^MT-")
oral <- subset(oral, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)


### Normalization
oral <- NormalizeData(oral)
oral <- FindVariableFeatures(oral, nfeatures = 4000)
# save(oral, file = "/home/aurora/scrna-seq/Rdata/preanalysis.RData")
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
# save(oral, file = "/home/aurora/scrna-seq/Rdata/scaledOral.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/scaledOral.RData")
oral <- RunPCA(oral, npcs = 50)


## Pre-integration UMAP
oral_preumap <- RunUMAP(oral, reduction = "pca", dims = 1:50)
preumap <- DimPlot(oral_preumap, reduction = "umap",
        group.by = "Group")
# save(oral_preumap, file = "/home/aurora/scrna-seq/Rdata/preUMAP.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/preUMAP.RData")


### Integration + Cell Annotation
oral <- IntegrateLayers(oral, method = CCAIntegration, 
                        orig.reduction = "pca", new.reduction = "integrated.cca")
# save(oral, file = "/home/aurora/scrna-seq/Rdata/integratedOral.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/integratedOral.RData")


### Post Integration UMAP
oral <- JoinLayers(oral)
oral_postumap <- RunUMAP(oral, reduction = "integrated.cca", dims = 1:50)
postumap <- DimPlot(oral_postumap, reduction = "umap",
                    group.by = "Group")
preumap+postumap
# save(oral_postumap, file = "/home/aurora/scrna-seq/Rdata/postUMAP.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/postUMAP.RData")


### Post Integration Clustering
oral <- FindNeighbors(oral, reduction = "integrated.cca", dims = 1:50)
oral <- FindClusters(oral, resolution = seq(from = 0.1, to = 0.3,  by = 0.1))
oral <- RunUMAP(oral, reduction = "integrated.cca", dims = 1:50)
clustree(oral)
# save(oral, file = "/home/aurora/scrna-seq/Rdata/analysised.RData")
load(file = "/home/aurora/scrna-seq/Rdata/analysised.RData")


### Define Marker Genes
DefaultAssay(oral) <- "RNA"
Idents(oral) <- "RNA_snn_res.0.3"
markers <- FindAllMarkers(oral,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.75)
significant <- markers[markers$p_val_adj < 0.2, ]
DimPlot(oral, group.by = c("Group", "RNA_snn_res.0.3"), label = T)
# save(significant, file = "/home/aurora/scrna-seq/Rdata/significant.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/significant.RData")


### Cell Annotation
neworal <- subset(oral, idents = c("4","11", "17", "18"),
                  invert = T)
neworal <- RenameIdents(neworal,
                        "0"="endothelial",
                        "6"="endothelial",
                        "1"="fibroblast",
                        "3"="fibroblast",
                        "2"="immune",
                        "5"="immune",
                        "8"="immune",
                        "9"="immune",
                        "10"="immune",
                        "13"="immune",
                        "14"="immune",
                        "15"="immune", 
                        "7"="epithelial",
                        "12"="epithelial",
                        "16"="other")
neworal$GeneralCellType <- neworal@active.ident
# save(neworal, file = "/home/aurora/scrna-seq/Rdata/neworal.RData")
# load(file = "/home/aurora/scrna-seq/Rdata/neworal.RData")
postanno <- DimPlot(neworal, group.by = "GeneralCellType", label = T)
preanno <- DimPlot(oral, group.by = "RNA_snn_res.0.3", label = T)
postanno+preanno


### plot 1B
load(file = "/home/aurora/scrna-seq/Rdata/neworal.RData")

pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,6)]

groupprop <- ggplot(as.data.frame(neworal@meta.data), aes(Group, fill=GeneralCellType))+
  geom_bar(position="fill")+
  theme(panel.background = element_rect(NA),
        axis.ticks = element_line(NA),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"),
        aspect.ratio = 3/1)+
  ylab("Proportion")+
  scale_fill_manual(values=pal)+
  guides(fill=F)+
  annotate("text", x=1, y=1.07, label = "16,292\ncells")+
  annotate("text", x=2, y=1.07, label = "28,942\ncells")+
  scale_x_discrete(labels=c("B", "G"), breaks=c("BM", "GM"))
groupprop

umap_celltype <- DimPlot(neworal, 
                         group.by = "GeneralCellType",
                         cols=pal,
                         pt.size = 0.4)+
  ylab("UMAP2")+
  xlab("UMAP1")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",
        plot.title = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0.005))+
  annotate("text", x=-10, y=-14, label="45,234 cells")
umap_celltype

ggarange(umap_celltype, groupprop, common.legend = T)



### supplemental A
load(file = "/home/aurora/scrna-seq/Rdata/preqcBM.RData")
load(file = "/home/aurora/scrna-seq/Rdata/preqcGM.RData")
load(file = "/home/aurora/scrna-seq/Rdata/postqcBM.RData")
load(file = "/home/aurora/scrna-seq/Rdata/postqcGM.RData")

qc <- c(preqc_BM, preqc_GM, postqc_BM, postqc_GM)
for ( i in 1:4){
  qc[[i]]<-VlnPlot(qc[[i]],
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   group.by = "GroupID",
                   pt.size = 0,
                   assay = "RNA")
}
ggarrange(qc[[1]], qc[[2]], qc[[3]], qc[[4]], nrow = 4)


### Supplemental C
preqc_BM$qc <- "preqc_BM"
postqc_BM$qc <- "postqc_BM"
BMmerge <- merge(preqc_BM, postqc_BM)
preqc_GM$qc <- "preqc_GM"
postqc_GM$qc <- "postqc_GM"
GMmerge <- merge(preqc_GM, postqc_GM)

BMcellnumb <- ggplot(as.data.frame(BMmerge@meta.data), 
       aes(x=GroupID,
           fill=factor(qc, levels = c("preqc_BM", "postqc_BM"))))+
  geom_bar(position = "dodge")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_fill_manual(labels=c("preQC", "postQC"), 
                      values = c("grey", "black"))

GMcellnumb <- ggplot(as.data.frame(GMmerge@meta.data), 
       aes(x=GroupID,
           fill=factor(qc, levels = c("preqc_GM", "postqc_GM"))))+
  geom_bar(position = "dodge")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c("grey", "black"))+
  guides(fill=F)

ggarrange(BMcellnumb, GMcellnumb, nrow=2, 
          common.legend = T, legend = "right")


### 

