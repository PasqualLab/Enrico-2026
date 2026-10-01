#######################################################
##### SAMPLE NORMALIZAZION AND HASHTAG ASSIGNMENT #####
#######################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)


########## LOAD OBJECTS ###########

sm <- readRDS("Cleaned_DS")

##########################################################
########## DIMRED, NORMALIZATION AND CLUSTERING ##########

# Normalization and clustering

DefaultAssay(sm) <- "RNA"

sm <- FindVariableFeatures(sm, nfeatures = 2500)
sm <- ScaleData(sm)
sm <- RunPCA(sm, verbose = T)
ElbowPlot(sm, 50) # defining PCA number

VariableFeaturePlot(object = sm)
VariableFeatures(sm)
top10 <- head(VariableFeatures(sm), 20)

LabelPoints(plot = VariableFeaturePlot(object = sm), points = top10, repel = TRUE)



# clustering with the first 30 PCs with 0.8 resolution
sm <- FindNeighbors(sm, dims = 1:30)
sm <- FindClusters(sm, resolution = 0.8, verbose = FALSE)
sm <- RunTSNE(sm, dims = 1:30, verbose = T)
DimPlot(sm, reduction = "tsne")

sm <- RunUMAP(sm, dims = 1:30)
DimPlot(sm, reduction = "umap")

########## ADT Normalization #########

DefaultAssay(sm) <- "ADT"
sm <- NormalizeData(sm, normalization.method = "CLR", margin = 2)

####### Hashtag identification across clusters #####

DefaultAssay(sm) <- "HTO"

sm@meta.data

# Batch effect evaluation for library and sample
DimPlot(sm, group.by = "hash.ID", reduction = "pca")
# No evident Batch effect by individual



# Number of retrieved cells cells for each time point
# in our experimental setup mice  1 to 6 were analyzed at 24h p.i., mice 7 to 12 at 72h p.i.


sum(sm@meta.data$hash.ID == "Hashtag-12-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-11-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-10-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-09-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-08-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-07-TotalSeqC")
# 16871 cells deriving from 72h samples
sum(sm@meta.data$hash.ID == "Hashtag-01-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-02-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-03-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-04-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-05-TotalSeqC") + sum(sm@meta.data$hash.ID == "Hashtag-06-TotalSeqC")
# 4389 cells deriving from 24h samples


###### Labeling cells by time-point origin #####

Boolean_Hash <- sm@meta.data$hash.ID %in% c("Hashtag-01-TotalSeqC", "Hashtag-02-TotalSeqC", "Hashtag-03-TotalSeqC","Hashtag-04-TotalSeqC","Hashtag-05-TotalSeqC","Hashtag-06-TotalSeqC")
table(Boolean_Hash)

# Boolean_Hash
# FALSE  TRUE 
# 16744  4424 

# cell time point labeling to save in the object
label_boolean <- ifelse(Boolean_Hash, "24h", "72h")
table(label_boolean)

# label_boolean
# 24h   72h 
# 4389 16871 

sm$TimepointGrouped <- label_boolean # Add the information to the obj

DimPlot(sm, group.by = "TimepointGrouped")
DimPlot(sm, group.by = "TimepointGrouped", reduction = "tsne")


##########################################
########## Save Seurat Object sm ##########
DefaultAssay(sm) <- "RNA"
DefaultAssay(sm)


saveRDS(sm, file = "Cleaned_Timepoint_DS")


