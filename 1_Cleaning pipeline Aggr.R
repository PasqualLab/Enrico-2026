
##############################################################
##### CLEANING PIPELINE FOR DC - T CELL CITE-SEQ SAMPLES #####
##############################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

raw_S <- Read10X("filtered_feature_bc_matrix")
summary(raw_S)

View(raw_S)
raw_S$`Gene Expression`
raw_S$`Antibody Capture` 

##### "s" is my main object to store GEX and ADT data #####
s <- CreateSeuratObject(raw_S$`Gene Expression`)
s



##### Create two assays to store ADT and HTO data #####
Hashtag_assay <-  CreateAssay5Object(counts = raw_S$`Antibody Capture`[1:12,]) #Hastags are stored in the first 12 row of the "Antibody Captured" dataset
s[["HTO"]] <- Hashtag_assay

ADT_assay <- CreateAssay5Object(counts = raw_S$`Antibody Capture`[13:36,]) #ADTs are stored in the next 24 row of the "Antibody Captured" dataset
s[["ADT"]] <- ADT_assay
s

rownames(s[["HTO"]])
rownames(s[["ADT"]])
s@meta.data

#####################################################################
##### Empty droplets, lowRNA highRNA cells and highMito removal #####

range(s$nCount_RNA) # mRNAs per cell

range(s$nFeature_RNA) # genes per cell


s <- PercentageFeatureSet(s, pattern = "^mt-", col.name = "percent.mt")
s@meta.data
range(s$percent.mt)


pm1 <- VlnPlot(s, features = c("nFeature_RNA")) +
  scale_y_continuous(breaks = seq(0, 6000, by = 500)) +
  geom_hline(yintercept = 5500, size = 0.25, color = "black") +
  geom_hline(yintercept = 500, size = 0.25, color = "blue")

pm2 <- VlnPlot(s, features = c("nCount_RNA")) +
  scale_y_continuous(breaks = seq(0, 40000, by = 2000)) +
  geom_hline(yintercept = 30000, size = 0.25)

pm3 <- VlnPlot(s, features = c("percent.mt")) +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  geom_hline(yintercept = 2, size = 0.25)

library(cowplot)
plot_grid(pm1, pm2, pm3, 
          labels = c("A", "B", "C"),
          ncol = 3)

pm4 <- FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "percent.mt")
pm5 <- FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pm4 | pm5

s 
#s pre-clean: 28232

s <- subset(s, subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & nCount_RNA < 30000 & percent.mt < 2)
s
# s post-clean: 26715

range(s$percent.mt)


##############################################################################
############ Normalization and scaling before doublets removal ###############

s <- NormalizeData(s)
s <- FindVariableFeatures(s, nfeatures = 2500)
s <- ScaleData(s)
s <- RunPCA(s)

ElbowPlot(s, 50) # defining PCA number
DimHeatmap(s, dims = 15:30, cells = 500, balanced = TRUE)
# With Elbowplot the inflection point is around PC 25 

# 80% rule:
variance <- (s@reductions$pca@stdev)^2
var80 <- sum(variance)*0.8
findpca <- which(cumsum(variance) >= var80)[1]
findpca 


# Uso le prime 26 PC per il clustering

s <- FindNeighbors(s, dims = 1:26)
s <- FindClusters(s)

s <- RunUMAP(s, dims = 1:26)
DimPlot(s, reduction = "umap")

######################################################################
##### Doublet removal with DoubletFinder #####

library(DoubletFinder) # https://github.com/chris-mcginnis-ucsf/DoubletFinder ; https://www.sciencedirect.com/science/article/pii/S2405471219300730
# DoubletFinder simulates artificial doublets from existing scRNA-seq data 
# by averaging the gene expression profiles of random pairs of cells

DefaultAssay(s) <- "RNA"
s
 
# three input parameters:
# 1) the number of expected real doublets
# 2) the number of artificial doublets (pN)
# 3) the neighborhood size (pK)

sweep.list <- paramSweep(s, PCs = 1:26, sct=FALSE) # Parameters tested: pN = 0.05-0.3, pK = 0.0005-0.3. 
sweep.stats <- summarizeSweep(sweep.list , GT = F) # Returns a table with key metrics for each tested pK 
optPKV <- find.pK(sweep.stats) # Find the optimal parameter for calculating the doublets

ggplot(optPKV, aes(pK, BCmetric, group = 1)) +
  geom_point() + 
  geom_line()

pK <- optPKV %>%
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
pK <- as.numeric(as.character(pK[[1]]))
pK


annot <- s@meta.data$seurat_clusters
homo <- modelHomotypic(annot) # estimates the proportion of homotypic doublets (same or very similar cell type)
exp.db <- round(0.08*nrow(s@meta.data)) # Calculate the expected fraction of doublets [the percentage is taken from 10x website based on the total number of loaded cells in the run]
exp.db.adj <- round(exp.db*(1- homo))  # Adjust the number of doublets on the homotypic doublets proportion

s <- doubletFinder(s, PCs = 1:26, pN = 0.25, pK = pK, nExp = exp.db.adj, reuse.pANN = NULL, sct = F) # pN is the proportion of cells considered for creating the artificial doublets (25% Default)
s@meta.data
table(s$DF.classifications_0.25_0.06_1947) 

# Doublets visualization
FeaturePlot(s, features = "pANN_0.25_0.06_1947")
DimPlot(s, reduction = "umap", group.by = "DF.classifications_0.25_0.06_1947")

# Doublets removal
s <- subset(s, subset = DF.classifications_0.25_0.06_1947 == "Singlet")
s

# 33732 features across 24768 samples within 3 assays 



#####################################################################
###################### Cleaninng by Hashtag  ############################

# Normalization
DefaultAssay(s) <- "HTO"
s<- NormalizeData(s, normalization.method = "CLR", margin = 2, assay = "HTO")
s

# HTODemux demultiplexes automatically assigning max hastag to every cell
s <- HTODemux(s, assay = "HTO", positive.quantile = 0.99)
s@meta.data
s

HTOHeatmap(sm, assay = "HTO")
# Cells without or with more than one hashtag are removed
s <- subset(s, subset = HTO_classification.global == "Singlet")
s
# 21260 cells with one and only one Hashtag


DimPlot(s, group.by = "seurat_clusters")



# Object cleaning 
s$pANN_0.25_0.06_1947 <- NULL
s$DF.classifications_0.25_0.06_1947 <- NULL
s$HTO_maxID <- NULL
s$HTO_secondID <- NULL
s$HTO_margin <- NULL
s$HTO_classification.global <- NULL
s$HTO_classification <- NULL


View(s@meta.data)

saveRDS(s, file = "Cleaned_DS")

rownames(s[["HTO"]])
rownames(s[["ADT"]])

