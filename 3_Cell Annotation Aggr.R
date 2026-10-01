##########################################
########## Exploratory analysis ##########
library(Seurat)

sm <- readRDS("Cleaned_Timepoint_DS")
sm
DefaultAssay(sm) <- "RNA"
# For the cell identification I used: 
# 1. automatic annotation with singleR on the  Immgen consortium dataset
# 2. manual annotation exploiting GEX and ADT data 


##################################################
######## Automatic cell type annotation ##########


library(SingleR)
# load ImmGenData reference
ref <- celldex::ImmGenData()
#View(as.data.frame(colData(ref)))

library(SingleCellExperiment)
# Since SingleR works on SingleCellExperiment object, I transform the Seurat obj in a SCE obj
sce <- as.SingleCellExperiment(sm, assay = "RNA")

auto.ann <- SingleR(sce, ref = ref, labels = ref$label.main) # Rough annotation
# auto.ann_fine <- SingleR(sce, ref = ref, labels = ref$label.fine) # fine annotation. Not used in the analysis

table(auto.ann$pruned.labels)
#View(table(auto.ann_fine$pruned.labels))

plotScoreHeatmap(auto.ann)
# plotScoreHeatmap(auto.ann_fine)


sm$singleR.labels <- auto.ann$labels[match(rownames(sm@meta.data), rownames(auto.ann))]
# sm$singleR.labels_fine <- auto.ann_fine$labels[match(rownames(sm@meta.data), rownames(auto.ann_fine))]
DimPlot(sm, group.by = "singleR.labels")



###########################################
###### Manual cell type annotation #####

DefaultAssay(sm) <- "RNA"
cDCs <- c("Itgax", "Zbtb46", "Cd209a", "H2-Ab1", "Dpp4") # common DC markers
cDC1 <- c("Itgax", "Cd8a", "Xcr1", "Clec9a", "Irf8", "Itgae") 
cDC2 <- c("Sirpa", "Itgam", "Irf4", "Cd4", "Mgl2", "Nod2", "Klf4", "Cx3cr1", "Esam")
Mo <- c("C5ar1", "Ccr2", "Ly6c1", "Sell", "Adgre1", "Fcgr1", "Mertk") # Monocye markers
mregDCs <- c("Cd200", "Cd63", "Icosl", "Btla", "Dll4", "Alcam", "Ccr7", "Cd5") # mature regulatory DCs 

Tcells <- c("Cd3e", "Cd8a", "Cd4", "Cd28", "Ctla4")
TF <- c("Tbx21", "Gata3", "Rorc", "Stat6", "Stat1", "Stat4", "Bcl6", "Batf") # Transcriptional Factor helper T cells

FeaturePlot(sm, features = cDCs)
FeaturePlot(sm, features = cDC1)
FeaturePlot(sm, features = cDC2)
FeaturePlot(sm, features = Mo)
FeaturePlot(sm, features = mregDCs)

FeaturePlot(sm, features = Tcells)
FeaturePlot(sm, features = TF)


##### Visualizing ADT information #####

DefaultAssay(sm) <- "ADT"
cDCs_ADT <- c("CD11c-TotalSeqC", "IA-IE-TotalSeqC", "Sirpa-TotalSeqC", "XCR1-TotalSeqC", "CD8a-TotalSeqC", "CD11b-TotalSeqC", "CD103-TotalSeqC", "CD207-TotalSeqC", "CD301b-TotalSeqC")
Mo_ADT <- c("Ly6C-TotalSeqC", "F4-80-TotalSeqC")
T_ADT <- c("CD3-TotalSeqC", "CD45.1-TotalSeqC", "CD4-TotalSeqC", "Va2-TotalSeqC", "Vb5-TotalSeqC", "CD62L-TotalSeqC", "CD44-TotalSeqC", "CD25-TotalSeqC", "CD69-TotalSeqC")
  
FeaturePlot(sm, features = cDCs_ADT)
FeaturePlot(sm, features = Mo_ADT)
FeaturePlot(sm, features = "Biotin-TotalSeqC")
FeaturePlot(sm, features = T_ADT)

###############################################################################
###### Creating two different object splitting T cells from Myeloid cells #####

UMAPPlot(sm, label = T)
DimPlot(sm, group.by = "TimepointGrouped")
Idents(sm) <- "seurat_clusters"
clusters <- sm$seurat_clusters

cell_type_annotation <- c(
  "0" = "Monocytes",
  "1" = "cDC2s",
  "2" = "OT-II",
  "3" = "mDCs",
  "4" = "OT-II",
  "5" = "OT-II",
  "6" = "OT-II",
  "7" = "mDCs",
  "8" = "OT-II",
  "9" = "cDC1s",
  "10" = "Endogenous T cells",
  "11" = "CD8+ DCs",
  "12" = "mDCs",
  "13" = "Doublets",
  "14" = "Ly6D+ cells",
  "15" = "B cells"
)

sm$cell_type <- unname(cell_type_annotation [as.character(sm$seurat_clusters)])

DimPlot(sm, group.by = "cell_type", label = T)


########################################################################
################ Add Biotin_levels column #############################

DefaultAssay(sm) <- "ADT"
adt_data <- GetAssayData(sm, layer = "data") 
Biotin_data <- adt_data["Biotin-TotalSeqC",]

sm$Biotin_data <- Biotin_data
sm@meta.data

#########################################################################
################### Further cleainig of sm object #######################

# B cells and doublets removal

sm.clean <- subset(sm, cell_type %in% c("B cells", "Doublets"), invert = T)
DefaultAssay(sm.clean) <- "RNA"
DimPlot(sm.clean, group.by = "cell_type", label = T)


#######################################################
########## Save Seurat Object sm.T and sm.dc ##########
saveRDS(sm.clean, file = "All_clusters")

