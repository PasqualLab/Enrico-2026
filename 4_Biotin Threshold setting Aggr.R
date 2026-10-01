library(Seurat)
library(ggplot2)

all <- readRDS("All_clusters")
UMAPPlot(all)
all@meta.data

##########################################################
############### Check Biotin level noise #################

FeaturePlot(all , features = "Biotin-TotalSeqC")

# I use as noise reference the Biotin lebels on Endogenous T cells
# since it is a population that I expect not having interact with OT-II SrtA+ 

EndoT <- subset(all, cell_type == "Endogenous T cells")
FeaturePlot(EndoT , features = "Biotin-TotalSeqC")

# I remove OT-II cells mislabeled in this cluster by CD45.1 recognition 

endo_data <- GetAssayData(EndoT, assay = "ADT")
cd45.1 <- endo_data["CD45.1-TotalSeqC",]
summary(cd45.1)
CD45.1_info <- ifelse(cd45.1>0, "positive", "negative")
EndoT$cd45 <- CD45.1_info

EndoT@meta.data

realEndo <- subset(EndoT, cd45 == "negative")
FeaturePlot(realEndo , features = "Biotin-TotalSeqC")
realEndo@meta.data

summary(EndoT$Biotin_data)
summary(all$Biotin_data)
summary(realEndo$Biotin_data)                  # the max contain the highest level of backround 
threshold <- quantile(realEndo$Biotin_data, probs = 0.95)   # I use as threshold the 95th percentile value 
                                                            # of the distribution of the Biotin signal inside a ctrl negative population
threshold
# 1.094818

# We set the threshold to 1 since it was not too far from the calculated noise threshold and we
# wanted to gain some positive cells for the 24h timepoint analysis with the scSeqCommDiff, that otherwise they would be too few

threshold = 1


#####################################################################
###################### Setting Biotin threshold ####################

all_data <- GetAssayData(all, assay = "ADT")
biotin_data <- all_data["Biotin-TotalSeqC",]
boolean <- biotin_data > threshold

biotin_threshold <- ifelse(boolean, "positive", "negative")
table(biotin_threshold)

all$Biotin_status <- biotin_threshold

View(all@meta.data)

DimPlot(all, group.by = "Biotin_status", cols = c("grey", "red"))

DimPlot(all, group.by = "TimepointGrouped", cols = c("blue", "grey"))



####################################################################################
###################### Check Biotin threshold in subpopulations ####################

mdc <- subset(all, cell_type == "mDCs")
h24 <- subset(mdc, TimepointGrouped == "24h") 
table(mdc$TimepointGrouped)
table(h24$Biotin_status) # check the biotin+ cells in the 24h timepoint mDCs for the scSeqCommDiff analysis

RidgePlot(all, features = "Biotin-TotalSeqC")
RidgePlot(EndoT, features = "Biotin-TotalSeqC")
VlnPlot(EndoT, features = "Biotin-TotalSeqC")



Idents(all) <- "cell_type"
RidgePlot(all, "Biotin-TotalSeqC")


DefaultAssay(all) <- "RNA"
all_d <- GetAssayData(all, assay = "ADT")

Idents(all) <- "cell_type"
DefaultAssay(all) <- "ADT"
VlnPlot(all, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") 

DefaultAssay(all) <- "RNA"

all@meta.data
all 

saveRDS(all, file = "All_clusters")
