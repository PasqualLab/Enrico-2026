library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)


#################################################################################
########################### mature DCs analysis #################################

mdc <- readRDS("mDCs_clusters")

DimPlot(mdc, group.by = "cell_type")

mdc
# 33732 features across 3487 samples within 3 assays

DefaultAssay(mdc) <- "RNA"
DimPlot(mdc, group.by = "seurat_clusters", label = T)

RidgePlot(mdc, features = "Biotin-TotalSeqC")

FeaturePlot(mdc, features = "Biotin-TotalSeqC")

mdc@meta.data

### Exploratory analysis of Biotin expression 

#Extract Biotin data from mdc object
Biotin_data <- FetchData(mdc, vars = "Biotin-TotalSeqC")[, 1]

VlnPlot(mdc, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = (quantile(Biotin_data, probs = 0.80)), linetype = "dashed", color = "red") +
  geom_hline(yintercept = (quantile(Biotin_data, probs = 0.20)), linetype = "dashed", color = "blue")

DimPlot(mdc, group.by = "Biotin_status", cols = c("grey", "blue", "red"))

###################################################################
###################### Manual Annotation ##########################

# Manual annotation exploiting the surface marker expression 

IDplotBarcode <- FeaturePlot(mdc, features = c("XCR1-TotalSeqC", "CD8a-TotalSeqC" ,"Sirpa-TotalSeqC", "CD11b-TotalSeqC", "CD301b-TotalSeqC", "CD103-TotalSeqC"))
IDplotGEX <- FeaturePlot(mdc, features = c("Xcr1", "Cd8a" ,"Sirpa", "Itgam", "Mgl2", "Itgae"))
IDplotBarcode | UMAPPlot(mdc, label = T) | IDplotGEX

Timeplot <- DimPlot(mdc, group.by = "TimepointGrouped")

Timeplot | UMAPPlot(mdc, label = T) | FeaturePlot(mdc, features = "Biotin-TotalSeqC")



# Dot plot fot mDCs cluster identification
DefaultAssay(mdc) <- "ADT"
ADT_feat <- rownames(mdc)
rownames(mdc)
ADT_mDCs <- ADT_feat[!ADT_feat %in% c("Va2-TotalSeqC", "Vb5-TotalSeqC", "CD45.1-TotalSeqC", "CD3-TotalSeqC", "CD62L-TotalSeqC", "CD44-TotalSeqC", "CD69-TotalSeqC","CD25-TotalSeqC", "Ly6G-TotalSeqC", "NK1.1-TotalSeqC", "B220-TotalSeqC")]

mdc_DotPlot <- DotPlot(mdc, features = ADT_mDCs) + RotatedAxis() +  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdYlBu"))) 
ggsave(filename = "mDCs_ADTs.svg", plot = p11, dpi = 300, width = 7, height = 5)
DefaultAssay(mdc) <- "RNA"
DimPlot(mdc, group.by = "seurat_clusters", label = T)



###############################################################################
################### DGE analysis across all mDCs clusters #####################

DefaultAssay(mdc) <- "RNA"
Idents(mdc)
DGE_mdc <- FindAllMarkers(mdc)
write.csv(DGE_mdc, file = "mDCs_AllvsAll_DEG.csv")



####################################################################################
######################## mDCs distribution at 24h and 72h ##########################

mdc24h <- subset(mdc, TimepointGrouped == "24h")
mdc72h <- subset(mdc, TimepointGrouped == "72h")

DimPlot(mdc24h, group.by = "seurat_clusters", label = T)
pie(table(mdc24h$seurat_clusters))

DimPlot(mdc72h, group.by = "seurat_clusters", label = T)
pie(table(mdc72h$seurat_clusters))

ratio24vs72 <- table(mdc24h$seurat_clusters) / table(mdc72h$seurat_clusters)


# Ratio of the mDCs between 24h and 72h
ggplot(
  data.frame(
    Cluster = factor(names(ratio24vs72)),
    Ratio = as.numeric(ratio24vs72)
  ),
  aes(x = Cluster, y = Ratio)
) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(
    x = "Seurat cluster",
    y = "24h / 72h ratio",
    title = "Cluster composition: 24h vs 72h"
  )
# Cluster 8 is the only one more abundant at 24h p.i. 


###############################################################################
######################## Analysis of interacting mDCs ##########################

### % of Biotin+ cell across clusters at 24h p.i.


Biotin_mdc24 <- FetchData(mdc24h, vars = c("Biotin_status", "seurat_clusters"))  

Biotin_mdc24

# Percentage of Biotin+ cells per cluster

Biotin_mdc24percent <- Biotin_mdc24 %>%
  group_by(seurat_clusters) %>%
  summarise(
    Total = n(),
    Biotin_pos = sum(Biotin_status == "positive"),
    Percentage = Biotin_pos / Total * 100,
    .groups = "drop"
  )

# Plot
ggplot(Biotin_mdc24percent,
       aes(x = factor(seurat_clusters), y = Percentage)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20)
  ) +
  theme_classic() +
  labs(
    title = "Biotin positivity at 24h",
    x = "Seurat cluster",
    y = "Biotin+ cells (%)"
  )



### % of Biotin+ cell across clusters at 72h p.i.

Biotin_mdc72 <- FetchData(mdc72h, vars = c("Biotin_status", "seurat_clusters"))  

Biotin_mdc72

# Percentage of Biotin+ cells per cluster

Biotin_mdc72percent <- Biotin_mdc72 %>%
  group_by(seurat_clusters) %>%
  summarise(
    Total = n(),
    Biotin_pos = sum(Biotin_status == "positive"),
    Percentage = Biotin_pos / Total * 100,
    .groups = "drop"
  )

# Plot
ggplot(Biotin_mdc72percent,
       aes(x = factor(seurat_clusters), y = Percentage)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20)
  ) +
  theme_classic() +
  labs(
    title = "Biotin positivity at 24h",
    x = "Seurat cluster",
    y = "Biotin+ cells (%)"
  )




### Cluster composition of Biotin+ DC at 24h p.i.


# calculate the % of clusters in all the biotin+ cells
Biotin_mdc24_composition <- Biotin_mdc24 %>%
  filter(Biotin_status == "positive") %>%
  group_by(seurat_clusters) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Percentage = Count / sum(Count) * 100
  )

#Histogram
ggplot(Biotin_mdc24_composition,
       aes(x = factor(seurat_clusters), y = Percentage)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20)
  ) +
  theme_classic() +
  labs(
    title = "Composition of Biotin+ DCs at 24h",
    x = "Seurat cluster",
    y = "Percentage of total Biotin+ DCs"
  )

# Pie 
ggplot(Biotin_mdc24_composition,
       aes(x = "", y = Count,
           fill = factor(seurat_clusters))) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  theme_void() +
  labs(
    title = "Biotin+ DC composition at 24h",
    fill = "Cluster"
  )



### Cluster composition of Biotin+ DC at 72h p.i.


# calculate the % of clusters in all the biotin+ cells
Biotin_mdc72_composition <- Biotin_mdc72 %>%
  filter(Biotin_status == "positive") %>%
  group_by(seurat_clusters) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Percentage = Count / sum(Count) * 100
  )

#Histogram
ggplot(Biotin_mdc72_composition,
       aes(x = factor(seurat_clusters), y = Percentage)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20)
  ) +
  theme_classic() +
  labs(
    title = "Composition of Biotin+ DCs at 24h",
    x = "Seurat cluster",
    y = "Percentage of total Biotin+ DCs"
  )

# Pie 
ggplot(Biotin_mdc72_composition,
       aes(x = "", y = Count,
           fill = factor(seurat_clusters))) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  theme_void() +
  labs(
    title = "Biotin+ DC composition at 24h",
    fill = "Cluster"
  )


#########################################################################
######################## Correlation analysis  ##########################
# mDC Biotin+
mdc72h_positive <- subset(mdc72h, Biotin_status == "positive")

# gene expression and biotin
gene_expr <- GetAssayData(mdc72h_positive, "RNA")
biotin <- Seurat::GetAssayData(mdc72h_positive,assay = "ADT",layer = "data")["Biotin-TotalSeqC", ]

# Spearman correlation
library(doParallel)
doParallel::registerDoParallel(cores = 5)
res <- foreach(i = 1:nrow(gene_expr), .combine = "rbind") %dopar%
{
    r <- cor.test(biotin, gene_expr[i,], method = "spearman")
    data.frame(corr = r$estimate, pval = r$p.value)
}
row.names(res) <- row.names(gene_expr)
res$gene <- row.names(res)

# remove NA, keep significant correlations (after FDR) and  order by correlation
res <- res[ !(is.na(res$corr)),]
res <- res[order(res$corr),]
res$order <- 1:nrow(res)
res$pval_adj <- p.adjust(res$pval, "fdr")
res <- res[res$pval_adj < 0.05,]

# retreive LR pairs from CellTalkDB
?scSeqComm::LR_pairs_Shao_2020_mouse
LRpairs <- scSeqComm::LR_pairs_Shao_2020_mouse
R <- unique(LRpairs$receptor) # receptors
L <- unique(LRpairs$ligand) # ligands

# filter for ligand and receptor genes and positive correlations
LRres <- res[res$gene %in% c(R, L),]

# plot
ggplot(LRres[LRres$corr>0,] , aes(x = order ,y = corr, color = gene %in% R))+ geom_point() +
  scale_color_manual(values = c("TRUE" = "#3288BD" , "FALSE" = "#F46D43"),
                     name = "Gene type",                     # legend title
                     labels = c("FALSE" = "Ligand", "TRUE" = "Receptor"))+
  ggrepel::geom_text_repel(
    data = LRres[LRres$corr>0,],
    aes(label = gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"), show.legend = FALSE
  ) + xlab("ligand/receptor genes") + ylab("Spearman coefficient") + ggtitle("mDCs_positive")+ theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid = element_blank())
