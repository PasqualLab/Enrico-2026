library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)

all <- readRDS("All_clusters")
all@meta.data
Idents(all) <- "seurat_clusters"
UMAPPlot(all, label= T)


#############################################################################################
############################ Myeloid cells deeper classification ############################

my <- subset(all, cell_type %in% c("Endogenous T cells", "OT-II"), invert = T)
UMAPPlot(my, label = T)
DimPlot(my, group.by = "Biotin_status", cols = c("grey", "red"))
DimPlot(my, group.by = "cell_type")

####Further cell annotation within myeloid clusters

# CD8a+ SIRPa+ DCs are described as tDCs here : Minutti CM et al. Distinct ontogenetic lineages dictate cDC2 heterogeneity. Nat Immunol. 2024 Mar;25(3):448-461. doi: 10.1038/s41590-024-01745-9. Epub 2024 Feb 13. PMID: 38351322; PMCID: PMC10907303.
# Veryfy the overlap with the transcriptional signature obtained from:
# Sulczewski FB et al. Transitional dendritic cells are distinct from conventional DC2 precursors and mediate proinflammatory antiviral responses. Nat Immunol. 2023 Aug;24(8):1265-1280. doi: 10.1038/s41590-023-01545-7. Epub 2023 Jul 6. PMID: 37414907; PMCID: PMC10683792.

tDCs_signature <- c(
  "Ngfr", "Rab3il1", "Alox5", "Cyp7b1", "Fkbp10", "Upk1b", "Cd209e", "Zbtb32", 
  "Trim16", "Ednrb", "Arhgap23", "Bhlhe22", "Cp", "Pros1", "Mrvi1", "Arhgef17", 
  "Aqp3", "Ccr1", "Fam78b", "Cd3g", "Cnn3", "Il2rb", "Chil5", "Tppp3", "Ly6i", 
  "Gimap3", "Cd200r1", "Crip2", "Pi16", "2510016D11Rik", "Dkkl1", "1500009L16Rik", 
  "P3h4", "Lmo1", "Ptgr1", "Emp1", "Lef1", "Slc34a1", "L1cam", "Themis", "Gm38244", 
  "Gm45552", "Bcl11b", "Cd3d", "Cd200r4", "Olr1", "Ppp1r3d", "Trbc1", "Trac", "Cd6", 
  "Gprc5c", "Rasd1", "Plvap", "Ckap4", "Cdr2l", "Gm49553", "Palld", "Tpcn1", "Hopx", 
  "Angpt1", "Otof", "Gm33973", "Dusp3", "S1pr1", "Itk", "Gpr156", "Paqr3", "Tiam2", 
  "Dpysl2", "Kcnk10", "Ccnd3", "Micu1", "Ccr2", "Klf12", "Zdhhc15", "Fgfr1", "Vim", 
  "Plscr3", "Gramd4", "Rwdd2a", "Tep1", "Stat5a", "Fry", "Ttc3", "Gstp3", "Gm12195", 
  "Tle1", "Mob3b", "Csad", "Sfxn3", "Rpusd2", "Bmp2k", "Fxyd5", "Ap2a2"
)

#Featureplot  with tDC signature
MStDCs_my <- AddModuleScore(my, features = list(tDCs_signature), name = "tDC_signature")
plot_tDCs <- FeaturePlot(MStDCs_my, features = "tDC_signature1", slot = "data", cols = rev(c("#F46D43","white","#3288BD")))


# Ly6D+ DCs defined as pre-DCs here: Dress RJ et al. Plasmacytoid dendritic cells develop from Ly6D+ lymphoid progenitors distinct from the myeloid lineage. Nat Immunol. 2019 Jul;20(7):852-864. doi: 10.1038/s41590-019-0420-3. Epub 2019 Jun 18. PMID: 31213723.
FeaturePlot(my, features = c("Ly6d", "Il7r", "Flt3", "Siglech", "Cd2", "Cd81", "Csf1r"))
# Absence of Csf1r (CD115) gene expression indicate a pre-pDC commitment
# Veryfy the overlap with the transcriptional signature 
prepDCs_signature <- c("Ly6d", "Il7r", "Flt3", "Siglech", "Cd2", "Cd81")
MSprepDCs_my <- AddModuleScore(my, features = list(prepDCs_signature), name = "prepDCs_signature")
plot_prepDCs <- FeaturePlot(MSprepDCs_my, features = "prepDCs_signature1", slot = "data", cols = rev(c("#F46D43","white","#3288BD")))



# Further identification of mature DCs with the mregDC signature: Maier B et al. A conserved dendritic-cell regulatory program limits antitumour immunity. Nature. 2020 Apr;580(7802):257-262. doi: 10.1038/s41586-020-2134-y. Epub 2020 Mar 25. Erratum in: Nature. 2020 Jun;582(7813):E17. doi: 10.1038/s41586-020-2326-5. PMID: 32269339; PMCID: PMC7787191.

mregDCs_markers <- c(
  "Slco5a1","Ly96","Ogfrl1","Kansl3","Tmem131","Tbc1d8","Map4k4","Nabp1",
  "Stat4","Stat1","Cflar","Icos","Nrp2","Ino80d","Adam23","Idh1","Mreg",
  "Dnajb2","Sp140","Sh3bp4","Agap1","Gpc1","Tnfrsf11a","Phlpp1","Nuak2",
  "Csrp1","Phlda3","Lad1","Cacna1s","Kif21b","Rgs1","Cep350","Rabgap1l",
  "Rcsd1","Pogk","Uap1","Slamf1","Pyhin1","Mfsd7b","G0s2","Dclre1c",
  "Frmd4a","Pfkfb3","Il2ra","Il15ra","Il1rn","Traf2","Tsc1","Ass1",
  "St6galnac6","Stxbp1","Traf1","Ggta1","Arl5a","Cytip","Ly75","Tank",
  "Gca","Cers6","Nostrin","Itga4","Dnajc10","Ube2l6","Arhgap1","Phf21a",
  "Traf6","Spred1","Bahd1","Pla2g4f","Vps39","B2m","Sema6d","Hdc",
  "Tmem127","Stard7","Prnp","Rassf2","Gpcpd1","Jag1","Flrt3","Dstn",
  "Ralgapa2","Srxn1","Procr","Src","Stk4","Slpi","Cd40","Zmynd8",
  "Rab22a","Ppdpf","Helz2","Arfrp1","Kdm6a","Chst7","Enox2","Renbp",
  "Emd","Gdi1","Apool","Rab9","Fabp5","Fabp4","Car13","Car2","Gyg",
  "Gnb4","Slc7a11","Rap2b","Slc33a1","Etfdh","Cd1d1","Etv3","Pmvk",
  "Slc27a3","Ints3","Npr1","S100a8","S100a9","Golph3l","Ankrd35",
  "Pias3","Nudt17","Polr3c","Rnf115","Pde4dip","Igsf3","Mab21l3",
  "Rhoc","Gnai3","Stxbp3","Vcam1","Plppr4","Fnbp1l","Synpo2","Dkk2",
  "Gbp5","Gbp2","Pnrc1","Glipr2","Nans","Galnt12","Nr4a3","Invs",
  "Akap2","Col27a1","Cdkn2a","Cdkn2b","Nsun4","St3gal3","Kdm4a",
  "Inpp5b","Mtf1","Gnl2","Ago3","Ago1","Rnf19b","Bsdc1","Marcksl1",
  "Epb41","Zdhhc18","Extl1","Clic4","Sh2d5","Iffo2","Plekhm2","Vps13d",
  "Tnfrsf1b","Kif1b","Spsb1","Tnfrsf9","Mmp23","Fam132a","Tnfrsf4",
  "Isg15","Krit1","Sri","Nub1","Tmem128","Sepsecs","Pi4k2b","Arap2",
  "Tbc1d1","Nipal1","Rufy3","Cxcl1","Anxa3","Bmp2k","Tmem150c",
  "Klhl8","Nudt9","Gbp8","Gbp9","Gbp4","Pxmp2","Adrbk2","Iscu","Ssh1",
  "Oasl1","Vsig10","Tesc","Hvcn1","Tmem120b","Sbno1","Asl","Rabgef1",
  "Gatsl2","Ncf1","Clip2","Wbscr27","Rasa4","Adap1","Fbxl18","Fscn1",
  "Cyth3","Smurf1","Lnx2","Tmem168","Tes","Ahcyl2","Strip2","Tmem140",
  "Zfp467","Herc6","Vopp1","Ndnf","Il12rb2","Cd8a","M1ap","Dok1",
  "Loxl3","Htra2","Mthfd2","Spr","Asprv1","Mxd1","Anxa4","Gfpt1",
  "Slc6a6","Ppp4r2","Wnk1","Mical3","Usp18","Apobec1","Eno2","Tapbpl",
  "Clec2i","Clec2d","Bcl2l14","Crebl2","Hebp1","Aebp2","Etnk1",
  "Caprin2","Fam60a","Amn1","Zfp524","Ceacam15","Zfp296","Relb",
  "Clptm1","Ceacam1","Sertad1","Plekhg2","Nfkbib","Kcnk6","Spint2",
  "Alkbh6","Sdhaf1","Hspb6","Nkg7","Spib","Atf5","Il4i1","Ftl1",
  "Rasip1","Ndnl2","Lrrk1","Arpin","Mex3b","Zfand6","Prkrir",
  "Stard10","Il18bp","Pgap2","Swap70","Sbf2","Rras2","Gga2","Il21r",
  "Mvp","Cd2bp2","Prr14","Fam53b","Adam8","Ptdss2","Stx11","Adgrg6",
  "Ncoa7","Marcks","Hsf2","Sgpl1","Dnajc12","Nrbf2","Ccdc6","Adora2a",
  "Gucd1","Gstt1","Icosl","Cstb","Cnn2","Efna2","Mob3a","Izumo4",
  "Gadd45b","Syn3","Dram1","Chpt1","Cdk17","Tmcc3","Plxnc1","Socs2",
  "Btg1","Tmtc2","Osbpl8","Tbc1d15","Rab21","Tmem19","Thap2","Zfc3h1",
  "Tmbim4","Rassf3","Ddit3","Ikzf4","Suox","Cd63","Rab20","Atp11a",
  "Lamp1","Polb","Nrg1","Ubxn8","Gsr","Mfhas1","Casp3","Nr2f6",
  "Fam129c","Fam32a","Large","Il15","Hook2","Nlrc5","Cpne2","Ccl22",
  "Cx3cl1","Gins3","Ndrg4","Nae1","Rrad","Nfat5","Hp","Mon1b","Cmc2",
  "Atmin","Crispld2","Fam89a","Itgb1","Arf4","Arhgap22","Ccser2",
  "Ghitm","Ddhd1","Ktn1","Ndrg2","Arhgef40","Psme2","Kpna3","Clu",
  "Epsti1","Sugt1","Tbc1d4","Tmtc4","Tmem123","Birc2","Birc3","Fut4",
  "Icam1","Zfp809","Vwa5a","Gramd1b","Ddx6","Phldb1","Scn2b","Fxyd2",
  "Dscaml1","Sik3","Zc3h12c","Tspan3","Peak1","Mpi","Sema7a","Stoml1",
  "Dennd4a","Rab8b","Rps27l","Lactb","Tpm1","Rora","Aqp9","Aldh1a2",
  "Gclc","Hmgn3","Mthfsl","Bcl2a1d","Bcl2a1a","Bcl2a1b","Rasa2",
  "Nmnat3","Pcbp4","Uba7","Bcl2a1c","Wdr48","Tmem158","Gatsl3","Myo1g",
  "Ramp3","Aftph","Peli1","Pex13","Pus10","Rel","Papolg","Ccdc88a",
  "Cpeb4","Il12b","Trim7","Psme2b","Sqstm1","Ltc4s","Tcf7","Irf1",
  "Pdlim4","Specc1","Pik3r5","Dnah2","Cxcl16","Zmynd15","Tm4sf5",
  "Eno3","Nup88","Txndc17","Ccl5","Mtmr4","Bzrap1","Mmd","Gngt2",
  "Nfe2l1","Mllt6","Arl5c","Ccr7","Smarce1","Stat5a","Stat3","Adam11",
  "Map3k14","Ace","Nup85","Gga3","Mif4gd","Mrpl38","Cyth1","Rptor",
  "Slc25a10","Net1","Klf6","Idi1","Gtpbp4","Gpr137b","Sox4","Cdkal1",
  "Serpinb1a","Serpinb6b","Serpinb9","Serpinb6a","Tubb2b","Psmg4",
  "Slc22a23","Hivep1","Cd83","Mylip","Atxn1","Kdm1b","Auh","Nfil3",
  "Cdc14b","Serinc5","Pik3r1","Pde4d","Map3k1","Gzma","Fam49a",
  "Nampt","Gdap10","Atxn7l1","Scin","Pnpla8","Fam177a","Sav1","Daam1",
  "Plek2","Zfp36l1","Ttll5","Ift43","Gtf2a1","Nrde2","Calm1",
  "Slc25a29","Traf3","Exoc3l4","Inf2","Gpr132","Itgb8","Capsl","Il7r",
  "Basp1","Ankrd33b","Laptm4b","Rnf19a","Ndrg1","Arc","Ppp1r16a",
  "Apol7c","Apol10b","C1qtnf6","Triobp","Micall1","Josd1","Gtpbp1",
  "Apobec3","Zc3h7b","Slc38a2","Adcy6","Cacnb3","Rhebl1","Tuba1a",
  "Prph","Lima1","Dip2b","Slc4a8","Grasp","Rogdi","Socs1","Snn",
  "Dnm1l","B3gnt5","Yeats2","Ccdc50","Pcyt1a","Rubcn","Sec22a",
  "Hspbap1","Pla1a","Cd80","Poglut1","Tmem39a","Arhgap31","Cd200",
  "Retnlg","Ift57","Cblb","Filip1l","Samsn1","Nrip1","Mir155hg",
  "Mx1","Rnaset2a","Mmp25","Rpl3l","Tmem8","H2-K1","H2-Eb2","H2-Q4",
  "H2-Q6","H2-Q7","H2-M2","Tmem63b","Mrpl14","Gtpbp2","Foxp4",
  "Rftn1","Stap2","Cd70","Trip10","Wash1","Arhgap28","Cdc42ep3",
  "Fbxo11","Foxn2","Crem","Cdh2","Gypc","Alpk2","Chka","Cdc42ep2",
  "Ehd1","Atl3","Cpsf7","Ms4a7","Nmrk1","Jak2","Insl6","Plgrkt",
  "Cd274","Pdcd1lg2","Ric1","Papss2","Fas","Pcgf5","Fgfbp3","March5",
  "Avpi1","Got1","Nfkb2","Gsto1","Dusp5","Vti1a","Tcf7l2",
  "Dclre1a","Atrnl1","PISD","Umad1"
)

#Featureplot only myeloid cells with mDCs signature
MS_mDCs_my <- AddModuleScore(my, features = list(mregDCs_markers), name = "mDCs_signature")
mDCs_plot <- FeaturePlot(MS_mDCs_my, features = "mDCs_signature1", slot = "data", cols = rev(c("#F46D43","white","#3288BD")))



#####################################################################################
################# Rename cell clusters with new identification ######################

DimPlot(all, group.by = "cell_type")
UMAPPlot(all)
all@meta.data

table(all$cell_type)

all$cell_type[all$cell_type == "Ly6D+ cells"] <- "pre-pDCs"
all$cell_type[all$cell_type == "CD8+ DCs"] <- "tDCs"

#Rename cluster 14 as cluster 13
all$seurat_clusters <- as.character(all$seurat_clusters)
all$RNA_snn_res.0.8 <- as.character(all$RNA_snn_res.0.8)

all$seurat_clusters[all$seurat_clusters == "14"] <- "13"
all$RNA_snn_res.0.8[all$RNA_snn_res.0.8 == "14"] <- "13"

all$seurat_clusters <- factor(
  all$seurat_clusters,
  levels = as.character(sort(as.numeric(unique(all$seurat_clusters))))
) # convert again cluster number as number keeping the numerical order 

all$RNA_snn_res.0.8 <- factor(
  all$RNA_snn_res.0.8,
  levels = as.character(sort(as.numeric(unique(all$RNA_snn_res.0.8))))
) # convert again cluster number as number keeping the numerical order 

Idents(all) <- "seurat_clusters"

table(all$cell_type)
UMAPPlot(all)
DimPlot(all, group.by = "cell_type")






###################################################################
########Exploratory analyses on different cell populations#########


#################### Myeloid cells ###############################
my <- subset(all, cell_type %in% c("Endogenous T cells", "OT-II"), invert = T)

DefaultAssay(my) <- "RNA"

my <- FindVariableFeatures(my, nfeatures = 2000)
my <- ScaleData(my)
my <- RunPCA(my, verbose = T)

DimHeatmap(my, dims = 20:30, cells = 500, balanced = TRUE)
DimPlot(my, reduction = "pca", group.by = "hash.ID")
ElbowPlot(my, 50)

variance <- (my@reductions$pca@stdev)^2
var80 <- sum(variance)*0.8
findpca <- which(cumsum(variance) >= var80)[1]
findpca 

# #Clustering
# # Clustering considering the first 25 PCs
my <- FindNeighbors(my, dims = 1:25)
my <- FindClusters(my, resolution = 0.8, verbose = FALSE)
my <- RunUMAP(my, dims = 1:25)
my <- RunTSNE(my, dims = 1:25, verbose = T)
p1 <- DimPlot(my, label = TRUE, reduction = "umap")
# DimPlot(my, label = TRUE, reduction = "tsne")

# DimPlot(my, group.by = "singleR.labels")
p2 <- DimPlot(my, group.by = "TimepointGrouped")

# DimPlot(my, group.by = "hash.ID")
p3 <- FeaturePlot(my, "Biotin-TotalSeqC")

p4 <- DimPlot(my, group.by = "Biotin_status", cols = c("grey", "red"))

p1 | p2 | p3 | p4

Idents(my) <- "cell_type"
RidgePlot(my, features = "Biotin-TotalSeqC")
VlnPlot(my, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = my@misc$Biotin_threshold, linetype = "dashed", color = "red") 

Idents(my) <- "seurat_clusters"
RidgePlot(my, features = "Biotin-TotalSeqC")
VlnPlot(my, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = my@misc$Biotin_threshold, linetype = "dashed", color = "red") 


FeaturePlot(my, "Biotin-TotalSeqC")

Idents(my) <- "seurat_clusters"

VlnPlot(my, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = my@misc$Biotin_threshold, linetype = "dashed", color = "red") 

UMAPPlot(my, label = T)
DimPlot(my, group.by = "cell_type")



#################### Monocytes ########################

mono <- subset(my, cell_type == "Monocytes")

mono <- FindVariableFeatures(mono, nfeatures = 2000)
mono <- ScaleData(mono)
mono <- RunPCA(mono, verbose = T)

ElbowPlot(mono, 50)
DimHeatmap(mono, dims = 20:30, cells = 500, balanced = TRUE)
DimPlot(mono, reduction = "pca", group.by = "hash.ID")

variance <- (mono@reductions$pca@stdev)^2
var80 <- sum(variance)*0.8
findpca <- which(cumsum(variance) >= var80)[1]
findpca 

# #Clustering
# # Clustering considering the first 30 PCs
mono <- FindNeighbors(mono, dims = 1:30)
mono <- FindClusters(mono, resolution = 0.8, verbose = FALSE)
mono <- RunUMAP(mono, dims = 1:30)
mono <- RunTSNE(mono, dims = 1:30, verbose = T)

UMAPPlot(mono, label = T)
FeaturePlot(mono, features = "Biotin-TotalSeqC")


VlnPlot(mono, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = my@misc$Biotin_threshold, linetype = "dashed", color = "red") 




###################################################################
########################### DCs only ##############################

dc <- subset(my, cell_type == "Monocytes", invert = T)

dc <- FindVariableFeatures(dc, nfeatures = 2000)
dc <- ScaleData(dc)
dc <- RunPCA(dc, verbose = T)

ElbowPlot(dc, 50)
DimHeatmap(dc, dims = 30:40, cells = 500, balanced = TRUE)
DimPlot(dc, reduction = "pca", group.by = "hash.ID")


variance <- (dc@reductions$pca@stdev)^2
var80 <- sum(variance)*0.8
findpca <- which(cumsum(variance) >= var80)[1]
findpca 

# #Clustering
# # Clustering considering the first 30 PCs
dc <- FindNeighbors(dc, dims = 1:30)
dc <- FindClusters(dc, resolution = 0.8, verbose = FALSE)
dc <- RunUMAP(dc, dims = 1:30)
dc <- RunTSNE(dc, dims = 1:30, verbose = T)
#p1 <- DimPlot(my, label = TRUE, reduction = "umap")
TSNEPlot(dc)
UMAPPlot(dc)

DimPlot(dc, group.by = "Biotin_status", cols = c("grey", "blue"))
RidgePlot(dc, features = "Biotin-TotalSeqC")

VlnPlot(dc, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = my@misc$Biotin_threshold, linetype = "dashed", color = "red") 



## dplyr to calculate the median Biotin levels per group

biotin_data <- FetchData(dc, vars = "Biotin-TotalSeqC")
biotin_data$cluster <- dc$seurat_clusters

median_biotin <- biotin_data %>%
  group_by(cluster) %>%
  summarise(median_expression = median(`adt_Biotin-TotalSeqC`, na.rm = TRUE))

median_biotin

DimPlot(dc, group.by = "cell_type")



####################  mDCs  ##########################

mdcs <- subset(dc, cell_type == "mDCs")

mdcs <- FindVariableFeatures(mdcs, nfeatures = 2000)
mdcs <- ScaleData(mdcs)
mdcs <- RunPCA(mdcs, verbose = T)

ElbowPlot(mdcs, 50)
DimHeatmap(mdcs, dims = 10:20, cells = 500, balanced = TRUE)
DimPlot(mdcs, reduction = "pca", group.by = "hash.ID")


variance <- (mdcs@reductions$pca@stdev)^2
var80 <- sum(variance)*0.8
findpca <- which(cumsum(variance) >= var80)[1]
findpca 

# #Clustering
# # Clustering considering the first 30 PCs
mdcs <- FindNeighbors(mdcs, dims = 1:30)
mdcs <- FindClusters(mdcs, resolution = 0.8, verbose = FALSE)
mdcs <- RunUMAP(mdcs, dims = 1:30)
#mdcs <- RunTSNE(mdcs, dims = 1:30, verbose = T)

#TSNEPlot(mdcs, label=T)
UMAPPlot(mdcs, label = T)

FeaturePlot(mdcs, features = "Biotin-TotalSeqC")
RidgePlot(mdcs, features = "Biotin-TotalSeqC")

VlnPlot(mdcs, features = "Biotin-TotalSeqC") +
  stat_summary(fun.y = median, geom='point', size = 20, colour = "red", shape = 95) +
  theme(legend.position = 'none') +
  geom_hline(yintercept = my@misc$Biotin_threshold, linetype = "dashed", color = "red") 


DimPlot(mdcs, group.by = "TimepointGrouped")


#############################################################
###################### Save objects #########################

saveRDS(my, file = "Myeloid_clusters")
saveRDS(mono, file = "Monocyte_clusters")
saveRDS(dc, file = "DC_clusters")
saveRDS(mdcs, file = "mDCs_clusters")


