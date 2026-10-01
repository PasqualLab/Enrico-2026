# Enrico-2026

# DC–T Cell CITE-seq Analysis Pipeline

This repository contains the R scripts used for the processing and analysis of CITE-seq data generated to investigate dendritic cell (DC)–CD4⁺ T cell interactions during T cell priming.

The dataset combines **gene expression (GEX)**, **antibody-derived tags (ADT)**, **sample hashtags (HTO)**, and **Biotin-TotalSeqC labeling** as a readout of recorded cellular interactions. The analysis was performed primarily using the `Seurat` framework.

The scripts are intended to be executed sequentially:

```text
1_Cleaning pipeline Aggr.R
        ↓
2_Timepoint Aggr.R
        ↓
3_Cell Annotation Aggr.R
        ↓
4_Biotin Threshold setting Aggr.R
        ↓
5_Objects reclustering by cell type Aggr.R
        ↓
6_mDCs_analysis Aggr.R
```

## 1. Data cleaning and quality control

**Script:** `1_Cleaning pipeline Aggr.R`

This script imports the Cell Ranger feature-barcode matrix and generates the initial Seurat object containing separate **RNA, HTO, and ADT assays**. HTOs correspond to the first 12 antibody-capture features, whereas the remaining antibody-derived features are stored in the ADT assay.

Initial quality control removes cells with:

- fewer than 500 detected genes;
- more than 5,500 detected genes;
- more than 30,000 RNA counts;
- mitochondrial RNA content ≥2%.

The RNA data are subsequently normalized, variable features are identified, and PCA is performed. The number of PCs used for downstream analysis is evaluated using both the elbow plot and cumulative explained variance.

Potential doublets are identified using **DoubletFinder**. The expected doublet rate is adjusted for the estimated proportion of homotypic doublets, and cells classified as doublets are removed.

Finally, the HTO assay is CLR-normalized and **HTODemux** is used to assign cells to individual samples. Cells lacking a unique hashtag assignment are excluded.

The resulting cleaned object is saved as:

```text
Cleaned_DS
```

The initial preprocessing and filtering strategy is implemented in the first script.  DoubletFinder-based filtering and HTO demultiplexing are subsequently performed before saving the cleaned dataset. 

## 2. Dimensionality reduction and time-point assignment

**Script:** `2_Timepoint Aggr.R`

The cleaned dataset is reprocessed for downstream analysis using RNA-based dimensionality reduction and clustering. Variable features are identified and PCA is followed by graph-based clustering at a resolution of `0.8`. Both UMAP and t-SNE representations are generated.

The ADT assay is independently normalized using CLR normalization.

Sample hashtags are then used to assign cells to their experimental time point:

```text
Hashtag 01–06 → 24 h
Hashtag 07–12 → 72 h
```

This information is stored in the metadata column:

```r
TimepointGrouped
```

The resulting object is saved as:

```text
Cleaned_Timepoint_DS
```

The clustering and dimensionality-reduction workflow uses 30 PCs, while hashtag identities provide the experimental time-point annotation. 

## 3. Cell-type annotation

**Script:** `3_Cell Annotation Aggr.R`

Cell identities are determined by combining **reference-based annotation**, RNA expression, and surface-protein expression.

An initial automated annotation is obtained using **SingleR** with the **ImmGen** reference dataset. The Seurat object is converted to a `SingleCellExperiment` object and SingleR predictions are stored in the Seurat metadata.

Cell identities are then manually refined using canonical transcriptional and surface markers. Marker sets are inspected for major populations including:

- conventional dendritic cells;
- cDC1;
- cDC2;
- mature DCs;
- monocytes;
- endogenous T cells;
- transferred OT-II cells.

Both RNA and ADT information are used during this annotation step.

Clusters are subsequently assigned to the major cell populations, including monocytes, cDC1s, cDC2s, mature DCs (mDCs), OT-II cells, endogenous T cells and additional DC populations. B cells and residual doublets are excluded from subsequent analyses. 

The normalized `Biotin-TotalSeqC` signal is additionally stored in the metadata as:

```r
Biotin_data
```

The annotated dataset is saved as:

```text
All_clusters
```

## 4. Definition of Biotin-positive interacting cells

**Script:** `4_Biotin Threshold setting Aggr.R`

Biotin labeling is used as the readout for recorded cellular interactions.

To estimate background Biotin signal, endogenous T cells are used as a negative-reference population. Cells showing CD45.1 signal are first excluded from this population to minimize contamination by transferred OT-II cells.

The 95th percentile of the Biotin distribution in the resulting negative-control population is calculated as an empirical estimate of background:

```text
95th percentile ≈ 1.095
```

For the downstream analysis, the operational threshold is set to:

```r
threshold <- 1
```

Cells are consequently classified as:

```text
Biotin-TotalSeqC > 1 → positive
Biotin-TotalSeqC ≤ 1 → negative
```

and the resulting classification is stored as:

```r
Biotin_status
```

The rationale and calculation of the background threshold are documented directly in the script. The threshold is then applied to the complete dataset to generate the `Biotin_status` classification used downstream.

## 5. Cell-type-specific reclustering

**Script:** `5_Objects reclustering by cell type Aggr.R`

To resolve heterogeneity within the broad cell populations identified above, selected populations are extracted and independently reclustered.

The analysis focuses primarily on the myeloid compartment. Additional annotation of specific DC populations is supported by published transcriptional signatures, including signatures associated with **transitional DCs (tDCs)**, **Ly6D⁺ pre-pDCs**, and **mature regulatory DCs (mregDCs)**.

Separate objects are generated for:

```text
Myeloid cells
├── Monocytes
└── Dendritic cells
    └── Mature DCs
```

Each population is independently processed through variable-feature selection, scaling, PCA, nearest-neighbor graph construction, clustering, and UMAP/t-SNE dimensionality reduction.

The mature DC population is specifically extracted from the DC compartment and reclustered using 30 PCs at a clustering resolution of `0.8`.

The following objects are saved for downstream analyses:

```text
Myeloid_clusters
Monocyte_clusters
DC_clusters
mDCs_clusters
```



## 6. Mature dendritic cell analysis

**Script:** `6_mDCs_analysis Aggr.R`

The final script performs a focused analysis of the reclustered mature DC compartment.

### Cluster characterization

Mature DC clusters are characterized using both transcriptional and surface-marker information. Expression of DC-associated markers such as XCR1, CD8α, SIRPα, CD11b, CD301b and CD103 is visualized across the UMAP.

ADT expression is additionally summarized using dot plots to facilitate cluster identification.

### Differential gene expression

Differentially expressed genes among all mature DC clusters are identified using:

```r
FindAllMarkers(mdc)
```

and exported as:

```text
mDCs_AllvsAll_DEG.csv
```



### Comparison between 24 h and 72 h

The mature DC dataset is divided according to experimental time point:

```r
mdc24h <- subset(mdc, TimepointGrouped == "24h")
mdc72h <- subset(mdc, TimepointGrouped == "72h")
```

Cluster composition is evaluated separately at 24 h and 72 h, and the relative representation of each cluster is compared between the two time points using the **24 h / 72 h ratio**.

### Analysis of interacting DCs

Biotin labeling is then used to investigate which mature DC populations participate in recorded DC–T cell interactions.

Two complementary measurements are calculated separately at 24 h and 72 h:

1. **Biotin positivity within each cluster** — the percentage of cells belonging to a given DC cluster that are Biotin⁺.
2. **Composition of the Biotin⁺ population** — the percentage contribution of each DC cluster to the total pool of Biotin⁺ mature DCs.

The first measurement estimates the frequency of recorded interactions within each DC population, whereas the second describes the cellular composition of the overall interacting DC compartment.

These analyses are performed independently at 24 h and 72 h and visualized using bar plots and pie charts. 

## Main R packages

The workflow relies primarily on:

```r
Seurat
SingleR
SingleCellExperiment
celldex
DoubletFinder
tidyverse
dplyr
ggplot2
patchwork
cowplot
RColorBrewer
scales
```

## Pipeline summary

Overall, the workflow progresses from **quality control and sample demultiplexing**, through **cell-type annotation and interaction classification**, to **cell-type-specific reclustering and characterization of interacting mature dendritic cell populations**.

```text
Cell Ranger feature-barcode matrix
              │
              ▼
       Quality control
              │
              ▼
       Doublet removal
              │
              ▼
       HTO demultiplexing
              │
              ▼
  24 h / 72 h assignment
              │
              ▼
     Cell-type annotation
              │
              ▼
   Biotin+/Biotin− assignment
              │
              ▼
 Myeloid/DC/mDC reclustering
              │
              ▼
   Mature DC characterization
              │
              ▼
24 h vs 72 h interaction analysis
```
