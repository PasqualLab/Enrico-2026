# Enrico-2026

# DC–T Cell Interaction CITE-seq Analysis

This repository contains the R scripts used for the processing and analysis of CITE-seq data generated to investigate dendritic cell (DC)–CD4⁺ T cell interactions during T cell priming.

The dataset combines gene expression (GEX), antibody-derived tags (ADT), sample hashtags (HTO), and Biotin-TotalSeqC labeling as a readout of recorded cellular interactions.

## Analysis Workflow

The analysis is organized as a sequential pipeline. Each script generates the Seurat object used by the following analysis step.

```text
1_Cleaning pipeline Aggr.R
        |
        v
2_Timepoint Aggr.R
        |
        v
3_Cell Annotation Aggr.R
        |
        v
4_Biotin Threshold setting Aggr.R
        |
        v
5_Objects reclustering by cell type Aggr.R
        |
        v
6_mDCs_analysis Aggr.R
```

## Scripts

| Step | Script | Description | Main output |
| --- | --- | --- | --- |
| 1 | `1_Cleaning pipeline Aggr.R` | Initial Seurat object construction, QC filtering, doublet removal, and HTO demultiplexing | `Cleaned_DS` |
| 2 | `2_Timepoint Aggr.R` | RNA clustering, dimensionality reduction, ADT normalization, and assignment of cells to the 24 h or 72 h experimental time points | `Cleaned_Timepoint_DS` |
| 3 | `3_Cell Annotation Aggr.R` | SingleR/ImmGen-assisted annotation followed by manual cell-type annotation using RNA and ADT markers | `All_clusters` |
| 4 | `4_Biotin Threshold setting Aggr.R` | Definition of the Biotin background threshold and classification of cells as Biotin-positive or Biotin-negative | updated `All_clusters` |
| 5 | `5_Objects reclustering by cell type Aggr.R` | Cell-type-specific subsetting and reclustering of myeloid cells, monocytes, DCs, and mature DCs | `Myeloid_clusters`, `Monocyte_clusters`, `DC_clusters`, `mDCs_clusters` |
| 6 | `6_mDCs_analysis Aggr.R` | Detailed analysis of mature DC clusters, differential gene expression, time-point composition, and Biotin-positive interacting DC populations | plots and `mDCs_AllvsAll_DEG.csv` |

## 1. Quality Control and Object Construction

The Cell Ranger feature-barcode matrix is imported using `Read10X`. Gene-expression counts are used to construct the main Seurat object, while antibody-capture features are separated into HTO and ADT assays.

Cells are retained using the following QC criteria:

| QC metric | Threshold |
| --- | --- |
| Detected genes | `500 < nFeature_RNA < 5500` |
| RNA counts | `nCount_RNA < 30000` |
| Mitochondrial RNA | `percent.mt < 2%` |

RNA data are normalized and subjected to PCA and initial clustering. Potential doublets are identified using **DoubletFinder**, with the expected doublet number corrected for the estimated homotypic doublet fraction. 

The HTO assay is subsequently CLR-normalized and `HTODemux` is used for sample demultiplexing. Cells without a unique hashtag assignment are removed.

## 2. Time-point Assignment

The cleaned dataset is reclustered using the first 30 PCs at a resolution of `0.8`, and UMAP and t-SNE representations are generated. The ADT assay is independently normalized using CLR normalization.

Hashtag identity is used to assign the experimental time point:

| Hashtags | Time point |
| --- | --- |
| Hashtag 01–06 | 24 h |
| Hashtag 07–12 | 72 h |

The resulting annotation is stored in the `TimepointGrouped` metadata column.

## 3. Cell-type Annotation

Cell identities are assigned using a combination of automated and manual annotation.

Automated annotation is performed with **SingleR** using the **ImmGen** reference dataset.

Predicted identities are subsequently refined using canonical RNA and ADT markers for:

- cDC1 and cDC2
- mature DCs
- monocytes
- OT-II cells
- endogenous T cells
- additional lymphoid and myeloid populations

Both gene-expression and surface-protein information are used for manual cluster identification.

The normalized `Biotin-TotalSeqC` signal is also added to the Seurat metadata as `Biotin_data`. B cells and residual doublets are removed before downstream analysis.

## 4. Biotin Threshold

Biotin-TotalSeqC is used as the readout of recorded cellular interactions.

Endogenous T cells are used as a negative-reference population to estimate background Biotin labeling. CD45.1-positive cells are first removed from this population to exclude contaminating transferred OT-II cells.

The 95th percentile of the Biotin distribution in this negative population is approximately `1.095`. An operational threshold of `1` is used for downstream analyses.

Cells are classified as:

```text
Biotin-TotalSeqC > 1  -> Biotin-positive
Biotin-TotalSeqC <= 1 -> Biotin-negative
```

The classification is stored in the `Biotin_status` metadata column.

## 5. Cell-type-specific Reclustering

The myeloid compartment is extracted from the complete dataset and further characterized using RNA expression, ADT expression, and published transcriptional signatures.

Additional signatures are used to support the identification of transitional DCs, Ly6D⁺ pre-pDCs, and mature regulatory DC populations.

Individual populations are subsequently extracted and independently reclustered:

```text
All cells
   |
   +-- Myeloid cells
          |
          +-- Monocytes
          |
          +-- Dendritic cells
                  |
                  +-- Mature DCs
```

For each population, variable-feature selection, scaling, PCA, graph construction, clustering, and dimensionality reduction are repeated.

The following objects are generated:

```text
Myeloid_clusters
Monocyte_clusters
DC_clusters
mDCs_clusters
```

The mature DC population is independently reclustered using the first 30 PCs.

## 6. Mature DC Analysis

The `mDCs_clusters` object is used for the final mature DC analysis.

Clusters are characterized using RNA and ADT expression, including XCR1, CD8α, SIRPα, CD11b, CD301b, and CD103.

Differentially expressed genes between mature DC clusters are identified using `FindAllMarkers`, with the results exported as:

```text
mDCs_AllvsAll_DEG.csv
```

The dataset is then separated into 24 h and 72 h populations to compare the relative abundance of individual mature DC clusters between the two time points.

Finally, Biotin-positive cells are analyzed to quantify:

- the percentage of Biotin⁺ cells within each mature DC cluster;
- the contribution of each cluster to the total Biotin⁺ mature DC population;
- differences in these distributions between 24 h and 72 h.

These analyses provide the basis for comparing the DC populations participating in recorded DC–T cell interactions at different stages of T cell priming.

## Main R Packages

The analysis uses the following R packages:

```text
Seurat
DoubletFinder
SingleR
SingleCellExperiment
celldex
tidyverse
dplyr
ggplot2
patchwork
cowplot
RColorBrewer
scales
```

## Main Seurat Objects

| Object | Description |
| --- | --- |
| `Cleaned_DS` | QC-filtered and hashtag-demultiplexed dataset |
| `Cleaned_Timepoint_DS` | Dataset with 24 h / 72 h annotation |
| `All_clusters` | Cell-type-annotated dataset with Biotin classification |
| `Myeloid_clusters` | Reclustering of the myeloid compartment |
| `Monocyte_clusters` | Reclustering of monocytes |
| `DC_clusters` | Reclustering of dendritic cells |
| `mDCs_clusters` | Reclustering of mature dendritic cells |
