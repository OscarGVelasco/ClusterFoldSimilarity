# ClusterFoldSimilarity
Calculate cluster similarity between clusters from different single cell datasets/batches/samples.


Installation
-----------------------------

1. The package binaries are available for download on github:
https://github.com/OscarGVelasco/ClusterFoldSimilarity/blob/main/ClusterFoldSimilarity.tar.gz


Example using pancreatic scRNA-Seq
-----------------------------
We will use **a set of single-cell RNA-Seq transcriptomic data from human pancreatic samples** included on the R Bioconductor package **scRNAseq** (Risso D, C. M. (2021). scRNAseq: Collection of Public Single-Cell RNA-Seq Datasets. R package version 2.8.0.)*
1) GSE84133 *(Baron et al., 2016)*, 2) GSE86469 *(Lawlor et al., 2017)*, 3) GSE81608 *(Xin et al.,2016)*, 4) ArrayExpress: E-MTAB-5061 *(Segerstolpe et al., 2016)*.

The four selected datasets contain metadata specifing the cell-type of each barcoded cell on the dataset.

``` r
library(Seurat)
library(scRNAseq)

pancreas_baron <- scRNAseq::BaronPancreasData(which = "human",ensembl = T)
table(colData(pancreas_baron)$label)
table(colData(pancreas_baron)$donor)

pancreas_lawlor <- scRNAseq::LawlorPancreasData()
table(colData(pancreas_lawlor)[,"cell type"])

pancreas_xin <- scRNAseq::XinPancreasData(ensembl = T)
assayNames(pancreas_xin) <- "counts"
table(colData(pancreas_xin)$cell.type)

pancreas_segers <- scRNAseq::SegerstolpePancreasData(ensembl = T)
assayNames(pancreas_segers)
```

# Create sce object

``` r
load("pancreas.sc.list.RData")
simi.table <- ClusterFoldSimilarity::cluster_fold_similarity(sce_list = pancreas.sc.list)
```

<object data="README_files/Figure_S1_supp.pdf" type="application/pdf" width="100%"> 
</object>
