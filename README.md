# G4 propensity in cf-RNA and Type 1 Interferon EigenGene in scRNAseq - Preeclampsia

## Methodology
Identification of G-quadruplex (G4) structures in the cell-free RNA transcripomic data in the genes having predictive power to differentiate pre-eclampsia from normal pregnancy [1,2]. We identified 3880 (15%) G4 regions across 674 genes that were found be predictors of the pregnancy (estimation for the G4 propensity is ongoing for ~9000 genes which were found to be correlated in [2]), shown in figure 2A,B.

Additionally, single cell rnaseq data from normal, early-onset preeclampsia and late-onset preeclampsia, we performed Seurat scRNA analysis, and Type 1 interferon pathway-related genes were found to be enriched across the groups (Figure 2C). 

Also, we performed cell type-specific differential expression analysis, and the top 500 genes were used for building machine learning models (random forests and XGBoost) with 5-fold cross-validation highlighting probable differential pathways across cell-type gene regulatory pathway enrichments.
