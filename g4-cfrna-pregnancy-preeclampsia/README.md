# Identification of G4 Propensity in the predictive or DE genes cf-RNA during normal/disease status of pregnancy

## Methodology
G4 propensity for a given gene was predicted using the G4Hunter algorithm which was developed for chromosome level estimation of G4 regions in DNA. We repurpose the algorithm for a gene and its all transcripts harboring exons only using Bioconductor packages from mapping gene symbol to ENTREZID to Transcript features to their exon level features. The exon features were merged to transcript features andits correspinding sequence. Using this G4Hunter algorithm was applied and all the predicted G4 regions and G4Score (G4 propensity) was estimated for all the regions. As suggested by the authors of G4Hunter, G4Score of greater than or equals to 1 were annotated as G4 for the genes differentially expressed or predictors of gestational age based on cell-free RNA (cfRNA)
       
##### Cell-Free RNA data
