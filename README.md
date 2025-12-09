# computational_cancer_Hackathon_2025
HNSCC transcriptome analysis and logica building


## R packages

Bioconductor
TCGAbiolinks
DESeq2
ggplot
## Data
https://portal.gdc.cancer.gov/projects/TCGA-HNSC

## Downloading
#Doownload the data from GDC data base
query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor")
)
GDCdownload(query)


raw_count_1 <- GDCprepare(query,summarizedExperiment = TRUE)


Clinical data

query_cli <- GDCquery(project = "TCGA-HNSC",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement"
                  )
GDCdownload(query_cli)
## Data filtering
- 
