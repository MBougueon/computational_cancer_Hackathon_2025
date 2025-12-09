library(dplyr)
library(ggplot2)
library(TCGAbiolinks)
library(DESeq2)


setwd("/home/mbougueon/UCL_post_doc_2426/Hackathon/")

query <- GDCquery( project = "TCGA-HNSC", data.category = "Transcriptome Profiling", 
                   data.type = "Gene Expression Quantification", 
                   workflow.type = "STAR - Counts", sample.type = c("Primary Tumor") )
GDCdownload(query)

raw_count_1 <- GDCprepare(query,summarizedExperiment = TRUE)

# counts and metadata
count_matrix2 <- assay(raw_count_1)
sample_metadata2 <- colData(raw_count_1)

hpv <- read_excel("hpvexp.xlsx", sheet = 'Table S1M')[,c('Patient Barcodes', 'HPV_status')]
hpv$id <- hpv$`Patient Barcodes`
colData(raw_count_1)$HPV_status <- hpv$HPV_status[match(colData(raw_count_1)$submitter_id, hpv$id)]

to_keep <- raw_count_1@colData$submitter_id %in% hpv$id
data_filtered2 <- raw_count_1[, to_keep]

head(data_filtered2@colData)

# counts and metadata
count_matrix2 <- assay(data_filtered2)
sample_metadata2 <- colData(data_filtered2)
sample_metadata2$HPV_status <- as.factor(sample_metadata2$HPV_status)

# plate id for batch effect correction
sample_metadata2$plate <- substr(sample_metadata2$barcode, 22, 25)
sample_metadata2$plate <- as.factor(sample_metadata2$plate)

# Construct a DESeqDataSet object
dds2 <- DESeqDataSetFromMatrix(countData=count_matrix2,
                               colData=sample_metadata2,
                               design=~plate + HPV_status)

# Filtering: removing rows with low gene counts
keep2 <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep2,]

vsd2 <- vst(dds2, blind = FALSE) 
plotPCA(vsd2, intgroup = "HPV_status")
