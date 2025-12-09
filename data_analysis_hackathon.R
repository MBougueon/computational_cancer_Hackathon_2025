library(dplyr)
library(ggplot2)
library(TCGAbiolinks)
library(DESeq2)
library(readxl)

#setwd("/home/mbougueon/UCL_post_doc_2426/Hackathon/")

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

# PCA plot
pca_data <- plotPCA(vsd2, intgroup = "HPV_status", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = HPV_status)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw()

dir.create("results", showWarnings = FALSE)
ggsave("results/PCA_plot.png", p, width = 6, height = 5, dpi = 300)

# ---------------------------
# Differential expression: HPV+ vs HPV-
# ---------------------------

# Run DESeq2
dds2 <- DESeq(dds2)

# Get results for HPV positive vs negative
res <- results(dds2, contrast = c("HPV_status", "positive", "negative"))

# Order by p-value (smallest first)
resOrdered <- res[order(res$pvalue), ]

# Remove rows with NA p-values (just to be safe)
resOrdered <- resOrdered[!is.na(resOrdered$pvalue), ]

# Save results for later / for report
saveRDS(resOrdered, file = "results/DE_results.rds")
write.csv(as.data.frame(resOrdered), file = "results/DE_results.csv")

# ---------------------------
# Volcano plot - NOTE FROM NAGWA: The volcano plot shows which genes differ between HPV-positive and HPV-negative HNSC tumors.
# The x-axis shows how much each gene changes, and the y-axis shows the significance.
#We found a large number of significantly upregulated genes in HPV-positive tumors, as shown by the blue points on the right side of the volcano.
#This suggests that HPV infection drives a distinct gene expression profile.
# ---------------------------

res_df <- as.data.frame(resOrdered)
res_df$logp <- -log10(res_df$pvalue)
res_df$significant <- res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) > 1

v <- ggplot(res_df, aes(x = log2FoldChange, y = logp, color = significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  xlab("Log2 fold change (HPV+ vs HPV-)") +
  ylab("-log10(p-value)") +
  ggtitle("Differential expression: HPV+ vs HPV-")

ggsave("results/Volcano_plot.png", v, width = 7, height = 5, dpi = 300)