#This is a tutorial to analyse spatial transcriptomics data from Geomx 
#using standR pipeline.
#majorly derived from the following tutorial pipeline.
#https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html
#Created for Statistical Genomics workshop 2025 conducted in NCBS.
#Author: Shrisruti Sriraman

#Install the packages once and load the libraries for subsequent runs.
#install.packages("tidyverse")
#install.packages("BH")
#install.packages("RSpectra")
#install.packages("ggpubr")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("BiocParallel")
#BiocManager::install("edgeR")
#BiocManager::install("limma")
#BiocManager::install("GenomicAlignments")
#BiocManager::install("EDASeq")
#BiocManager::install("SpatialExperiment")
#BiocManager::install("standR")
#BiocManager::install("vissE")
#BiocManager::install("scater")

#required libraries
library(standR)
library(SpatialExperiment)
library(edgeR)
library(limma)
library(tidyverse)
library(dplyr)
library(readxl)
library(Matrix)
library(scater)
library(vissE)
library(ggplot2)
library(ggrepel)
library(ggalluvial)
library(ggpubr)
library(GSEABase)
library(msigdb)

################################################################################
### I. Loading data          
################################################################################
#1. Load count matrix file (gene read counts)
count_file <- read_xlsx(
  path = "data/Export3_BiologicalProbeQC.xlsx", 
  sheet = "TargetCountMatrix") %>% as.data.frame()
head(count_file)[ ,1:5]

#2. Load the sample annotation file (sample/patient metadata file)
sample_annofile <- read_xlsx(
  path="data/Export3_BiologicalProbeQC.xlsx",
  sheet = "SegmentProperties") %>% as.data.frame()
head(sample_annofile)[, 1:5]

#use table() to understand the sample statistics.
table(sample_annofile$Type)
table(sample_annofile$`CD3+`)
table(sample_annofile$CD3)
table(sample_annofile$`CD20+`)
table(sample_annofile$CD20)
table(sample_annofile$`SMA+`)
table(sample_annofile$CD11c)

#3. Load the feature annotation file (gene annotation file)
feature_annofile <- read_xlsx(
  path="data/Export3_BiologicalProbeQC.xlsx",
  sheet = "TargetProperties") %>% as.data.frame()
head(feature_annofile)[, 1:5]

################################################################################
### II. Create Spatial Experiment object using readGeoMx()          
################################################################################
spe <- readGeoMx(count_file, sample_annofile, feature_annofile)
spe

#Let us understand the data stored in spe
dim(spe)
metadata(spe) |> names()
assayNames(spe)
colData(spe)
rowData(spe)
################################################################################
### III. Quality Control          
################################################################################
#A. Gene level QC
#default settings: min_count = 5 and sample_fraction = 0.9
#Question: What do you think is happening here?
spe <- addPerROIQC(spe, rm_genes = TRUE)
dim(spe)
metadata(spe) |> names()

#Question: If we change the sample_fraction to 0.5, 
#will we expect more or less genes to be removed?
plotGeneQC(spe, ordannots = "regions", col = regions, point_size = 2)

#B. ROI level QC
#The two metrics considered here are low cell count and library sizes.
#Question:
#1. Which annotation file has information about the cell count?
#2. In DSP, how are the cells counted?
#3. What does library size tell us about the sample quality?

plotROIQC(spe, x_threshold = 150, color = SlideName)
qc <- colData(spe)$AOINucleiCount > 150
table(qc)

spe <- spe[ ,qc]
dim(spe)

plotROIQC(spe,  x_axis = "AOISurfaceArea", x_lab = "AreaSize",
          y_axis = "lib_size", y_lab = "Library size", col = SlideName)
################################################################################
### IV. Dimension reduction          
################################################################################
set.seed(100)
spe <- scater::runPCA(spe)
pca_results <- reducedDim(spe, "PCA")
drawPCA(spe, precomputed = pca_results, col = Type)

#Question: Are the samples from the 5 different slides also different?
# Replace ? with the appropriate colname to identify the slide based difference.
drawPCA(spe, precomputed = pca_results, col = SlideName)

#PCA paired plots
plotPairPCA(spe, col = Type, precomputed = pca_results, n_dimension = 4)

#UMAP
#Question: Do these graph inform us more about the samples?
set.seed(100)
spe <- scater::runUMAP(spe, dimred = "PCA")
plotDR(spe, dimred = "UMAP", col = Type)

################################################################################
### V. Normalization          
################################################################################
#standR has different normalization options - TMM, RPKM, TPM, CPM, UQ and sizefactor
spe_tmm <- geomxNorm(spe, method = "TMM")

#Question: Where is the normalized factor stored in the spe object?
plotRLExpr(spe_tmm, assay = 2, color = SlideName) + ggtitle("TMM")
################################################################################
### VI. Batch correction          
################################################################################
#Question: Has the normalization reduced the batch effect? 
#Test it out using the following codes.
set.seed(100)

spe_tmm <- scater::runPCA(spe_tmm)
pca_results_tmm <- reducedDim(spe_tmm, "PCA")
plotPairPCA(spe_tmm, precomputed = pca_results_tmm, color = Type)

#Question: What did you observe?

#Batch correction using RUV
#Obtain the negative control genes from the data.
#Here top_n refers to the list of least variable genes 
spe <- findNCGs(spe, batch_name = "SlideName", top_n = 300)
metadata(spe) |> names()

for(i in seq(3)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "Type", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = Type, 
                    title = paste0("k = ", i)))
}

#Question: Which k value separates better?
spe_ruv <- geomxBatchCorrection(spe, factors = "Type", 
                                NCGs = metadata(spe)$NCGs, k = 4)

set.seed(100)
spe_ruv <- scater::runPCA(spe_ruv)
pca_results_ruv <- reducedDim(spe_ruv, "PCA")
plotPairPCA(spe_ruv, precomputed = pca_results_ruv, 
            color = SlideName, title = "RUV4, k = 2", n_dimension = 4)

#Batch correction using limma
spe_lrb <- geomxBatchCorrection(spe,
                                batch = colData(spe)$SlideName, method = "Limma",
                                design = model.matrix(~Type, data = colData(spe)))
plotPairPCA(spe_lrb, assay = 2, color = Type, title = "Limma removeBatch")

#Evaluate the two methods and choose the best one
spe_list <- list(spe, spe_ruv, spe_lrb)
plotClusterEvalStats(spe_list = spe_list,
                     bio_feature_name = "Type",
                     batch_feature_name = "SlideName",
                     data_names = c("Raw","RUV4","Limma"))

#RLE plots
plotRLExpr(spe_ruv, assay = 2, color = SlideName) + ggtitle("RUV4")
plotRLExpr(spe_lrb, assay = 2, color = SlideName) + ggtitle("Limma removeBatch")

#Which one do we select?
################################################################################
### VII. Differential expression analysis          
################################################################################
#DEG analysis using limma-voom pipeline

#Include the weights generated from the batch correction as covariates.
colData(spe_ruv)[,seq(ncol(colData(spe_ruv))-1, ncol(colData(spe_ruv)))] |>
  head()

#Create DGE object
dge <- SE2DGEList(spe_ruv)

#Compute the design matrix
design <- model.matrix(~0 + Type + ruv_W1 + ruv_W2 , data = colData(spe_ruv))
colnames(design)
colnames(design) <- gsub("^Type","",colnames(design))
colnames(design) <- gsub(" ","_",colnames(design))
colnames(design)

#Contrast matrix
contr.matrix <- makeContrasts(
  BvT = B_cell_zone - T_cell_zone,
  levels = colnames(design))

#Question: Create a contrast matrix between Germinal center and Medulla 
#and save it to \contr.matrix_1 

#Filter out the low expressed genes
keep <- filterByExpr(dge, design)
table(keep)

#Question: How many genes are removed?

rownames(dge)[!keep]
dge_all <- dge[keep, ]

#BCV check
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)
plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))
bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)
#Question: What can you comment about the high BCV genes?

highbcv <- bcv_df$BCV > 0.8
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)

#Limma-voom pipeline for DEG analysis
v <- voom(dge_all, design, plot = TRUE) 
fit <- lmFit(v)
fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit_contrast, robust = TRUE)
results_efit<- decideTests(efit, p.value = 0.05)
summary_efit <- summary(results_efit)
summary_efit

#Question: Obtain the top 20 up and down regulated DEGs.
#Use log2FC > 2 and  < -2 with p < 0.05 to filter the DEGs. 
#Also, check what happens if we decrease the FC cutoff.

#Visualization
de_results_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf)
de_genes_toptable_BvT <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

#Question: Create plots for different FC cutoff.
de_results_BvT %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                     ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
  ggplot(aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = de_genes_toptable_BvT %>% 
                    mutate(DE = ifelse(logFC > 2 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <2 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("B cell zone vs. T cell zone in Lymph node (limma-voom)") +
  scale_color_manual(values = c("blue","gray","red")) +
  theme(text = element_text(size=15))
################################################################################
### VIII. Enrichment analysis          
################################################################################
#We will use MSigDB for the enrichment analysis. Since the data is quite large, a RData object is provided.
#You can skip the commented out line and directly load the RData object.

#Load the Msigdb genesets
#msigdb_hs <- getMsigdb(version = '7.2')
#msigdb_hs <- appendKEGG(msigdb_hs)
#sc <- listSubCollections(msigdb_hs)
#gsc <- c(subsetCollection(msigdb_hs, c('h')),
#         subsetCollection(msigdb_hs, 'c2', sc[grepl("^CP:",sc)]),
#         subsetCollection(msigdb_hs, 'c5', sc[grepl("^GO:",sc)])) %>%
#  GeneSetCollection()

#Filter the genesets that contain less than 5 genes.
#fry_indices <- ids2indices(lapply(gsc, geneIds), rownames(v), remove.empty = FALSE)
#names(fry_indices) <- sapply(gsc, setName)

#gsc_category <- sapply(gsc, function(x) bcCategory(collectionType(x)))
#gsc_category <- gsc_category[sapply(fry_indices, length) > 5]

#gsc_subcategory <- sapply(gsc, function(x) bcSubCategory(collectionType(x)))
#gsc_subcategory <- gsc_subcategory[sapply(fry_indices, length) > 5]
#fry_indices <- fry_indices[sapply(fry_indices, length) > 5]
#names(gsc_category) = names(gsc_subcategory) = names(fry_indices)

#Use fry to get the enriched pathways
#fry_indices_cat <- split(fry_indices, gsc_category[names(fry_indices)])
#fry_res_out <- lapply(fry_indices_cat, function (x) {
#  limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)
#})

#post_fry_format <- function(fry_output, gsc_category, gsc_subcategory){
#  names(fry_output) <- NULL
#  fry_output <- do.call(rbind, fry_output)
#  fry_output$GenesetName <- rownames(fry_output)
#  fry_output$GenesetCat <- gsc_category[rownames(fry_output)]
#  fry_output$GenesetSubCat <- gsc_subcategory[rownames(fry_output)]
#  return(fry_output)
#}

#fry_res_sig <- post_fry_format(fry_res_out, gsc_category, gsc_subcategory) %>%
#  as.data.frame() %>%
#  filter(FDR < 0.05) 

load("fry_res_sig.RData")

#Barplot of the pathways enriched in upregulated genes.
fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Up") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "red") +
  theme_bw() +
  coord_flip() +
  ggtitle("Up-regulated")

#Barplot of the pathways enriched in downregulated genes.
fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Down") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw() +
  coord_flip() +
  ggtitle("Down-regulated")
################################################################################

