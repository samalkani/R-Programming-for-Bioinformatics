# 1. Gene Expression Analysis - Microarrays

# Gene expression analysis is a cornerstone of bioinformatics and molecular biology,
# providing insights into cellular functions, regulations, and underlying mechanisms 
# of diseases. With the advancement of genomic technologies, microarrays have emerged
# as a powerful technique for measuring the expression levels of thousands of genes 
# simultaneously. This chapter aims to guide readers through the process of analyzing 
# microarray data using R, a widely-used programming language in bioinformatics, 
# emphasizing the tools and techniques that facilitate gene expression analysis.

# 1. A. Overview of Microarrays

# Microarrays are a high-throughput technology that allows researchers to analyze 
# gene expression profiles across various conditions or treatments. Each microarray 
# consists of thousands of DNA probes, each corresponding to a specific gene. 
# The basic principle involves hybridizing labelled RNA or cDNA samples to the 
# probes, where the fluorescent signal intensity at each spot indicates the level 
# of gene expression.

# 1. A. I. Types of Microarrays

# * __cDNA Microarrays__: Composed of known cDNAs and used for comparative analysis 
# of gene expression.

# * __Oligonucleotide  Microarrays__: Contain short, synthetic oligonucleotides 
# that provide high specificity and sensitivity for gene expression profiling.

# * __SNP Arrays__: Used to detect single nucleotide polymorphisms rather than 
# gene expression.

# 1. B. Data Acquisition and Pre-processing.

# Before diving into analysis, the first step is acquiring the raw microarray data.
# This may involve downloading data from repositories such as the __Gene Expression 
# Omnibus (GEO)__ or __The European Nucleotide Archive (ENA)__. Once data is 
# obtained, pre-processing is crucial.

# 1. B. I. Loading Data into R

# R provides several packages to facilitate the import and management of microarray
# data. Below is a code snippet demonstrating how to load data from GEO data set 
# using the "GEOquery" package.

# Install and load necessary packages
# install.packages("BiocManager")
# BiocManager::install("GEOquery")
library(GEOquery)
library(Biobase)

# Load GEO dataset
geo_dataset <- getGEO("GSE57820", GSEMatrix = TRUE)
ls()

# Re-assign dataset name to something shorter
geo = geo_dataset
geo

# Title of study
geo[["GSE57820_series_matrix.txt.gz"]]@experimentData@title

# Abstract
geo[["GSE57820_series_matrix.txt.gz"]]@experimentData@abstract

# Split geo_dataset into pdata, fdata and expression data
pdata_geo=geo[["GSE57820_series_matrix.txt.gz"]]@phenoData
head(pdata_geo)
pdata_geo <- pData(pdata_geo)
head(pdata_geo)
fdata_geo = geo[["GSE57820_series_matrix.txt.gz"]]@featureData
head(fdata_geo)
expression_data = geo[["GSE57820_series_matrix.txt.gz"]]@assayData[["exprs"]]
head(expression_data, n = 5)

# 1. B. II. Quality Control

# Quality control is essential to ensure that the data is reliable and suitable 
# for analysis. This step involves identifying and removing low-quality samples 
# or outliers. Common tools for visualizaing data quality include box plots, 
# density plots and principal component analysis (PCA).

# Boxplot
boxplot(expression_data, main = "Quality Control: Boxplot of Expression Data")

# PCA for outlier detection
pca_result <- prcomp(t(expression_data), scale. = TRUE)
plot(pca_result$x, main = "PCA of Expression Data")

# 1. B. III. Normalization

# Normalization adjusts the data to account for systematic biases and ensures 
# that differences in expression levels are biological rather than technical 
# artifacts. Several methods exist, with the "RMA" (Robust Multi-array Average) 
# algorithm being a commonly used approach for background correction, normalization, 
# and summarization.

# Normalizing data using RMA method
# BiocManager::install("affy")
library(affy)
# normalized_data <- rma(geo_dataset[[1]])

# Error: unable to find an inherited method for function ‘probeNames’ for signature 
# ‘object = "ExpressionSet"’

# The error you are seeing:
  
# Error: unable to find an inherited method for function ‘probeNames’ for signature
# ‘object = "ExpressionSet"’ indicates that the rma function is expecting an input 
# of a different class (usually AffyBatch), but you are providing an object of class 
# ExpressionSet. The rma function from the affy package does not directly operate 
# on an ExpressionSet; it expects raw CEL files read using ReadAffy() or similar
# functions.

# Here’s a structured way to solve this:

# Download Raw CEL Files (preferred for rma)
# Example GEO accession
# geo_id <- "GSE57820"

# Download the GEO dataset including raw CEL files
# geo_data <- getGEO(geo_id, GSEMatrix = FALSE)

# Check if raw CEL files exist:
# names(geo_data)   # Should include CEL files under supplementary files

# Read CEL Files Before RMA rma works on AffyBatch objects, not ExpressionSet:
# Suppose CEL files are downloaded to current directory
# cel_files <- list.celfiles(".", full.names = TRUE)
# affy_data <- ReadAffy(filenames = cel_files)

# Normalize using RMA
# normalized_data <- rma(affy_data)

# If You Only Have ExpressionSet

# If you only have an ExpressionSet already, RMA has likely already been applied to the data.
# Use exprs() to access normalized expression values:

# Access expression matrix
normalized_data <- exprs(geo_dataset[[1]])
head(normalized_data)

# Applying rma() again is not necessary and will error out because ExpressionSet is not compatible with rma.

# 1. C. Differential Expression Analysis

# Once the data is pre-processed and normalized, the focus shifts to identifying
# differentially expressed genes (DEGs) across experimental conditions. The "limma" 
# package is a popular choice for this purpose, utilizing linear modelling and 
# empirical Bayes methods

# 1. C. I. Design Matrix

# Before performing differential expression analysis, define the experimental design
# with a design matrix that captures the conditions being compared.

# Create a design matrix

# GSM1394594 GSM1394595 GSM1394596 GSM1394597 GSM1394598 GSM1394599 = treatment samples

# GSM1394600 GSM1394601 GSM1394602 GSM1394603 GSM1394604 GSM1394605 = negative control samples


# Sample data setup (replace with your actual data)
sample_names <- paste0("Sample", 1:12)
sample_names
sample_ids <- geo[["GSE57820_series_matrix.txt.gz"]]@phenoData@data[["geo_accession"]]
sample_ids
titles <- geo[["GSE57820_series_matrix.txt.gz"]]@phenoData@data[["title"]]
titles
groups <- factor(rep(c("Treatment", "Control"), each = 6))
groups
pheno_data <- data.frame(SampleName = sample_names, Sample_id = sample_ids, Title = titles, Group = groups)
pheno_data
rownames(pheno_data) <- pheno_data$SampleName # set sample names as row names
rownames(pheno_data)

# Using the formula ~Group Default reference is first level, "Control"
mod <- model.matrix(~ pheno_data$Group, data = pdata_geo)
mod

# 1. C. II. Fit Model and Conduct Analysis

# After setting up the design matrix, fit the linear model for each gene and 
# apply statistical tests to identify DEGs.

# Load the Libraries
library(limma)

# Fitting the linear model
fit <- lmFit(normalized_data, mod)
fit <- eBayes(fit)
head(fit)
summary(fit)

# Extracting DEGs
results <- topTable(fit, adjust.method = "BH", sort.by = "P", number = 100)
head(results)

# 1. C. III. Result Interpretation - Volcano Plot

# The output from the differential expression analysis includes log-fold changes, 
# p-values and adjusted p-values. Using a volcano plot, one can visually represent 
# the results, highlighting significantly altered genes.

# Volcano Plot
volcano_data <- data.frame(rownames(results), logFC = results$logFC, results$P.Value, -log10(results$P.Value))
head(volcano_data) 
plot(volcano_data$logFC, volcano_data$negLogPval, pch=20, main="Volcano Plot", 
     xlab = "Log Fold Change", ylab = "Statistical Significance (-log10)")
abline(h = -log10(0.05), col="red")

# 1. C. III. Result Interpretation - Enhanced Volcano Plot

# Load Library
library(EnhancedVolcano)

# Volcano Plot
EnhancedVolcano(volcano_data, 
                lab = volcano_data$rownames.results., 
                x = "logFC",
                y = "results.P.Value", 
                pCutoff = 0.05, 
                FCcutoff = 1.0,
                # Change axis limits
                xlim = c(-1.05, 1.5),     # X-axis: log2 fold change from -5 to +5
                ylim = c(0, 11.6),      # Y-axis: -log10(p-value) from 0 to 20
                title = "Volcano Plot",
                subtitle = 'Differential Expression Results')

# 1. D. Functional Enrichment Analysis

# Once DEGs are identified, functional enrichment analysis can elucidate biological 
# pathways and processes significantly associated with the list of altered genes. 
# Tools like "clusterProfiler" facilitate pathway analysis, such as Gene Ontology 
# (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG).

# Install "clusterProfiler" and "BSgenome.Hsapiens.UCSC.hg38" packages
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.egGO")
# BiocManager::install("topGO")
# BiocManager::install("Rgraphviz")

# Install gprofiler2 if not already installed
# if (!requireNamespace("gprofiler2", quietly = TRUE)) {
#   install.packages("gprofiler2")
# }


# Load libraries
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gprofiler2)
library(topGO)
library(Rgraphviz)

# Example: Illumina probe IDs (replace with your own list)
illumina_ids <- c(rownames(results))

# Convert Illumina IDs to Entrez Gene IDs for human (hsapiens)
# The 'target' parameter specifies the output namespace
# The 'organism' parameter specifies the species
conversion_result <- gconvert(
  query = illumina_ids,
  organism = "hsapiens",
  target = "ENTREZGENE_ACC",  # Entrez Gene ID
  filter_na = TRUE            # Remove rows with no match
)

# Display the conversion table
head(conversion_result, n = 2)

# Performing GO analysis
go_results <- enrichGO(gene = conversion_result$name, OrgDb = org.Hs.eg.db, 
                       keyType = "SYMBOL", pAdjustMethod = "BH")

plotGOgraph(go_results)

# 1. E. Visualization

# Effective visualization aids in understanding the results of gene expression
# analysis. In addition to the volcano plot, heatmaps, and clustering dendrograms
# can illustrate expression patterns across samples.

# Load library
library(pheatmap)

# Heatmap Generation
heatmapData <- normalized_data[rownames(results[results$adj.P.Val < 0.05,]),]

pheatmap(heatmapData, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")

# By integrating tools and techniques for data pre-processing, normalization, 
# differential expression analysis, functional enrichment and visualization, 
# researchers can derive meaningful biological insights from gene expression 
# profiles. As bioinformatics continues to evolve, leveraging R and its extensive 
# Bioconductor ecosystem will remain essential for unlocking the mysteries of gene
# regulation in health and disease.

# 2. Pre-processing Microarray Data with Affy and Oligo

# This chapter focuses on the pre-processing of microarray data using two popular 
# R packages: __affy__ and __oligo__. By the end of this chapter, you will have 
# a solid understanding of how to pre-process microarray data, enhancing the 
# reliability and interpretability of your analyses

# 2. A. Overview of Microarray Data

# Microarray experiments involve measuring the expression levels of thousands of 
# genes across different samples. The raw output from these experiments typically 
# consists of intensity values for each probe in the microarray. Different platforms 
# and technologies often yield raw data in formats such  as CEL files (for 
# Affymetrix arrays) or TXT files (for Illumina arrays). Standard pre-processing
# steps help transform these raw data into normalized gene expression values that
# are ready for downstream analysis.

# 2. B. Pre-processing Workflow

# A typical microarray data pre-processing workflow includes the following steps.

# * __Importing raw data__

# * __Quality assessment__

# * __Background correction__

# * __Normalization__

# * __Summarization__

# * __Filtering and transformation__

# This chapter will delve into each of these steps, detailing how to implement 
# them using __affy__ and __oligo__.

# 2. C. Importing Raw Data

# 2. C. I. Using Affy.

# The __Affy__ package is designed for processing Affymetrix micro-array data.
# The first step is to load the required libraries and import the CEL files:

# Install "oligo" if not already installed
BiocManager::install("oligo")
BiocManager::install("illuminaio")

# Load necessary libraries
library(affy)
library(oligo)

# Set working directory to where CEL files are located
setwd("D:/Ajay Files/R Programming for Bioinformatics/Chapter 7/Microarrays")

# Read in the CEL files
rawData_affy <- ReadAffy()
rawData_affy

# 2. C. II. Using oligo

# The __oligo__ package is versatile and can handle almost any microarray data, 
# including Illumina and Agilent platforms. For instance, to import Illumina data, 
# you can use:

# Load necessary libraries
library(oligo)
library(illuminaio)

# Read in the IDAT files for Illumina arrays
# rawData_Ill <- read.idat("D:/Ajay Files/R Programming for Bioinformatics/Chapter 7/IDAT example files")

# Make sure to replace the path with the actual location of your data files.

# 2. D. Quality Assessment

# Quality control is essential in any micro-array experiment. We can visualize 
# the data using boxplots which help identify any outliers:

# Create boxplots for quality assessment
boxplot(rawData_affy, main = "Boxplot of Raw Data")

# 2. D. I. Affy Quality Metrics

# For Affymetrix arrays, you can use the "affyqc()" function to generate quality
# metrics:

# Generate quality metrics 
# affyQCMetrics <- affyQC(rawData_affy) # affyQC() doesn't work!

# Alternative
# BiocManager::install("affyPLM")
library(affyPLM)

fit <- fitPLM(rawData_affy)
fit
NUSE(fit)
RLE(fit)

# 2. D. II. Oligo Quality Metrics

# For Illumina and other arrays, the __oligo__ package provides similar tools, 
# but you may need to evaluate metrics specific to your platform.

# 2. E. Background Correction

# Background correction is a technique to adjust for any inherent noise in the 
# micro-array experiment. The __affy__ package provides several algorithms for 
# background correction.

# 2. E. I. Perform background correction on Affymetrix micro-array

# Run the MAS 5.0 algorithm, which performs background correction, normalization, 
# and summarization
bgCorrectedData_affy <- mas5(rawData_affy)
bgCorrectedData_affy # Expression set created

# 2. E. II. Perform background correction on Illumina micro-array

# In __oligo__, you can use the quantile normalization method to perform 
# background correction:

# Perform background correction and normalization
# bgCorrectedData_Ill <- normalize(rawData_Ill)

# 2. F. Normalization

# Normalization is crucial for reducing systematic biases across arrays. The 
# __affy__ package offers various normalization methods like RMA (Robust Multi-
# array Average)

# Normalize the data using RMA 
# The rma function is available on oligo and affy packages so you have to specify 
# which package to get the rma function from. The rma function converts affy batch 
# data into an Expression set but you also have to specify you want to normalize 
# the data. You can't normalize an expression set (bgCorrectedData_affy) using rma,
# you have to use an affy batch data set(rawData_affy).
normalizedData <- affy::rma(rawData_affy, normalize = TRUE)
normalizedData

# With _oligo__, normalization is typically integrated into the data importing 
# process. Ensure to follow the specific normalization method recommended by 
# your array type.

# 2. G. Summarization

# Summarization transforms probe-level data to a gene-expression measure. In the 
# __affy__ package, we can achieve this through "rma()" or "mas5()" depending on 
# your analysis needs:

# Summarize the data to gene-level 
summarizedData <- exprs(normalizedData)
head(summarizedData, n = 2)

# In __oligo__, summarization is handled internally during normalization:

# Extract summarized data for further analysis
# geneExpressionData <- exprs(normalizedData)

# 2. H. Filtering and Transformation

# Once you have a summarized expression matrix, filtering lowly-expressed genes 
# can enhance the quality of your downstream analysis. It often involves setting 
# a threshold for expression levels:

# Dimensions of summarized data prior to filtering
dim(summarizedData)

# Filtering
threshold = 3
filteredData <- summarizedData[rowMeans(summarizedData) > threshold,]
head(filteredData, n = 2)

# Dimensions of summarized data after filtering
dim(filteredData)

# It's essential to choose a threshold based on biological relevance and 
# experimental design.

# We explored the systematic process of pre-processing microarray data using 
# __affy__ and __oligo__ in R. We walked through the steps of importing data,
# performing quality assessment, background correction, normalization, 
# summarization and filtering of gene expression data. Each of these steps is 
# essential to ensure that the data analyzed is of high quality and can yield 
# meaningful biological insights.

# 3. Differential Expression Analysis with limma

# The objective of this analysis is to identify genes whose expression levels 
# significantly differ between defined experimental groups, aiding in the 
# understanding of biological processes, disease mechanisms, and potential 
# therapeutic targets.

# In the realm of bioinformatics, one of the most widely used tools for differential
# expression analysis is the "limma" packages, which stands for "linear models 
# for microarray data". Originally developed for microarray data, limma has been
# successfully adapted for RNA-seq data, particularly through the use of voom 
# transformation. This chapter covers the basic concepts and step-by-step 
# procedures for conducting differential expression analysis using limma in R, 
# providing practical examples and insights for researchers in the field.

# 3. A. Overview of the limma Package

# The limma package is implemented in R, offering a suite of functions for the 
# analysis of gene expression data. It employs linear modelling to assess 
# differential expression and is particularly noted for its robust statistical 
# framework, allowing for the inclusion of complex designs and the incorporation
# of multiple testing corrections. The core functionalities of limma can be best
# summarized as follows:

# * __Modelling__: Fit linear models to expression data

# * __Hypothesis Testing__: Test for differential expression using empirical 
# Bayes methods.

# * __Multiple Testing Adjustment__: Control the false discovery rate (FDR) and 
# provide adjusted p-values.

# 3. B. Installation and Loading of limma

# To get started with limma, you'll first need to install the package from 
# Bioconductor, as it is not available on CRAN. The installation process can be 
# accomplished via the following commands in R:

# Install BiocManager if you haven't already.
# install.packages("BiocManager")

# Install limma
# BiocManager::install("limma")

# Once installed, load the package into your R session:
library(limma)

# 3. C. Data Preparation

# The initial step in any differential expression analysis is data preparation.
# For this example, we will consider a hypothetical RNA-seq dataset. The dataset
# might be structured with genes as rows and samples as columns, often stored in 
# a data frame or matrix format. Here's a basic structure of how your data might 
# look:

# Random seed for reproducibility
set.seed(123)

# Create 100 rows X 10 columns count matrix with negative binomial distribution
gene_counts <- matrix(rnbinom(n = 1000, mu = 10, size = 1),
                      nrow = 100, ncol = 10)
rownames(gene_counts) <- paste0("Gene_", 1:100)
colnames(gene_counts) <- paste0("Sample_", 1:10)

gene_counts

# Creating a data frame to hold sample information
sample_info <- data.frame(Group = rep(c("Control", "Treatment"), each = 5))
sample_info 

# 3. D. Normalization and Transformations

# Before fitting models, it's critical to normalize the count data to account for 
# differences in sequencing depth and composition. For RNA-seq data, the voom 
# function is used:

# Load edgeR library
library(edgeR) # This is often  required for voom

# Creating a DGEList object for normalization
dge <- DGEList(counts = gene_counts)
dge

# Applying voom transformation
v <- voom(dge, design = model.matrix(~ Group, data = sample_info))
v

# 3. E. Linear Model Fitting

# Once the data is normalized, we can fit a linear model to the transformed data.
# The design matrix describes the structure of the experimental groups, which 
# will be used in the modelling process.

# Fit the linear model
fit <- lmFit(v, design = model.matrix(~ Group, data = sample_info))
fit

# 3. F. Empirical Bayes Moderation

# To improve the estimation of the variance and its reliability, we apply 
# empirical Bayes moderation:

fit <- eBayes(fit)
fit

# 3. G. Identifying Differentially Expressed Genes

# After fitting the model and applying empirical Bayes methods, we can create a 
# table of differentially expressed genes (DEGs) using the "topTable" function. 
# This function allows for the extraction of genes based on adjusted p-values:

result <- topTable(fit, coef = 2, adjust.method = "BH", number = Inf)
head(result)

# In this case, "coef = 2" indicates that we want to compare the Treatment group
# against the control group. The results can be filtered based on significance
# levels, typically a threshold of FDR < 0.05 is used to declare differential 
# expression

# 3. H. Visualization

# Visualization is crucial in understanding the results of differential expression 
# analyses. Common techniques include volcano plots and heatmaps.

# Install package
# BiocManager::install("EnhancedVolcano")

# Load Library
library(EnhancedVolcano)

# Volcano Plot
EnhancedVolcano(result, 
                lab = rownames(result), x = "logFC",
                y = "P.Value", pCutoff = 0.05,
                title = "Volcano Plot of Differential Expression")

# Heatmaps can be created to visualize the expression patterns of DEGs across 
# samples:

# load Library
library(pheatmap)

# Select top DEGs for heatmap
top_genes <- head(row.names(result[result$P.Value < 0.05,]), 20)
top_genes

# Heatmap
pheatmap(v$E[top_genes,], clustering_distance_rows = "correlation", 
         clustering_distance_cols = "euclidean")

# This chapter introduced the fundamental concepts and procedures for conducting 
# differential expression analysis using the limma package in R. Through a 
# systematic approach encompassing data preparation, normalization, model fitting, 
# hypothesis testing and visulization, researchers can leverage limma to glean 
# biological insights from gene expression data.


