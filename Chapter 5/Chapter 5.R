# Chapter 5: Data Visualization for Bioinformatics

# Bioinformatics involves the analysis of vast amounts of diverse data, including
# genomic sequences, protein structures, metabolic pathways and gene expression
# profiles. Effective visualization techniques can simplify the interpretation of 
# these data sets, revealing patterns, trends and relationships that are critical 
# for biological discoveries and innovations. In this chapter, we will explore 
# various data visualization techniques relevant to bioinformatics, the tools and 
# libraries available for creating visualizations, and case studies demonstrating 
# the application of these techniques in real-world scenarios.

# 1. Importance of Data Visualization in Bioinformatics

# 1. A. Understanding Multidimensional Data.

# Biological data is often high-dimensional and complex. Genomic data, for instance,
# many involve thousands of genes with various expression levels across multiple 
# conditions or time points. Visualization aids in reducing this complexity by 
# providing a graphical representation that allows researchers to explore and 
# identify significant trends and anomalies. Tools like Principal Components 
# Analysis (PCA) plots and t-SNE (t-distributed Stochastic Neighbor Embedding) 
# offer intuitive representations that help on dimensionality reductions, making
# it easier to visualize relationships among different data sets.

# 1. B. Communicating Results

# The Communicating Results

# The communication of results to stakeholders, including researchers, clinicians, 
# and the public, is vital for advancing scientific knowledge and translational 
# medicine. Well-designed visualizations can convey complex bioinformatics findings 
# effectively, allowing for clearer interpretation and engagement.

# Visualization can transform statistical outputs into more accessible formats, 
# making it easier for non-experts to grasp essential discoveries.

# 1. C. Hypothesis Generation

# Visual exploration of data can lead to new hypotheses and questions that might 
# not arise through traditional analytical methods alone. For instance, cluster 
# heat maps can visually demonstrate co-expression patterns of genes in different 
# samples, suggesting potential regulatory relationships that merit further 
# investigation.

# 2. Common Data Visualization Techniques in Bioinformatics

# 2. A. Heat maps

# Heat maps are widely used in bioinformatics to represent data matrices, such 
# as gene expression levels across different samples. They allow for the visualization
# of hierarchies and patterns through color gradients, enabling quick identification 
# of clusters or outliers. Common libraries like Seaborn in Python or ggplot2 in 
# R facilitate the creation of informative heat maps with ease.

# 2. B. Scatter Plots

# Scatter plots are effective for visualizing correlations between two continuous 
# variables. In genomics, scatter plots can illustrate the relationship between
# gene expression levels in different conditions, aiding in the identification
# of potential biomarkers or therapeutic targets. Additionally, bubble plots can 
# enhance scatter plots by adding a third variable through the size and color of 
# the points.

# 2. C. Box Plots

# Box plots are useful for comparing distributions of gene expression across 
# different categories, such as treatment groups or phenotypes. They succinctly
# represents measures of central tendency and variability, revealing insights 
# into data distribution and identifying potential outliers.

# 2. D. Network Graphs

# Biological systems are inherently interconnected, and visualizing these 
# relationships is crucial for understanding complex interactions. Network graphs 
# can depict protein-protein interactions, gene regulatory networks, and metabolic 
# pathways, helping researchers visualize and explore the interconnectedness of 
# biological entities.

# 2. E. Genomic Visualizations

# Genomic data can be visualized through various methods, including genome browsers,
# which allow for the exploration of genomes at different scales. Tools such as IGV
# (Integrative Genomics Viewer) provide interactive interfaces for examining 
# sequence alignments, variant calls and annotations in a user-friendly format.

# 3. tools and libraries for Data Visualization.

# A wide array of tools and libraries are available for generating visualizations 
# in bioinformatics:

# __Python Libraries__:

# * __Matplotlib__: A foundational library for creating static, animated and 
# interactive visualizations in Python.

# * __Seaborn__: Built on Matplotlib, it provides a high-level interface for 
# drawing attractive and informative statistical graphics.

# * __Ploty__: Facilitates the creation of interactive plots and dashboards, 
# making data explorations more engaging.

# __R Packages__:

# * __ggplot2__: A powerful and flexible visualization package based on the 
# Grammar of Graphics, widely used in the bioinformatics community.

# * __ComplexHeatmap__: Specifially designed for creating complex heatmaps and 
# annotations.

# * __Visualizations Platforms__:

# * __Tableau__: A business intelligence tool that allows users to create 
# interactive, shareable dashboards. Although not bioinformatics-specific, it can
# be leveraged for biological data visualization.

# * __Cytoscape__: A platform for visualizing complex networks and integrating
# these with any type of attribute data, commonly used in systems biology.

# 4. Case Studies

# 4. A. Case Study 1: The Genomic Expression of Cancer

# In a study analyzing the expression levels of thousands of genes in tumor versus
# normal tissues, researchers used heatmaps to display differential expression patterns.
# Clustering techniques revealed distinct groups of genes associated with specific
# cancer types, which were further investigated for potential diagnostic and 
# therapeutic applications.

# 4. B. Case Study 2: Protein Interaction Networks

# Researchers investigating the interactions of proteins related to Alzheimer's 
# disease used network graphs to visualize the complex interplay between various 
# proteins. This visualization facilitated the identification of key regulatory 
# proteins that may serve as targets for new therapeutic strategies.

# 5. Best practices for Effective Visualization

# * __Clarity__: Ensure that visualizations clearly convey the intended message.
# Use appropriate scales, labels and legends to avoid misinterpretation.

# * __Simplicity__: Avoid over-complicating visualizations. Focus on conveying 
# the primary finding without unnecessary embellishments.

# * __Interactivity__: Where possible, provide interactive visualizations that 
# allow users to delve deeper into the data and explore different aspects.

# * __Consistency__: Use consistent color schemes and styles across related 
# visualizations that allow users to delve deeper into the data and explore 
# different aspects.

# By leveraging various visualization techniques and tools, researchers can 
# uncover biological insights, communicate their findings effectively, and engage
# in hypothesis generation. As the field of bioinformatics continues to evolve, 
# so too will the methods and technologies for data visualization, enabling deeper
# understanding of biology at every level.

# 6. Creating Effective Plots for Biological Data with ggplot2.

# This chapter will guide you through the fundamentals of using ggplot2 for 
# biological data, including how to prepare your data, create basic plots, 
# customize visualizations, and incorporate best practices for data presentation.

# 6. A. Understanding ggplot2 Basics

# Before we dive into creating plots, it's crucial to understand how ggplot2 
# organizes and structures plotting. At its core, ggplot2 builds plots layer by 
# layer, allowing for immense flexibility. The basic structure of a ggplot 
# command follows this syntax:

# ggplot(data, aes(x = ..., y = ...,)) + geom_(type)(parameters)

# __data__: The data frame containing your biological data

# __aes__: Aesthetic mappings that specify how data should be mapped to visual 
# properties (like x and y coordinates, colors, sizes).

# __geom__: The geometric object that represents the data (e.g., points, lines, 
# bars). To start using ggplot2, you must first install and load the package:

# 6. A. I. Example

# Load Libraries
library(devtools)
library(Biobase)
library(tidyverse)
library(ggplot2)

# 6. B. Preparing Your Biological Data

# Before Visualizing your data, preparation is key. Often, biological data sets 
# can be messy or not structured in a way that's immediately useful for plotting. 
# The folowing steps are essential

# 6. B. I. Data Cleaning

# Ensure your data set does not contain NA values or erroneous entries. Use 
# functions such as 'na.omit()' or 'tidyr::drop_na()'. Check your data types and 
# transform them as needed, using 'dplyr::mutate()' for conversions.

# Load in data from a connection
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
ls()
bot = bottomly.eset
bot
pdata_bot=pData(bot)
head(pdata_bot)
fdata_bot = featureData(bot)
head(fdata_bot)
edata = exprs(bot)
head(edata, n=1)

# Check genomic data for NAs
is.na(edata)[1,]

# Change expression set (matrix) to a data frame
edata1 = as.data.frame(edata)
edata2 = as_tibble(edata)
edata3 = data.frame(edata)

# Log transform edata
log_edata = log2(edata1 +1)
head(log_edata)

# 6. B. II. Data transformation

# Sometimes, you may need to transform your data for analysis. This might involve
# calculating means, variances or creating new variables. The 'dplyr' package is 
# excellent for data manipulation:

# Load "dplyr" library
library(dplyr)

# Data Summarization

# Group the 21 samples into the two different strains
c57BL_6J = group_by(log_edata[,1:10])
DBA_2J = group_by(log_edata[,11:21])
genes = row.names(log_edata)

# Split each of the samples within the two strains into separate lists and concatenate them together
sample_1 = c(c57BL_6J[,1], c57BL_6J[,2], c57BL_6J[,3], c57BL_6J[,4], c57BL_6J[,5], 
             c57BL_6J[,6], c57BL_6J[,7], c57BL_6J[,8], c57BL_6J[,9], c57BL_6J[,10],
             DBA_2J[,1], DBA_2J[,2], DBA_2J[,3], DBA_2J[,4], DBA_2J[,5], DBA_2J[,6],
             DBA_2J[,7], DBA_2J[,8], DBA_2J[,9], DBA_2J[,10], DBA_2J[,11])

# Un-list the large list with 21 elements
sample_1 = unlist(sample_1, use.names=FALSE)

# Create an atomic vector
sample_1 = as.vector(sample_1)

# Create a Sample Data Frame with the two strains as separate categories within the Group variable
gene_data <- data.frame(
  Group = c(rep("C57BL_6J", each = 365360), rep("DBA_2J", each = 401896)),
  genes = rep(genes, times = 21),
  Expression = sample_1)

head(gene_data)
tail(gene_data)

# Calculate the mean gene expression by the two different strains
mean_expression_by_group <- gene_data %>%
  group_by(Group) %>%
  summarize(
    Mean_Expression = mean(Expression, na.rm = TRUE), # Calculate mean, ignoring NAs
    Num_Samples = n() # Optional: count samples in each group
  )

print(mean_expression_by_group)

# 6. C. Basic Plots with ggplot2

# With your data prepared, you can start creating plots. This section covers the 
# creation of some common plot types suitable for biological data.

# 6. C. I. Scatter Plots

# Scatter plots are used to examine relationships between two continuous variables
# For example, you could plot gene expression between strain c57BL_6J and DBA_2J

# Data Visualization
library(ggplot2)

# NB - had to truncate the gene expression levels for the DBA_2J strain because it had one additional sample vs c57BL_6J
ggplot() + geom_point(aes(x = gene_data$Expression[1:365360], y = gene_data$Expression[365361:730720])) + 
  labs(title = "Gene Expression of  strain c57BL_6J vs strain DBA_2J", 
       x = "Gene Expression of strain c57BL_6J", y = "Gene Expression of strain DBA_2J")

# 6. C. II. Bar plots

# Bar plots are ideal for comparing categorical data. For instance, comparing 
# average weights of different species:

species_wt <- read.csv("Species_Weight.csv", header = TRUE, stringsAsFactors = FALSE)
species_wt

ggplot(species_wt, aes(x = Species, y = Weight.kg., fill = Species)) +
  geom_bar(stat = "summary", fun = "mean") +
  labs(title = "Average Weight by Species", x = "species",
       y = "Average Weight(kg)") + theme_classic()

# 6. C. III. Box plots

# Boxplots are useful to visualize the distribution of a continuous variable 
# across different categories. They can effectively display variance and outliers:

ggplot(gene_data, aes(x = Group, y = Expression)) + geom_boxplot(outlier.colour = "red") + 
  labs(title = "Distrubution of gene expression levels by strain", x = "Strain", 
       y = "Gene Expression") + theme_light()

# 6. D. Customizing Your Plots

# Customizing your plots can enhance their clarity and aesthetic appeal. Here 
# are some common customizations:

# 6. D. I. Themes

# ggplots2 comes with several built-in themes that can modify the overall 
# appearance of your plots. Use "theme_*()" functions to easily adjust your plot's 
# aesthetics:

ggplot(gene_data, aes(x = Group, y = Expression)) + geom_boxplot(outlier.colour = "red") + 
  labs(title = "Distrubution of gene expression levels by strain", x = "Strain", 
       y = "Gene Expression") + theme_linedraw()

# 6. D. II. Color and Fill

# Using distinct colors helps differentiate categories. You can specific colors 
# directly in the "aes()" function:

ggplot(species_wt, aes(x = Species, y = Weight.kg., fill = Species, colour = Species)) +
  geom_bar(stat = "summary", fun = "mean") +
  labs(title = "Average Weight by Species", x = "species",
       y = "Average Weight(kg)") + theme_classic()

# 6. D. III. Labels and titles

# Always include informative titles and axis labels to guide the audience. Use the 
# "labs()" function or "ggtitle()" for titles.

ggplot(gene_data, aes(x = Group, y = Expression, fill = Group, color = Group)) + geom_boxplot(outlier.colour = "red") + 
  labs(title = "Distrubution of gene expression levels by strain", x = "Strain", 
       y = "Gene Expression") + theme_linedraw()

# 7. Best Practices in Scientific Visualization

# When visualizating biological data, adherence to best practices enhances data 
# interpretation:

# * __Keep it Simple__: Avoid clutter and focus on the data.

# * __Use Appropriate Scales__: Use logarithmic scales if data spans multiple 
# orders of magnitude.

# * __Incorporate Legends Wisely__: Ensure the legend is informative and doesn't
# obscure data.

# * __Test Different Presentations__: Sometimes, a different plot type might 
# better convey your message.

# By preparing data appropriately, mastering basic plots, customizing aesthetics, 
# and following best practices, you can create informative and visually appealing 
# graphics that enhance your research storytelling. ggplot2 is a versatile tool 
# that, once mastered, can significantly improve the presentation of your biological
# data analyses. Happy plotting!

# 8. Visualizing Genomic Data with Gviz and Other Tools

# One of the key challenges in genomics is the effective visualization of complex
# genomic data, which is essential for interpretation and communication of research
# findings. In this chapter, we will explore the powerful tools available in R, 
# focusing on the Gviz package, as well as other visualization techniques that 
# complement Gviz

# 8. A. Understanding Genomic Data.

# Genomic data can comprise a wide variety of information, including gene 
# annotations, expression levels, variant data can help researchers identify 
# patterns, test hypotheses, and communicate their findings more clearly.
# For instance, visual representations of genomic data can highlight gene 
# structures, expressional differences between conditions, and the localization 
# of genomic features.

# 8. B. Getting Started with Gviz

# 8. B. I. Installation and Setup

# To begin using Gviz, you need to have R and Bioconductor installed. If you 
# haven't already set it up, you can find installation instructions on the 
# Bioconductor website. After setting up Bioconductor, Gviz can be installed with
# the following commands:

# if(!requireNamespace("BiocManager", quietly = TRUE))
#    installed.packages("BiocManager")

# BiocManager::install()

# BiocManager::install(c("SparseArray", 'GenomicRanges', 'S4Arrays'), 
# type = "source", force = TRUE)

# BiocManager::install("Gviz")

# 8. B. II. Basic Concepts of Gviz

# Gviz is a flexible R package designed specifically for visualizing genomic data.
# It allows users to create a variety of plots, including gene models, ideograms, 
# and coverage plots, while integrating different layers of genomic information.

# 8. B. III. Creating Your First Plot

# Let's walk through the steps to create a simple genomic plot using Gviz. 
# The following example shows how to visualize gene annotation data.

# Load Library
library(Gviz)

# Create a GeneRegionTrack
geneTrack <- GeneRegionTrack("D:/Ajay Files/R Programming for Bioinformatics/Chapter 5/athal_genes.gtf")

# Create a Genome Axis Track
genomeTrack <- GenomeAxisTrack()

# Create the plot
plotTracks(list(genomeTrack, geneTrack), from = 1, to = 100000)

# 8. B. IV. Customizing Your Plots

# Gviz offers extensive customization options for enhancing the aesthetics of 
# your plots. This includes changing colors, labels, and track heights. For example, 
# the following code snippet shows how to customize the colors and add a title.

plotTracks(list(genomeTrack, geneTrack), from = 1, to = 100000,
           main = "Genomic Visualization", background.fill = "lightgrey",
           col = "steelblue")

# 9. Integrating Other Visualization Tools

# While Gviz provides a comphrehensive approach to genomic data visualization, 
# other R packages also complement Gviz well. Here are a few notable examples:

# 9. A. ggplot

# "ggplot" is a versatile plotting system that can be used for a variety of data 
# visualization tasks, including genomic data. Combining Gviz with ggplot2 can 
# enhance the interactivity and customize plots beyond the capabilities of Gviz 
# alone. For example, you can extract data from Gviz, manipulate it using "dplyr",
# and visualize it using "ggplot2"

# 9. B. Complex Heatmap

# "ComplexHeatmap" is another powerful R package for visualizing complex data 
# matrices, specifically tailored for high-dimensional genomic data. It allows 
# users to create annotated heatmaps that can include various track annotations, 
# providing a great way to visualize gene expression data alongside genomic 
# features.

# 9. C. Plotly

# "plotly" is focused on creating interactive plots and can be effectively used 
# with Gviz outputs. By converting static Gviz plots into interactive web visuals, 
# researchers can explore their data in more detail.

# 10. Case Studies

# 10. A. Case Study 1: Differential Expression Analysis

# In a typical analysis pipeline, one might start with statistical tests to 
# identify differentially expressed genes (DEG's). After obtaining a set of DEG's,
# Gviz can be used to visualize expression levels across different conditions. 
# The combination of Gviz with heatmaps and volcano plots provides a comprehensive 
# understanding of the data.

# 10. B. Case Study 2: Variant Visualization

# When working with genomic variants, Gviz can visualize SNPs or structural 
# variants alongside gene annotations and other genomic features. By overlaying 
# variant data on gene tracks, researchers can identify potential functional 
# impacts of variants in research questions related to disease.

# Visualizing genomic data is a crucial component of genomics research. With the 
# Gviz package and other visualization tools in R, researchers can create informative 
# and compelling visual representations of complex genomic information. By integrating 
# multiple visualization packages, you can take your analyses to the next level, 
# leveraging the strengths of each library to produce comprehensive, interactive 
# and insightful genomic visualizations







