# R Programming for Bioinformatics

# Chapter 1 Introduction to R for Bioinformatics

# 1. Installing Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE))
   installed.packages("BiocManager")BiocManager::install()

# 2. Install specific bioinformatics Packages
BiocManager::install("GenomicRanges")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")

# 3. Installing Additional Useful Packages

# 3. A. Tidyverse - A collection of R packages for Data Science
install.packages("tidyverse")

# 3. B. Data Visualization
install.packages("ggplot2")

# 3. C. Data Manipulation and Transformation
install.packages("dpyr")

# 3. D. Interactive Visualizations
install.packages("plotly")

# 4. Integrate version control with Git into your Bioinformatics workflow
# [https://git-scm.com/](https://git-scm.com/)

# 5. Testing Your Environment

# 5. A. Loading Packages
library(GenomicRanges)
library(DESeq2)

# 5. B. Check Versions
packageVersion("GenomicRanges")
packageVersion("DESeq2")
