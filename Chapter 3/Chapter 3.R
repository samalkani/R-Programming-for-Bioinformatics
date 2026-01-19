# R Programming for Bioinformatics

# Chapter 3 - Working with Biological Data Structures in R

# 1. Overview

# A comprehensive overview of key biological data structures used in R will be 
# covered here. Including the following:

# * Sequences
# * Phylogenetic Trees
# * Genomic data
# * Practical examples and Coding Techniques

# 1. A. Biological Data Types

# * __Sequences__ - DNA, RNA and protein sequences are fundamental data types in 
# molecular biology. They are typically represented as strings of nucleotides or
# amino acids

# * __Phylogenetic Trees__ - These trees depict evolutionary relationships among 
# species or genes. They are structured into branching diagrams and can be 
# analyzed to study evolution and speciation.

# * __Genomic Data__ - This includes data derived from whole genomes, such as 
# variants identified through sequencing, gene expression levels, and annotations

# * __Transcriptomic Data__ - Involving RNA sequencing data, it often requires 
# processing and normalization before analysis.

# * __Proteomic Data__ - Represents protein sequences, structures, or expression 
# data, useful for understanding biological functions and pathways.

# 1. B. Basic Data Structures in R

# The most prominent data structures provided by R include;

# * __Vectors__ - One-dimensional arrays that hold data of the same type. For 
# example, a vector to hold DNA nucleotide sequences.

# * __Lists__ - A collection of components that can contain different types of 
# data (including other lists), making them suitable for complex biological data.

# * __Data Frames__ - Two-dimensional, tabular structures where each column can 
# hold different types of data. Data frames are particularly useful for storing 
# tabular genomic data, such as gene expression values across different samples.

# * __Matrices__ - Two-dimensional arrays that can contain only one data type. 
# They are often used for numerical data, including statistical analysis of gene 
# expression.

# * __S3 and S4 Classes__ - The object-oriented programming systems allow for 
# encapsulation and management of more complex biological data types (like 
# phylogenetic trees)

# 1. C. Working with Biological Sequences

# The "Biostrings" package from the Bioconductor project used in R is useful for 
# managing nucleotide and protein sequences.

# 1. C. I. Example - Analyzing DNA Sequences.

# Installing Bioconductor
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   installed.packages("BiocManager")BiocManager::install()

# Install specific bioinformatics Packages
# BiocManager::install("GenomicRanges")

# Loading Library
library(Biostrings)

# Create a DNA sequence
dna_seq <- DNAString("AGCTTAGCAGTA")
print(dna_seq)

# Reverse Complement
rev_comp <- reverseComplement(dna_seq)
print(rev_comp)

# Find a specific motif
motif <- "AGC"
match <- matchPattern(motif, dna_seq)
print(match)

# Tasks central to many biological analyses performed here include;

# * Creating a DNA sequence
# * Computing its reverse complement
# * Searching for specific motif

# 1. D. Working with Phylogenetic Trees

# Phylogenetic trees can be visualized and manipulated using the "ape" package 
# in R. This package offers tools for tree interpretation and modification.

# 1. D. I. Example - Creating and Plotting a Phylogenetic Tree.

# Install and load the "ape" package
install.packages("ape")

# Load Library
library(ape)

# Create a simple phylogenetic tree
tree <- read.tree(text = "((A,B),(C,D));")
plot(tree)

# Add branch labels
nodelabels()

# The above code performs the following;

# * A basic phylogenetic tree
# * Phylogenetic tree visualization
# * Addition of branch labels

# 1. E. Handling Genomic Data with "Genomic Ranges"

# Genomic data analysis often requires the organization and manipulation of vast 
# ranges of genomic coordinates. The "Genomic Ranges" package helps manage such 
# data efficiently, allowing for the representation of ranges on the genome.

# 1. E. I. Example - Working with "Genomic Ranges".

# Install and load "Genomic Ranges" package
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   installed.packages("BiocManager")BiocManager::install()

# BiocManager::install("GenomicRanges")

# Load Library
library(GenomicRanges)

# Create Genomic Ranges
gr <- GRanges(seqnames = Rle(c("chr1", "chr1", "chr2")), 
              ranges = IRanges(start = c(1, 100, 150), end = c(50, 200, 300)),
              strand = Rle(c("+", "-", "+")))
# Show the genomic ranges
print(gr)

# The example above performs the following;

# * Definition of genomic ranges
# * Essential for overlaying different types of genomic data
# * e.g. mutations, gene annotations and expression profiles

# Ability to handle sequences, trees, and genomic ranges effectively is 
# foundational for research across various fields, from genetics to ecology.

# 2. Understanding Genomic Ranges, Summarized Experiment and Related Classes.

# Among these, the __Genomic Ranges__ and __Summarized Experiment__ classes play
# a pivotal role in supporting the storage, retrieval, and analysis of genomic 
# features and experimental results. This chapter delves into the structure and 
# functions of these classes, providing a foundation for understanding their 
# applications in genomic data analysis.

# 2. A. The Genomic Ranges Package

# 2. A. I. Overview

# The __Genomic Ranges__ package in Bioconductor offers data structures to 
# represent genomic intervals (i.e., regions of DNA, RNA, or other genomic features) 
# and facilitates their manipulation. It is designed for efficient handling and 
# analysis of genomic data, leveraging R's capabilities to operate on large data sets.

# 2. A. II. Key Classes

# The central class in the __Genomic Ranges__ package is "GRanges", which 
# represents a set of genomic ranges. Let's explore its components:

# __Seqnames__: The chromosome or sequence name on which the range is located 
# (e.g. "chr1", "chr2").

# __Ranges__: The actual start and end positions of the genomic features.

# __Strand__: The DNA strand ("+", "-", or "*") to indicate which side of the 
# double helix the feature lies.

# __Metadata__: Additional annotations or information related to the genomic 
# features (e.g., gene names, descriptions).

# 2. A. III. Creating GRanges Objects

# Creating "GRanges" Object is straightforward. Here's an example to illustrate:

# Define the genomic ranges
seqnames <- c("chr1", "chr1", "chr2") 
start <- c(100, 500, 200)
end <- c(200, 600, 300)
strand <- c("+", "-", "*")

# Create the GRanges object
gr <- GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end),
              strand = strand)

# Inspect the GRanges object
print(gr)

# This code creates;

# * "GRanges" object
# * Containing 3 genomic intervals
# * Each associated with specific chromosomes and strands

# 2. A. IV. Common Operations

# The __Genomic Ranges__ package provides numerous functions for manipulating 
# genomic ranges, such as;

# * __Find Overlaps__: Identify ranges that overlap with others using the 
# "find Overlaps" function.

# * __Merge__: Combine overlapping genomic ranges with the "reduce" function.

# * __Counts__: Count how many ranges fall within specified genomic intervals.

# These operations can be conducted efficiently, making it feasible to analyze large data sets.

# 2. B. Summarized Experiment Class

# 2. B. I. Overview.

# While "GRanges" is focused on genomic features, the __Summarized Experiment__
# class serves as a data container for storing high-dimensional genomic data, 
# such as count matrices derived from RNA-seq experiments. This class allows 
# researchers to manage both the experiment-level metadata and the genomic features
# associated with the data.

# 2. B. II. Structure

# A __Summarized Experiment__ object contains:

# * __Assays__: A list of matrices storing genomic data, typically with rows 
# representing genomic features and columns representing samples.

# * __Row Data__: Metadata associated with the rows (features), often in the form 
# of a "Data Frame".

# * __Column Data__: Metadata associated with the columns (samples), also in a 
# "Data Frame".

# 2. B. III. Creating Summarized Experiment Objects

# Here's how you can create a __Summarized Experiment__ object.

# Load Library
library(SummarizedExperiment)

# Sample Data (rpois = 100 data points - Poisson Dist.(Pd), Lambda = mean of Pd, nrow = 10 rows)
counts <- matrix(rpois(100, lambda = 10), nrow = 10)

# 10 features, 10 samples (gene & sample names)
row_data <- DataFrame(gene_id = paste0("gene", 1:10))

col_data <- DataFrame(sample_id = paste0("sample", 1:10))

# Create the Summarized Experiment object
se <- SummarizedExperiment(assays = list(counts = counts), rowData = row_data, 
                           colData = col_data)

# Inspect the object
print(se)

# This creates;

# * Summarized Experiment object
# * Includes matrix of counts
# * Associated metadata

# 2. B. IV Accessing Data

# The __Summarized Experiment__ class implements several methods for accessing 
# its components:

# * __"assay()"__: Extracts the assay data.
# * __"rowData() and colData()"__: Retrieve feature and sample metadata, 
# respectively

# Using these methods, one can easily handle complex experimental data, providing 
# a structured interface for analysis.

# 3. Combining GRanges and Summarized Experiment

# One of the powerful features of R is the ability to seamlessly integrate different
# data structures. In genomic studies, you often encounter situations where you need
# to analyze genomic feature alongside experimental data. This is where combining 
# "GRanges" with "Summarized Experiment" becomes invaluable.

# 3. A. Adding Genomic Ranges to Summarized Experiment

# To associate genomic ranges with experimental data, you can use "rowRanges()" 
# to assign "GRanges" object to your __Summarized Experiment__:

# Define the genomic ranges
seqnames <- c("chr1", "chr1", "chr2", "chr1", "chr1", "chr2", "chr1", "chr1", "chr2", "chr2") 
start <- c(100, 500, 200, 300, 400, 600, 700, 800, 900, 1000)
end <- c(200, 600, 300, 400, 500, 700, 800, 900, 1000, 1100)
strand <- c("+", "-", "*", "+", "-", "*", "+", "-", "*", "-")

# Create the GRanges object
gr <- GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end),
              strand = strand)

# Inspect the GRanges object
print(gr)

# Assign GRanges to Summarized Experiment object
rowRanges(se) <- gr
print(rowRanges(se))

# This way, the genomic context of the features is embedded directly within the 
# __Summarized Experiment__ object, allowing for integrated analyses.

# The __Genomic Ranges__ and __Summarized Experiment__ classes in R provide 
# researchers with powerful tools for managing and analyzing genomic data. The 
# ability yo represent genomic features through "GRanges" and store high-dimensional
# data in "Summarized Experiment" facilitates complex analyses and aids in deriving
# meaningful biological insights. Understanding and employing these classes can 
# significantly enhance the conduct of genomic studies, making R a versatile 
# choice for bioinformatics tasks

# 4. Importing and Managing Biological Data Sets in R

# In contemporary biological research, the proficiency to handle, analyze, and 
# visualize large data sets is indispensible. R, a powerful programming language
# and environment, provides many tools for importing and managing biological data 
# sets. This chapter focuses on how to effectively work with different forms of 
# biological data in R, emphasizing best practices for importing, cleaning, and
# maintaining data sets for subsequent analysis.

# 4. A. Understanding Biological Data Types

# Before diving into the mechanisms of importing biological data sets, it's crucial
# to understand the types of data commonly encountered in biological research. 
# Typically, biological data sets can be categorized into the following:

# * __Genomic Data__: Sequencing data, SNP data, or gene expression data (often 
# in formats like FASTA, VCF, or TXT).

# * __Proteomic Data__: Data pertaining to protein expressions, often stored in 
# formats like CSV or Excel spreadsheets.

# * __Metabolomic Data__: Data capturing metabolic profiles, usually structured 
# in tables.

# * __Ecological Data__: Information concerning species distributions, environmental
# variables, etc.

# Each data type may come in a distinct format, and R provides specific packages
# tailored to each of these formats

# 4. B. Importing Data into R

# 4. B. I. Using Base R Functions

# R offers a variety of base functions for importing data. Here are some common 
# approaches:

# 4. B. I. Using Base R Functions - CSV

# * __CSV Files__: The "read.csv()" function is used to load CSV files. The "header" 
# argument specifies whether the first row contains variable names.

data1 <- read.csv("Genomic_Ranges_CSV.csv", header = TRUE, stringsAsFactors = FALSE)
data1

# 4. B. I. Using Base R Functions - XLSX

# * __Excel Files__: Though R does not have built-in functionalities for Excel, 
# additional packages like "readxl" can be utilized.

library(readxl)
data2 <- read_excel("Genomic_Ranges_XLSX.xlsx")
data2

# 4. B. I. Using Base R Functions - TEXT

# * __Text Files__: The "read.table()" function can be employed for importing 
# tab-delimited text files

data3 <- read.table("Genomic_Ranges_TSV.txt", header = TRUE, sep = "\t", 
                   stringsAsFactors = FALSE )
data3

# 4. B. II. Utilizing Specialized Packages

# For specific data types, specialized packages are often more efficient:

# * __Bioconductor__: This is a repository of tools for the analysis of genomic 
# data and can handle formats like BAM, VCF and more. The "Biostrings" package, 
# for example, facilitates the import and manipulation of DNA sequences.

# 4. B. II. Utilizing Specialized Packages - "Biostrings" Package

# Load Library
library(Biostrings)

# Read DNA sequence data - FASTA file
dna_sequences <-readDNAStringSet("dna.example.fasta")

# Elements of S4 object of class "DNAStringSet"
dna_sequences@pool
dna_sequences@ranges
dna_sequences@elementType
dna_sequences@elementMetadata
dna_sequences@metadata

# * __tidyverse__: The "readr" package from the tidyverse collection provides 
# user-friendly functions such as;

# "read_csv()" and "read_tsv()" which are optimized for speed and data integrity.

# 4. B. II. Utilizing Specialized Packages - "readr" package

# Load library
library(readr)

# Read csv using a function from the "readr" package
data4 <- read_csv("Genomic_Ranges_CSV.csv")
data4

# 5. Data Cleaning and Management

# Once data is imported, cleaning and manageing the data set is crucial to ensure 
# accuracy in analysis. The following steps are foundational:

# 5. A. Exploring the Data

# Utilize functions like "head()", "str" and "summary()" to get an overview of the imported data
head(data4)
str(data4)
summary(data4)

# 5. B. Handling Missing Values

# Missing data is a common issue in biological data sets. Use functions like 
# "na.omit()" or "is.na()" to identify and handle missing values.

# Read in data set with missing values
data5 <- read.csv("Genomic_Ranges_CSV_NA.csv", header = TRUE, stringsAsFactors = FALSE)
data5

# Clean data
data_cleaned <- na.omit(data5)
data_cleaned

# 5. C. Transformation and Normalization

# Data may require transformation (logarithmic, square root) or normalization 
# (scaling to a 0 - 1 range). This can be accomplished using base R or the "dplyr" 
# package.

# Load library
library(dplyr)

# Read in data set with numeric data
data6 <- read.csv("Genomic_Ranges_CSV_numeric.csv", header = TRUE, stringsAsFactors = FALSE)
data6

# Normalization of numerical columns across the data set
data_normalized <- data6 %>% mutate(across(where(is.numeric), ~ scale(.)))
data_normalized

# 5. D. Formatting and Reshaping Data

# Biological data sets may be wide or long format. The "tidyverse" provides 
# "pivot_longer()" and "pivot_wider()" functions that facilitate this transformation.

# Load Library
library(tidyr)

# Read in data set 
data6 <- read.csv("Genomic_Ranges_CSV_numeric_ID.csv", header = TRUE, stringsAsFactors = FALSE)
data6

# Conversion of wide data set to a long data set
data_long <- pivot_longer(data6, cols = c("start", "end"), 
                          names_to = "variable", values_to = "value")
data_long

# Conversion of a long data set back to a wide data set
data_wide <- pivot_wider(data_long, names_from = variable,
                         values_from = value)
data_wide

# 6. Managing Large Data Sets

# When working with large biological data sets, memory limitations can become an 
# issue. R offers several strategies:

# 6. A. Data Table

# The "data.table" package provides an enhanced version of data frames that is 
# more efficient in terms of memory and speed.

# Load Library
library(data.table)

# Read in data set 
data7 <- fread("Literature search - complete set.csv")
data7

# 6. B. Database Connections

# For extremely large data sets, consider using R to connect to databases (e.g. 
# SQLite, MySQL) allowing for efficient data management without loading entire 
# data sets into memory.

# Load Library
# library(DBI)

# con <- dbConnect(RSQLite::SQLite(), "my_database.sqlite")
# data <- dbReadTable(con, "table_name")

# 7. Best Practices in Data Management

# To ensure the integrity and reproducibility of analyses, adhere to the following 
# best practices:

# __Document Everything__: Keep a detailed log of all data manipulations and 
# analyses, either with comment in the R script or in a separate README file.

# __Regular Backups__: Store data sets and scripts in version-controlled systems
# (e.g. Git) to track changes and recover older versions as necessary.

# __Automate Processes__: Automate repetitive tasks using functions or scripts 
# to minimize human error and save time.

# Importing and managing biological data sets in R requires a comprehensive 
# understanding of data types, the appropriate use of packages and functions, 
# and meticulous data management practices. By mastering these techniques, 
# researchers can significantly enhance the robustness and reliability of their 
# analyses, paving the way for insightful discoveries in the biological sciences

