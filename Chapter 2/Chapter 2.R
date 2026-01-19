# R Programming for Bioinformatics

# Chapter 2 - Introduction to R for Bioinformatics

# 1. Getting Started with R - Downloading R
# [The comprehensive R Archive  Network (CRAN)](https://cran.r-project.org/)

# 2. Getting started with R Studio - Download R studio
# [RStudio website](https://www.rstudio.com/products/rstudio/download/)

# 3. Assign values to variables
my_variable <- 42
my_variable

# 4. Data Types - Several Types

# * Vectors - Basic data structure representing a sequence of elements

# * Data Frames - 2-Dimensional structure akin to a table, suitable for storing datasets

# * Lists - Collections of items that can be of different types

# * Functions - Utilize built-in functions to perform operations

mean_value <- mean(c(1, 2, 3, 4, 5)) # Calculate a mean of a vector
mean_value

# 5. Data Manipulation

# Data preprocessing is critical in bioinformatics for accurate analysis. 
# The "dplyr" package part of the "tidyverse" offers powerful tools for data 
# manipulation

# 5. A. Reading Data using R functions - CSV file
data <- read.csv("Literature search - complete set.csv")
head(data, n=1)

# 5. B. Reading Data using R functions - Excel file
install.packages("readxl")
library(readxl)
data1 <- read_excel("Literature search - complete set.xlsx")
head(data1, n=1)

# 5. C. Data Cleaning

# Involves handling missing values, filtering data sets and tranforming variables. 
# Key functions in "dplyr" include;

# * filter() - Select rows based on conditions

# * select() - Choose specific columns

# * mutate() - Create new variables or modify existing ones.

# 5. C. I. Load Libraries
library(devtools)
library(Biobase)
library(org.Hs.eg.db)
library(DBI)
library(AnnotationDbi)
library(BiocGenerics)

# 5. C. II. Load in data from a connection
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

# 5. C. III. Load Libraries
library(tidyverse)
library(dplyr)
library(tibble)

# 5. C. IV. Change expression set to a data frame
edata1 = as.data.frame(edata)
edata2 = as_tibble(edata)
edata3 = data.frame(edata)

# 5. C. V. Filter out data for NA's and Mutate data to create a new variable log expression (only first column used)
clean_data <- edata1 %>% filter(!is.na(edata1)[,1] & !is.na(edata1)[,2] &
                                !is.na(edata1)[,3] & !is.na(edata1)[,4]) %>% mutate(
                                  log_SRX033480 = log(edata1[,1] + 1)) %>% mutate(
                                  log_SRX033488 = log(edata1[,2] + 1)) %>% mutate(
                                  log_SRX033481 = log(edata1[,3] + 1)) %>% mutate(
                                  log_SRX033489 = log(edata1[,4] + 1)) 

head(clean_data, n=1)

# 6. Data Summarization
summary_data <- clean_data %>% group_by(fdata_bot@data) %>% summarize(
  mean_log_SRX033480 = mean(log_SRX033480),
  mean_log_SRX033488 = mean(log_SRX033488),
  mean_log_SRX033481 = mean(log_SRX033481),
  mean_log_SRX033489 = mean(log_SRX033489))
  

head(summary_data)

# 7. Data Visualization
library(ggplot2)

ggplot(clean_data[,22:25], aes(x = clean_data[,22], y = clean_data[,25])) + geom_point() + 
  labs(title = "Gene Expression of sample log_SRX033480 vs log_SRX033480", 
       x = "Gene Expression sample log_SRX033480", y = "Gene Expression sample log_SRX033480")

# 8. Advanced Visualization - Heat map
install.packages("pheatmap")
library(pheatmap)
pheatmap(as.matrix(clean_data[,22:25]), cluster_rows = TRUE, cluster_cols = TRUE, 
         main = "Heatmap of Gene Expressions")

# 9. Statistical tests - t - test

# 9. A. Identifying the low genes
low_genes = rowMeans(edata1) < 5
table(low_genes)

# 9. B. Filter out low gene values, counts < 5
filt_edata = filter(edata1,!low_genes)
dim(filt_edata)

# 9. C. Transpose Summary Data Set
filt_edata_t = t(filt_edata)
dim(filt_edata_t)

# 9. D. Log2 Transform
log2_filt_edata_t = log2(filt_edata_t + 1)

# 9. E. Make gene1 and gene2 two numeric vectors
gene_ENSMUSG00000000001 <- as.numeric(log2_filt_edata_t[,1])
gene_ENSMUSG00000000001 <- na.omit(gene_ENSMUSG00000000001)
gene_ENSMUSG00000000003 <- as.numeric(log2_filt_edata_t[,2])
gene_ENSMUSG00000000003 <- na.omit(gene_ENSMUSG00000000003)

# 9. F. Create a data frame
genes12 <- data.frame(
  gene = rep(c("gene_ENSMUSG00000000001", "gene_ENSMUSG00000000003"), each = 21),
  gene_expression = c(gene_ENSMUSG00000000001, gene_ENSMUSG00000000003)
)

# 9. G. Print all data
print(genes12)

# 9. H. Compute t-test
T_test_result <- t.test(gene_expression ~ gene, data = genes12, var.equal = TRUE)
T_test_result

# 10. Bioinformatics Packages

# * Bioconductor - Tools for analyzing Genomic Data

# * DESeq2 - Differential Gene Analysis

# * EdgeR - RNA-seq data analysis

# 10. A. Installation of Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE))
  installed.packages("BiocManager")BiocManager::install()

# 10 . B. Install specific bioinformatics Packages
BiocManager::install("DESeq2")

# 11. Basic Syntax, Data Types and Structures

# 11. A. R Syntax Basics

# R code is structured in expressions; variable assignment, executed as commands 
# and combined to perform complex operations

# 11. A. I. Variables and Assignment <- or =
my_variable <- 10
my_variable
another_variable = "Hello, R!"
another_variable

# 11. A. II. Comments
# This is a single-line comment
my_variable <- 10 # assigning a value 10 to my_variable

# 11. A. III. Basic Functions - Calculating the mean of a numeric vector
mean_value <- mean(c(1, 2, 3, 4, 5))

# 11. A. IV. Control Structures - "if", "for", and "while" loops - allow conditional execution

# If-else statement
if (my_variable > 5){print("Greater than 5")
  
}else{print("5 or less")
    
}

# For Loop
for(i in 1:5){print(i)
  
}

# 11. B. Data Types in R

# 11. B. I. Numeric - represents numbers - integers & real numbers

# Numeric
num1 <- 42
num1

# Numeric 
num2 <- 3.14
num2

# 11. B. II. Character - store text strings
text <- "This is a string"
text

# 11. B. III. Logical - Boolean values of TRUE and FALSE
is_valid <-TRUE
is_valid

# 11. B. IV. - Factor - represents categorical data, storing values as levels
categories <- factor(c("Low", "Medium", "High"))
categories

# 11. B. V Date and Time - handles date and time by using "Date" and "POSIXct" classes

# Current Date
today <- Sys.Date()
today

# Current Date and Time
timestamp <- Sys.time()
timestamp

# 11. C. Data Structures in R

# There are several data structures used in R for storing and manipulating data 
# effectively for efficient data analysis

# 11. C. I. Vectors - One dimensional array that can hold a sequence of elements 
# of one data type

# Numeric Vector
numeric_vector <- c(1, 2, 3, 4, 5)
numeric_vector

# Character Vector
character_vector <- c("A", "B", "C")
character_vector

# 11. C. II. Lists - hold elements of different data types e.g. vectors, data 
# frames & other lists.
my_list <- list(name = "Alice", age = 30, scores = c(90, 85, 88))
my_list

# 11. C. III. Matrices - 2D-Array with elements of the same data type

# 3X3 matrix
my_matrix <- matrix(1:9, nrow = 3)
my_matrix

# 11. C. IV. Data Frames - 2D-structure with elements of different data types e.g spreadsheets/tables
my_data_frame <- data.frame(
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 35),
  Height =c(5.5, 6.0, 5.8)
)
my_data_frame

# 11. C. V. Arrays - Hold data in more than two dimensions, all entries are of 
# the same data type

# 2X2X2 array
my_array <- array(1:8, dim = c(2, 2, 2))
my_array

# 12. Writing Functions and Scripts for Bioinformatics Tasks

# Efficient reproducible and automated data analysis methods are critical in 
# bioinformatics. One of the most effective ways to handle such tasks is by writing 
# functions and scripts tailored to bioinformatics applications.

# 12. A. Understanding the Basics of Scripting

# 12. A. I. Choosing an appropriate Programming Language.

# __Python__ - advantages include readability and extensive libraries for 
# scientific computing

# __R__ - strong in statistical analysis and visualization and numerous packages 
# tailored for biostatistics and genomics.

# __Perl__ - strong in text processing and rapid prototyping

# __Bash/Shell Scripting__ - Useful in automating command-line tasks and 
# processing large data sets

# 12. A. II. Basic Concepts

# Strong understanding of basic programming concepts like variables, data types, 
# control structures (if statements, loops) and functions

# 12. B. Writing Functions

# Functions are reusable pieces of code that perform a specific task and allows 
# development of modular code which makes maintenance, reading and debugging of 
# code easy.

# 12. B. I. Defining Functions - using Python
# def gc_content(sequence):
#   g_count = sequence.count('G')
#   c_count = sequence.count('C')
#   total_count = len(sequence)
#   return (g_count + c_count) / total_count * 100
  
# 12. B. II. Function Parameters and Return Values - returning multiple values using 
# tuples
# def nucleotide_counts(sequence):
#   a_count = sequence.count('A')
#   t_count = sequence.count('T')
#   g_count = sequence.count('G')
#   c_count = sequence.count('C')
#   return a_count, t_count, g_count, c_count

# 12. B. III. Handling Exceptions - handling errors ensures robust function creation
# - using try-except blocks to manage exceptions.
# def safe_gc_content(sequence):
#   try:
#       return gc_content(sequence) 
#   except ZeroDivisionError:
#       return 0

# 13. Writing Scripts

# Scripts are collections of functions and commands that can be executed to 
# perform a bioinformatics task

# 13. A. Structuring a Script

# Good organization enhances readability. 

# * Documentation
# * Imports
# * Function definitions
# * Execution code

# 13. A. I. Documentation

# """This script calculates the GC content of a DNA sequence."""

# 13. A. II. Import necessary libraries
# install.packages("reticulate")
# library(reticulate)
# py_install("Bio")
# from Bio.Seq import Seq

# 13. A. III. Function to calculate GC content
# def gc_content(sequence):

# 13. A. IV. Function implementation
# pass

# 13. A. V. Main execution block
# if __name__=="__main__":
#   dna_seq = Seq("AGCTAGCTAGCTAGCTAGCTGCA")
# print(f"GC Content:{gc_content(dna_seq)}%")

# 13. B. Command-Line Arguments

# For scripts that accept input directly from the command line use "argparse" 
# library in Python to handle this.

# import argparse

# def main():
#   parser = argparse.ArgumentParser(description = "Calculate GC content.")
#   parser.add_argument("sequence", type = str, help = "Input DNA sequence")
#   args = parser.parse_args()
  
# print(f"GC Content: {gc_content(args.sequence)}%")
# if name__ == "__main__":
#   main()

# 14. Automating Workflows

# In bioinformatics, workflows often involve executing several scripts or tools 
# in sequence. Scripting languages can automate these tasks with the help of 
# workflow management systems.

# * Snakemake,
# * Nextflow

# A single script encapsulates complex processes and that script calls other 
# scripts or tools while preserving dependencies and execution order.

# 14. A. Integrating External Tools

# Using subprocesses allows you to integrate command-line tools directly into 
# Python scripts.

# import subprocess

# def run_blast(query, db):
#   result = subprocess.run(['blastn', '-query', query, '-db', db]), 
#   capture_output = True, text = True)

#   return result.stdout

# 14. B. Logging

# Keeping track of what has been done is crucial for reproducibility and debugging.
# Utilizing logging captures script activity

# import logging

# logging.basicConfig(level=logging.INFO) def main():
#   logging.info("Starting GC content calculation")

# more code

# 15. Best Practices

# * __Modularity__: Write small, specific functions that do one thing well.

# * __Documentation__: Comment your code and write docstrings for your functions 
# to explain their purpose and usage.

# * __Version Control__: Use systems like Git to manage changes and collaborate 
# efficiently.

# * __Testing__: Implement unit tests to ensure your functions work as expected

# Writing functions and scripts is essential for efficient data analysis in 
# bioinformatics. Understanding how to structure and implement your scripts will
# enhance your ability to analyze biological data and share your findings 
# reproducibly.





  









