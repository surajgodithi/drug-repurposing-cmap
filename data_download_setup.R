library(GEOquery)
library(tidyverse)
library(recount3)

setwd("F:/drug_repurposing")
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

print("Downloading supplementary files for GSE142279...")

# Download all supplementary files
getGEOSuppFiles("GSE142279", makeDirectory = TRUE, baseDir = "data/")

# List what files we got
supp_files <- list.files("data/GSE142279/", full.names = TRUE)
print("Downloaded files:")
print(supp_files)

# Look for common RNA-seq file patterns
count_files <- supp_files[grepl("count|matrix|expression|fpkm|tpm|norm", 
                                supp_files, ignore.case = TRUE)]
print("Potential expression data files:")
print(count_files)

library(R.utils)

# Decompress the .gz file
gunzip("data/GSE142279/GSE142279_FPKM.xls.gz", 
       destname = "data/GSE142279/GSE142279_FPKM.xls", 
       overwrite = TRUE)

# Now read the Excel file
library(readr)
expression_data <- read_delim("F:/drug_repurposing/data/GSE142279/GSE142279_FPKM.xls", 
                              delim = "\t")

# Check what we got
print("Data loaded successfully!")
print(dim(expression_data))
print(head(expression_data[,1:5]))


# Look at all column names to confirm the pattern
print("All sample names:")
print(colnames(expression_data)[-1])  # Exclude the ID column

# Check the sample pattern
sample_names <- colnames(expression_data)[-1]
conditions <- ifelse(grepl("_T$", sample_names), "Tumor", "Normal")
print("Sample conditions:")
print(table(conditions))

# Look at a few more rows to understand the data
print("Sample of expression data:")
print(expression_data[1:10, 1:6])

# Check for missing values
print("Any missing values?")
print(sum(is.na(expression_data)))
