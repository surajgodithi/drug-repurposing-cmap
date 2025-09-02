# =============================================================================
# cmap_preperation.R
# Cancer Drug Repurposing - GSE142279
# Converting gene signatures for Cmap analysis
# =============================================================================

library(tidyverse)
library(biomaRt)

if (!exists("top_up")) {
  top_up_df <- read.csv("results/cmap_upregulated_genes.csv")
  top_down_df <- read.csv("results/cmap_downregulated_genes.csv")
  
  up_ensembl <- top_up_df$ensembl_id
  down_ensembl <- top_down_df$ensembl_id
} else {
  up_ensembl <- top_up$ensembl_id
  down_ensembl <- top_down$ensembl_id
}

print(paste("Converting", length(up_ensembl), "upregulated genes"))
print(paste("Converting", length(down_ensembl), "downpregulated genes"))

#Converting Enseml ID to gene symbols

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# converting upregulated genes
up_conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = up_ensembl,
  mart = ensembl
)

# converting downregulated genes
down_conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = down_ensembl,
  mart = ensembl
)

print("Conversion results:")
print(paste("Upregulated: ", nrow(up_conversion), "genes converted"))
print(paste("Downregulated: ", nrow(down_conversion), "genes converted"))

#Cleaning

up_symbols <- up_conversion$hgnc_symbol[up_conversion$hgnc_symbol != ""]
down_symbols <- down_conversion$hgnc_symbol[down_conversion$hgnc_symbol != ""]

print("Final gene lists for CMap:")
print(paste("Upregulated symbols:", length(up_symbols)))
print(paste("Downregulated symbols:", length(down_symbols)))

print("Sample upregulated genes:")
print(head(up_symbols, 20))
print("Sample downregulated genes:")
print(head(down_symbols, 20))

#Save for Cmap

write.table(up_symbols, "results/cmap_up_gene_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(down_symbols, "results/cmap_down_gene_symbols.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.csv(up_conversion, "results/upregulated_gene_conversion.csv", row.names = FALSE)
write.csv(down_conversion, "results/downregulated_gene_conversion.csv", row.names = FALSE)
