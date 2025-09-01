# =============================================================================
# differential_expression_analysis.R
# Cancer Drug Repurposing - GSE142279
# Differential Expression Analysis: Tumor vs Normal
# =============================================================================

library(limma)
library(edgeR)
library(tidyverse)
library(pheatmap)

#Preparing the data and design matrix

sample_names <- colnames(expression_data)[-1]
sample_annotation <- data.frame(
  sample_id = sample_names,
  condition = ifelse(grepl("_T$", sample_names), "Tumor", "Normal"),
  tissue_type = case_when(
    grepl("^RC", sample_names) ~ "Right_Colon",
    grepl("^LC", sample_names) ~ "Left_Colon",
    grepl("^RE", sample_names) ~ "Rectum"
  ),
  patient_id = gsub("_[NT]$", "", sample_names),
  stringsAsFactors = FALSE
)

print("Sample annotation:")
print(head(sample_annotation))
print(table(sample_annotation$condition, sample_annotation$tissue_type))

expr_matrix <- as.matrix(expression_data[, -1])
rownames(expr_matrix) <- expression_data$ID
colnames(expr_matrix) <- sample_names

print(paste("Expression matrix dimensions:", nrow(expr_matrix), "x", ncol(expr_matrix)))

#Filter low expressed genes

log_expr <- log2(expr_matrix + 0.1)

min_samples <- ceiling(0.25 * ncol(expr_matrix))
keep_genes <- rowSums(expr_matrix > 1) >= min_samples

print(paste("Genes before filtering:", nrow(expr_matrix)))
print(paste("Genes after filtering:", sum(keep_genes)))

expr_filtered <- log_expr[keep_genes, ]
genes_filtered <- expression_data$ID[keep_genes]


# Matrix for analysis

condition <- factor(sample_annotation$condition, levels =c("Normal","Tumor"))
patient <- factor(sample_annotation$patient_id)

design <- model.matrix(~ patient + condition)
print("Desiign matrix dimensions:")
print(dim(design))


# Differential experession analysis

print("Running differntial expression analysis")

fit <- lmFit(expr_filtered, design) #linear
fit2 <- eBayes(fit) #Empiracal Bayes

results <- topTable(fit2, coef = "conditionTumor", number = Inf, sort.by = "P")

print("Top 10 differentially expressed genes:")
print(head(results, 10))


# Add annotations and filtering genes

results$ensembl_id <- rownames(results)

sig_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]

print(paste("Total significant genes (padj < 0.05, |logFC| > 1):", nrow(sig_genes)))
print(paste("Upregulated in tumor:", sum(sig_genes$logFC > 1)))
print(paste("Downregulated in tumor:", sum(sig_genes$logFC < -1)))

#Preparing for Cmap

up_genes <- sig_genes[sig_genes$logFC > 1, ]
up_genes <- up_genes[order(up_genes$logFC, decreasing = TRUE), ]

down_genes <- sig_genes[sig_genes$logFC < -1, ]
down_genes <- down_genes[order(down_genes$logFC, decreasing = FALSE), ]

top_up <- head(up_genes, 150)
top_down <- head(down_genes, 150)

print("=============================================================================")
print("GENE SIGNATURES FOR CMAP:")
print(paste("Top upregulated genes:", nrow(top_up)))
print(paste("Top downregulated genes:", nrow(top_down)))
print("=============================================================================")

#Saving results

write.csv(results, "results/all_deg_results.csv", row.names = TRUE)
write.csv(sig_genes, "results/significant_genes.csv", row.names = TRUE)

write.csv(top_up[, c("ensembl_id", "logFC", "adj.P.Val")], 
          "results/cmap_upregulated_genes.csv", row.names = FALSE)
write.csv(top_down[, c("ensembl_id", "logFC", "adj.P.Val")], 
          "results/cmap_downregulated_genes.csv", row.names = FALSE)

# Save just gene IDs for CMap (they prefer simple lists)
write.table(top_up$ensembl_id, "results/cmap_up_gene_ids.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(top_down$ensembl_id, "results/cmap_down_gene_ids.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

print("Results saved!")

#Visualization 

library(ggplot2)

volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.6, size = 1) +
  geom_point(data = sig_genes, aes(color = ifelse(logFC > 0, "Up", "Down")), 
             size = 1.2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.7) +
  labs(title = "Volcano Plot: Colorectal Tumor vs Normal",
       x = "log2 Fold Change", 
       y = "-log10 Adjusted P-value",
       color = "Direction") +
  theme_minimal()

print(volcano_plot)
ggsave("figures/volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)


