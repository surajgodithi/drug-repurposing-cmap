# Drug Repurposing with Connectivity Map (CMap)

This project explores **drug repurposing** using gene expression data.  
The goal is to identify **differentially expressed genes (DEGs)** between **tumor vs normal samples** and then use the **Connectivity Map (CMap / CLUE)** to find drugs that could potentially *reverse* the cancer gene signature.  

---

## 🚀 Project Overview
1. **Dataset:**  
   - Source: GEO dataset [GSE142279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142279)  
   - Samples: 40 colorectal tissue samples (20 tumor, 20 matched normal)  
   - Perfect paired design → tumor and normal from the same patient  

2. **Analysis:**  
   - Performed differential expression analysis (tumor vs normal)  
   - Identified **2,882 significant DEGs**  
     - **1,306 upregulated in tumor**  
     - **1,576 downregulated in tumor**  
   - Results highly significant (adjusted p-values as low as 10^-16)  

3. **Visualization:**  
   - Volcano plot shows clear separation of upregulated (red) and downregulated (blue) genes:  

   ![Volcano Plot](figures/volcano_plot.png)

4. **Next Steps:**  
   - Extract **top 150 upregulated** and **top 150 downregulated** genes  
   - Submit to **Connectivity Map (CMap/CLUE)**  
   - Identify candidate drugs that may reverse the tumor gene signature  

---

## 📦 Tools & Packages

### Bioconductor
- limma — linear modeling for differential expression
- edgeR — count-based RNA-seq DE with negative binomial models
- GEOquery — download/parse GEO series/platforms
- recount3 — access uniformly processed RNA-seq data

### CRAN
- tidyverse — dplyr, ggplot2, readr, tibble, etc.
- pheatmap — heatmaps
- R.utils — utilities (e.g., file ops, gunzip helpers)

### Visualization
- `ggplot2`, `EnhancedVolcano`, `pheatmap`, `ComplexHeatmap`  

---

## 🧪 Results: CMap Analysis

After identifying **2,882 significant differentially expressed genes (DEGs)** in colorectal tumor vs normal samples (1,306 upregulated, 1,576 downregulated), we submitted the **top 150 upregulated** and **top 150 downregulated** genes to the **Connectivity Map (CMap, Touchstone 1.0)**.  

CMap compares our tumor “signature” against thousands of **drug-induced gene expression profiles** to identify compounds that could *reverse* the cancer state.

### 🔹 Key Findings
- Several drug classes showed **strong negative connectivity scores (τ ≈ –90 to –100)**, indicating they may reverse the colorectal tumor gene signature.
- The most consistent signals included:

| Drug Class              | Representative Compounds         | Mechanism of Action (MoA)                  |
|--------------------------|----------------------------------|---------------------------------------------|
| **PKC activators**       | Prostratin, Ingenol, PMA         | Activate protein kinase C signaling |
| **HDAC inhibitors**      | Vorinostat, Trichostatin A       | Epigenetic regulators (histone deacetylase inhibition) |
| **EGFR inhibitors**      | Gefitinib, Erlotinib             | Block EGFR tyrosine kinase signaling |
| **PI3K/mTOR inhibitors** | Rapamycin, PI-103                | Target survival/proliferation signaling |
| **Topoisomerase inhibitors** | Camptothecin, Etoposide     | Interfere with DNA replication |
| **Aurora kinase inhibitors** | Alisertib                   | Disrupt mitotic checkpoints |
| **SRC inhibitors**       | Dasatinib                        | Block Src family tyrosine kinases |
| **MEK inhibitors**       | Trametinib                       | Inhibit MAPK/ERK pathway |
| **Proteasome inhibitors**| MG-132                           | Prevent protein degradation |
| **HSP90 inhibitors**     | Geldanamycin                     | Destabilize oncogenic proteins |

### 🔹 Interpretation
- **PKC activators** and **HDAC inhibitors** emerged as the **strongest reversal classes**, suggesting that altering transcriptional and epigenetic states may counteract colorectal tumor programs.  
- Multiple independent **kinase inhibitors (EGFR, PI3K, MEK, SRC, Aurora)** also scored highly, consistent with the known dependence of colorectal cancer on signaling pathways controlling cell growth and survival.  
- Together, these results provide a **shortlist of candidate compounds and drug classes** for further validation in colorectal cancer models.

---

### 📊 Visual Summary
- **Volcano Plot:** Clear separation of upregulated (red) vs downregulated (blue) genes.  
- **CMap Results:** Multiple drug classes with strong negative connectivity, led by PKC activators and HDAC inhibitors.  

This pipeline demonstrates how **transcriptomic signatures** from patient-matched tumor/normal samples can be leveraged to identify **drug repurposing opportunities** using CMap.


## ⚙️ Reproducibility
All analysis is done in R (v4.x).  
Required packages can be installed with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "GEOquery", "recount3"))
install.packages(c("tidyverse", "pheatmap", "R.utils"))
