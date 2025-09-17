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
- **limma** — linear modeling for differential expression  
- **edgeR** — count-based RNA-seq DE with negative binomial models  
- **GEOquery** — download/parse GEO series/platforms  
- **recount3** — access uniformly processed RNA-seq data  

### CRAN
- **tidyverse** — dplyr, ggplot2, readr, tibble, etc.  
- **pheatmap** — heatmaps  
- **R.utils** — utilities (e.g., file ops, gunzip helpers)  

### Visualization
- **ggplot2**, **EnhancedVolcano**, **pheatmap**, **ComplexHeatmap**  

---

## 📊 Differential Expression Results
- **2,882 significant DEGs** detected (adj. p < 0.05)  
  - **Upregulated genes (n=1,306):** higher expression in tumor  
  - **Downregulated genes (n=1,576):** lower expression in tumor  
- Volcano plot confirms clear separation of tumor vs normal expression profiles.  
- Gene lists prepared for CMap submission (top 150 up + 150 down).  

---

## 🧪 Detailed CMap Results Analysis

### Top 10 Drug Candidates (Connectivity Scores)

Based on our colorectal cancer gene signature (139 upregulated, 143 downregulated genes), CMap identified the following compounds with the strongest negative connectivity scores (indicating potential to reverse the tumor signature):

| Rank | Compound | Mechanism of Action | τ Score | Clinical Status |
|------|----------|-------------------|---------|----------------|
| 1 | PKC activator (class) | Protein kinase C activation | -98.45 | Multiple compounds in development |
| 2 | Prostratin | PKC activator | -97.36 | Research compound |
| 3 | Avrainvillamide-analog-3 | Nucleophosmin inhibitor | -97.43 | Experimental |
| 4 | Ingenol | PKC activator | -97.37 | **FDA-approved** (ingenol mebutate) |
| 5 | Phorbol-12-myristate-13-acetate | PKC activator | -96.58 | Research tool compound |
| 6 | ON-01910 (Rigosertib) | PLK inhibitor | -96.43 | **Clinical trials** (Phase III) |
| 7 | Scoulerine | Adrenergic receptor antagonist | -96.21 | Natural alkaloid |
| 8 | Vinorelbine | Tubulin inhibitor | -95.17 | **FDA-approved** chemotherapy |
| 9 | MLN-8054 | Aurora kinase inhibitor | -94.80 | Clinical development discontinued |
| 10 | SB-216763 | Glycogen synthase kinase inhibitor | -94.68 | Research compound |

### Key Findings

**PKC Activators Dominate**: Four of the top 10 hits (ranks 1, 2, 4, 5) target protein kinase C signaling, suggesting PKC pathway modulation strongly opposes colorectal cancer gene expression patterns.

**Clinically Relevant Hits**: 
- **Vinorelbine** (τ = -95.17) is an FDA-approved chemotherapy agent already used in various cancers
- **Ingenol mebutate** (τ = -97.37) is FDA-approved for actinic keratosis treatment
- **Rigosertib** (τ = -96.43) reached Phase III trials for myelodysplastic syndromes

**Novel Targets**: Nucleophosmin inhibition (Avrainvillamide-analog-3) and GSK inhibition (SB-216763) represent less explored pathways in colorectal cancer.

### Statistical Context
- **Total compounds screened**: ~1,300 in CMap Touchstone database
- **Significance threshold**: τ > |90| considered strong connectivity
- **Hit rate**: 10/1,300 compounds (0.77%) showed strong negative connectivity
- **Gene coverage**: 92% of submitted upregulated genes, 96% of downregulated genes were recognized by CMap

### Biological Interpretation

The strong enrichment of **PKC activators** suggests that colorectal tumors may have suppressed PKC signaling pathways. PKC proteins regulate:
- Cell differentiation and apoptosis
- Tumor suppressor functions
- Cell cycle checkpoints

**Mechanistic hypothesis**: Colorectal cancers may evade growth control by downregulating PKC-dependent tumor suppressor mechanisms. PKC activators could potentially restore these protective pathways.

### Clinical Implications

**Immediate repurposing candidates**:
1. **Vinorelbine** - Already approved chemotherapy, could be tested in colorectal cancer models
2. **Ingenol derivatives** - Topical agent that might be reformulated for systemic use

**Drug development priorities**:
1. **PKC modulator optimization** - Multiple hits suggest this pathway merits focused drug development
2. **Rigosertib evaluation** - Failed compound in other cancers may have colorectal cancer activity

### Limitations

**Computational predictions require validation**:
- CMap uses immortalized cell lines, not primary colorectal tissue
- Connectivity scores don't predict therapeutic windows or toxicity
- Gene expression reversal may not translate to tumor growth inhibition
- Patient tumor heterogeneity not captured in this analysis

**Next steps for validation**:
1. Test top compounds in colorectal cancer cell lines
2. Validate gene expression changes in vitro
3. Assess growth inhibition and apoptosis induction
4. Evaluate compounds in patient-derived xenograft models

### Conclusion

This analysis identified **PKC pathway modulation** as a priority therapeutic strategy for colorectal cancer, with **vinorelbine** and **ingenol derivatives** representing immediate repurposing opportunities. The computational predictions provide a focused starting point for experimental validation rather than definitive therapeutic recommendations.

---

## 🔮 Roadmap
- [x] Download and preprocess GEO dataset  
- [x] Run DEG analysis (tumor vs normal)  
- [x] Generate volcano plot + DEG summary  
- [x] Submit signature to CMap (Touchstone 1.0)  
- [x] Summarize top compounds and MoAs  
- [ ] Annotate FDA-approved vs experimental compounds  
- [ ] Cross-check hits with PRISM viability data  
- [ ] Add literature review for top candidates  
- [ ] Propose validation experiments  

---

## ⚙️ Reproducibility
All analysis performed in **R (v4.x)**.  

Install required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "GEOquery", "recount3"))
install.packages(c("tidyverse", "pheatmap", "R.utils"))
