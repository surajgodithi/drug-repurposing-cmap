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

## Literature Validation of Top 3 Compounds

### 1. PKC Activators in Colorectal Cancer

**Current Evidence:**
- PKCδ-selective activators show promise in colon cancer therapy, with Roy-Bz potently inhibiting proliferation of colon cancer cells by inducing PKCδ-dependent mitochondrial apoptotic pathways
- PKCβ is the major PKC isoform expressed by colon cancer, and specific mutations can lead to increased growth rate and cell viability
- However, many clinical trials of PKC inhibitors in cancers showed no significant clinical benefits, limiting cancer therapeutic strategies targeting PKC alone

**Assessment:** PKC activators show preclinical promise but clinical translation remains challenging. The predominance of PKC activators in your CMap results aligns with research showing PKC pathway dysregulation in colorectal cancer.

### 2. Vinorelbine in Colorectal Cancer

**Current Evidence:**
- A prospective multicentre phase II clinical study of vinorelbine in BRAF V600E mutated metastatic colorectal cancer did not show signs of clinical activity despite encouraging preclinical data
- Ongoing clinical trials are testing vinorelbine specifically in advanced BRAF-like colon cancer
- Vinorelbine shows activity in gastrointestinal cancers as a tubulin inhibitor that binds tubulin and inhibits microtubule formation

**Assessment:** Vinorelbine has been tested in colorectal cancer with mixed results. While it failed in BRAF-mutated patients, ongoing trials suggest continued interest in specific molecular subtypes.

### 3. Rigosertib (PLK Inhibitor) in Colorectal Cancer

**Current Evidence:**
- Rigosertib shows potent anti-tumor responses in colorectal cancer by inhibiting RAS signaling pathway and suppressing tumor cell proliferation
- Phase I trials demonstrate rigosertib can be administered with tolerable toxicity profile and shows evidence of antitumor activity in advanced solid tumors
- Clinical progress has been hampered by lack of understanding of its mechanism of action, as it's considered a multi-target inhibitor

**Assessment:** Rigosertib shows the strongest literature support for colorectal cancer activity among your top hits, with documented efficacy against RAS-mutated colorectal cancer cells.

---

## Cross-Reference with Existing CRC Therapeutic Targets

### Established CRC Therapeutic Targets:
- **EGFR pathway:** Cetuximab, panitumumab (anti-EGFR antibodies)
- **VEGF pathway:** Bevacizumab (anti-angiogenic)
- **BRAF pathway:** Encorafenib + cetuximab combinations
- **MSI-H tumors:** Immune checkpoint inhibitors
- **HER2 amplification:** Trastuzumab + pertuzumab
- **KRAS G12C:** Sotorasib, adagrasib

### Alignment with Your Findings:
**Strong alignment:**
- **Rigosertib** directly targets RAS signaling, complementing approved KRAS inhibitors
- **Aurora kinase inhibitors** (MLN-8054) target mitotic pathways relevant to rapidly dividing cancer cells

**Novel pathways identified:**
- **PKC modulation** represents an underexplored therapeutic avenue in CRC
- **Glycogen synthase kinase inhibition** (SB-216763) may affect Wnt/β-catenin signaling, a key CRC pathway

**Limited alignment:**
- **Vinorelbine** (tubulin inhibitor) doesn't target established CRC molecular drivers
- **Nucleophosmin inhibitors** have limited precedent in CRC therapy

---

## Study Limitations and Future Directions

### Methodological Limitations

**1. Cell Line vs. Patient Tissue Discrepancy**
- CMap uses immortalized cancer cell lines that may not recapitulate primary tumor biology
- Gene expression patterns in cell culture may differ from patient tumors
- Drug responses in vitro often don't translate to clinical efficacy

**2. Gene Expression vs. Functional Outcomes**
- Reversing gene expression signatures doesn't guarantee tumor growth inhibition
- Drug connectivity scores don't predict therapeutic windows or toxicity
- Multiple resistance mechanisms may override gene expression changes

**3. Database and Technical Limitations**
- CMap Touchstone 1.0 represents a subset of available compounds (~1,300)
- Gene symbol conversion led to ~8% gene loss (282→291 genes analyzed)
- Single time-point analysis doesn't capture dynamic drug responses

**4. Biological Context Limitations**
- Colorectal cancer molecular subtypes (MSI vs MSS, KRAS vs BRAF mutations) not considered
- Tumor microenvironment effects not captured in cell line data
- Patient-to-patient heterogeneity not addressed

### Future Experimental Directions

**Immediate Validation Steps:**
1. **Test top 3 compounds** (rigosertib, vinorelbine, PKC activators) in colorectal cancer cell lines (HCT116, SW480, LS174T)
2. **Validate gene expression changes** using qPCR for top 20 upregulated/downregulated genes
3. **Assess growth inhibition** using MTT/CellTiter-Glo viability assays
4. **Confirm mechanism of action** through pathway-specific assays

**Advanced Validation:**
1. **Patient-derived organoid models** to better recapitulate primary tumor biology
2. **Combination therapy screening** since single-agent PKC targeting has shown limited clinical success
3. **Molecular subtype-specific testing** (KRAS-mutant vs BRAF-mutant vs microsatellite unstable)
4. **In vivo xenograft validation** for most promising compounds

**Clinical Translation Considerations:**
1. **Biomarker development** to identify patients most likely to respond
2. **Optimal dosing and scheduling** based on pharmacokinetic studies
3. **Combination strategies** with established CRC therapeutics
4. **Resistance mechanism studies** to predict and overcome treatment failure

### Conclusion

This computational drug repurposing analysis successfully identified biologically plausible therapeutic candidates for colorectal cancer, with **rigosertib showing the strongest literature support**. The predominance of PKC activators suggests an underexplored therapeutic pathway worthy of further investigation. However, the study limitations emphasize that these findings represent promising starting points for experimental validation rather than clinical recommendations.

The methodology demonstrates the power of connectivity mapping for hypothesis generation while highlighting the critical need for experimental validation in appropriate disease models.

---

## ⚙️ Reproducibility
All analysis performed in **R (v4.x)**.  

Install required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "GEOquery", "recount3"))
install.packages(c("tidyverse", "pheatmap", "R.utils"))
