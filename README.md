# deseq
Differential Expression Analysis using DESeq2 This repository contains an R script for performing differential gene expression analysis using DESeq2. The workflow includes downloading RNA-seq raw counts and metadata, preprocessing the data, running DESeq2 for differential expression analysis, and visualizing results with various plots.
# Differential Expression Analysis using DESeq2  

This repository contains an R script for performing differential gene expression analysis using DESeq2. The workflow includes downloading RNA-seq raw counts and metadata, preprocessing the data, running DESeq2 for differential expression analysis, and visualizing results with various plots.  

## Steps Covered  

### 1. Downloading Data  
- Raw counts from [EBI Gene Expression Atlas](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5243)  
- Experimental metadata  

### 2. Data Wrangling  
- Formatting count data with gene IDs as row names  
- Cleaning and restructuring metadata  
- Converting genotype information into a factor  

### 3. Differential Expression Analysis  
- Creating a `DESeqDataSet`  
- Filtering out low-count genes  
- Running DESeq to analyze differential gene expression  

### 4. Spot-checking Expression of Key Genes  
- Extracting and plotting expression levels for genes of interest (e.g., `ZEB1`)  

### 5. Visualization  
- **MA Plot:** Highlights differentially expressed genes  
- **Volcano Plot (EnhancedVolcano):** Shows significance vs. fold change of genes  

## Key Commands  

```r
# Generate an MA plot
plotMA(res)

# Create a volcano plot for significant gene expression changes
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')
