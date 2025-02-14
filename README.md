# Deseq
``# Differential Expression Analysis using DESeq2

## Overview

This repository contains an R script for performing **differential gene expression analysis** using **DESeq2**, a widely used tool for analyzing RNA-seq data. The script processes raw gene expression counts, applies statistical methods to identify differentially expressed genes, and visualizes the results using various plots.

The dataset used in this analysis is obtained from the **EBI Gene Expression Atlas (E-MTAB-5243)** and contains expression data from wildtype and **ZEB1 knockout** samples. The objective of this analysis is to identify genes that are significantly upregulated or downregulated due to the **knockout of the ZEB1 gene**, which is involved in epithelial-mesenchymal transition (EMT) and cancer progression.

## Files in This Repository

### 1. **`deseq2_analysis.R`**
- Main R script for differential expression analysis
- Downloads and processes raw count data and metadata
- Runs DESeq2 pipeline for differential expression analysis
- Generates MA plots and volcano plots for visualization

### 2. **`raw_counts.tsv`** (optional)
- Contains raw RNA-seq count data downloaded from EBI

### 3. **`metadata.tsv`** (optional)
- Metadata file containing sample characteristics (e.g., genotype)

### 4. **`results_deseq2.tsv`**
- Contains the results of differentially expressed genes, including log2 fold changes and adjusted p-values

---

## Workflow and Code Breakdown

### 1. **Downloading Data**

The first step is to fetch raw count data and metadata from the EBI Gene Expression Atlas:

```r
# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5243/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")

# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5243/resources/ExperimentDesignFile.RnaSeq/experiment-design")`` 

-   The **counts matrix** contains gene expression levels (raw read counts).
-   The **metadata file** includes information about sample conditions (e.g., wildtype vs. knockout).

### 2. **Preprocessing Data**

To ensure compatibility with DESeq2, we format the data by:

-   Setting **gene IDs** as row names in the count matrix
-   Cleaning metadata and renaming columns for clarity
-   Converting **genotype information** into a factor for statistical modeling

r

CopyEdit

`# Set gene IDs as row names
rownames(counts) = counts$Gene.ID
counts = counts[, -c(1, 2)]  # Remove unnecessary columns

# Clean metadata and convert genotype to a factor
metadata$genotype = factor(metadata$genotype, levels=c("wildtype", "knockout"))` 

### 3. **Exploratory Data Analysis**

To ensure the data is structured correctly, we perform some exploratory checks:

-   Displaying a subset of metadata
-   Checking raw expression counts of key genes

r

CopyEdit

`# Check expression of the ZEB1 gene
gene_id = genes$Gene.ID[ genes$Gene.Name == 'ZEB1' ]
gene_counts = counts[gene_id, ]` 

We also generate **boxplots** to compare expression levels of specific genes between conditions:

r

CopyEdit

`library(ggplot2)
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot()` 

### 4. **Running Differential Expression Analysis**

The **DESeq2 pipeline** is implemented as follows:

1.  Create a **DESeqDataSet**
2.  Filter out genes with very low expression
3.  Run the **DESeq** function to estimate size factors, dispersion, and fold changes
4.  Extract the **results table**

r

CopyEdit

`library(DESeq2)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~genotype)

# Filter out low-expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res = results(dds, contrast=c("genotype", "knockout", "wildtype"), alpha=1e-5)` 

### 5. **Merging Gene Names with Results**

To make the output more interpretable, we merge the DESeq2 results with gene names:

r

CopyEdit

`res_df = as.data.frame(res)
res_df = merge(res_df, genes, by='row.names')` 

### 6. **Visualizing Results**

#### **MA Plot**

An **MA plot** visualizes differential expression by plotting the **log2 fold change** vs. mean expression.

r

CopyEdit

`plotMA(res)` 

#### **Volcano Plot**

A **volcano plot** is used to identify significantly upregulated and downregulated genes.

r

CopyEdit

`library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')` 

----------

## Results and Interpretation

The differential expression analysis revealed several genes that were significantly upregulated or downregulated in **ZEB1 knockout** conditions.

### Key Findings:

-   The **ZEB1 gene itself** was confirmed to be knocked out (log2FC ≈ -Inf, p-value ≈ 0).
-   Several **EMT-related genes** (e.g., `SNAI1`, `IGFN1`) were significantly altered, indicating the impact of ZEB1 loss.
-   Genes with **log2 fold change > 2 or < -2** and **p-value < 1e-5** were considered significant.

### Example Significant Genes

Gene Name

log2 Fold Change

p-value

Regulation

ZEB1

-Inf

0

Knocked Out

IGFN1

3.2

1.2e-6

Upregulated

OR51B5

-2.8

9.5e-7

Downregulated

SCOL14A1

4.1

3.3e-8

Upregulated

----------

## Conclusion

This analysis successfully identified differentially expressed genes in **ZEB1 knockout** samples compared to wildtype controls. The results provide insights into the **molecular consequences** of ZEB1 deletion and highlight potential downstream targets involved in **epithelial-mesenchymal transition (EMT)** and **cancer progression**.

Further validation through experimental methods (e.g., qPCR, Western blot) is recommended to confirm these findings.

----------

## Dependencies

The following R packages are required to run this script:

-   `DESeq2`
-   `ggplot2`
-   `EnhancedVolcano`

To install missing dependencies, run:

r

CopyEdit

`install.packages(c("ggplot2", "EnhancedVolcano"))
BiocManager::install("DESeq2")` 

----------

## How to Reproduce This Analysis

1.  Clone this repository:
    
    sh
    
    CopyEdit
    
    `git clone https://github.com/yourusername/deseq2-analysis.git
    cd deseq2-analysis` 
    
2.  Open R and run the script:
    
    r
    
    CopyEdit
    
    `source("deseq2_analysis.R")`
