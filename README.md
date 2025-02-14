# RNA-seq Differential Expression Analysis Documentation

## Table of Contents
- [Overview](#overview)
- [Data Preparation](#data-preparation)
- [DESeq2 Analysis](#deseq2-analysis)
- [Visualization](#visualization)
- [Key Commands Reference](#key-commands-reference)

## Overview
This workflow demonstrates RNA-seq differential expression analysis using DESeq2, focusing on comparing gene expression between knockout and wildtype samples.

## Data Preparation

### Download Raw Data
R
# Download count data
counts = read.delim("path/to/counts/file")

# Download metadata
metadata = read.delim("path/to/metadata/file")


### Data Wrangling
R
# Set gene IDs as row names
rownames(counts) = counts$Gene.ID

# Extract gene information
genes = counts[, c("Gene.ID", "Gene.Name")]
counts = counts[, -c(1, 2)]

# Prepare metadata
rownames(metadata) = metadata$Run
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
colnames(metadata) = c("genotype")

# Clean genotype names
metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
metadata$genotype[metadata$genotype == 'ZEB1 knockout'] = 'knockout'

# Convert to factor
metadata$genotype = factor(metadata$genotype, levels=c("wildtype", "knockout"))


## DESeq2 Analysis

### Initialize DESeq Dataset
R
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~genotype
)


### Filter and Run Analysis
R
# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)

# Get results
res = results(dds, 
    contrast = c("genotype", "knockout", "wildtype"),
    alpha = 1e-5
)


### Process Results
R
# Convert to data frame
res_df = as.data.frame(res)

# Merge with gene names
res_df = merge(res_df, genes, by='row.names')


## Visualization

### MA Plot
R
plotMA(res)


### Volcano Plot
R
library(EnhancedVolcano)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue'
)


## Key Commands Reference

| Command | Description | Example |
|---------|-------------|---------|
| read.delim() | Read tab-delimited files | counts = read.delim("file.txt") |
| DESeqDataSetFromMatrix() | Create DESeq dataset | dds <- DESeqDataSetFromMatrix(countData, colData, design) |
| DESeq() | Run differential expression analysis | dds <- DESeq(dds) |
| results() | Extract results | res = results(dds, contrast=c("group1", "group2")) |
| plotMA() | Generate MA plot | plotMA(res) |
| EnhancedVolcano() | Create volcano plot | EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue') |

### Important Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| design | Formula for experimental design | Required |
| contrast | Groups to compare | Required |
| alpha | Significance level | 0.1 |
| independentFiltering | Independent filtering | TRUE |

## Notes
- Always check data quality before analysis
- Filter low-count genes to reduce noise
- Use appropriate multiple testing correction
- Validate key findings with visualization

## Dependencies
- R (>= 4.0.0)
- DESeq2
- ggplot2
- EnhancedVolcano
