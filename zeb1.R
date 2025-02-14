# Download data
# ------------------------------------------------------------------------------

# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5243/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)

# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5243/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)


# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------

# DESeq expects the counts to have gene IDs as row names
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)

# Remove unused columns (gene ID and gene name)
genes = counts[, c("Gene.ID", "Gene.Name")]
counts = counts[, -c(1, 2)]
head(counts)

# DESeq expects the metadata matrix to have sample IDs in the rownames
head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

# Only keep columns of interest
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
# Look at metadata to see how the variables change with respect to each other
metadata

# Rename column
colnames(metadata) = c("genotype")
metadata

# Remove spaces in names to avoid DESeq warnings
metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
metadata$genotype[metadata$genotype == 'ZEB1 knockout'] = 'knockout'
metadata

# Turn genotype into a factor
metadata$genotype = factor(metadata$genotype, levels=c("wildtype", "knockout"))
metadata$genotype


# Spot check expression for knockout gene SNAI1
# ------------------------------------------------------------------------------

gene_id = genes$Gene.ID[ genes$Gene.Name == 'ZEB1' ]
gene_counts = counts[gene_id, ]
gene_counts

gene_data = cbind(metadata, counts=as.numeric(gene_counts))
gene_data

library(ggplot2)
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot()


library(DESeq2)

# Run DESeq
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~genotype)

# Ignore genes with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)

# Compare expression
res = results(dds, contrast=c("genotype", "knockout", "wildtype"), alpha=1e-5)
res
# Sidenote: "~" is not a DESeq specific operator
head(iris)
model = lm(Petal.Width ~ Petal.Length, iris)
plot(iris$Petal.Length, iris$Petal.Width)
abline(model)

# Merge gene name into data frame so can compare to GXA UI using gene names
res_df = as.data.frame(res)
head(res_df)
head(genes)
res_df = merge(res_df, genes, by='row.names')
head(res_df)

genes_to_check = c("IGFN1", "ZEB1", "OR51B5", "SCOL14A1")
res_df[res_df$Gene.Name %in% genes_to_check, ]

# Visualization
# ------------------------------------------------------------------------------

# MA plot
plotMA(res)
# Volcano plot

library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')
