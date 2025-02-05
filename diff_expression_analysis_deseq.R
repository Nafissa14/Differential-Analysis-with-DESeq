
# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)

library(DESeq2)
# ------------------------------------------------------------------------------
# Download data
# ------------------------------------------------------------------------------

# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts) # counts for each gene 

# Download metadata
# To add information about the samples 
# The experimental design file shows which sample corresponds to the gene knockout and which one is the reference 
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")

head(metadata)


# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------

# DESeq expects the counts to have gene IDs as row names( to avoid having 1,2,3,..6 as rownames)
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)


# Remove unused columns (gene ID and gene name)
# DESeq only wants to see columns that represents samples
# Let us remove the other columns 
## Keep a copy of the columns to suppress
genes = counts[, c("Gene.ID", "Gene.Name")]
genes

## Suppress the columns 1 and 2
counts = counts[, -c(1, 2)]
head(counts) # genes as rows and samples as the columns 

# DESeq expects the metadata matrix to have sample IDs in the rownames
head(metadata)
rownames(metadata)=metadata$Run
head(metadata) #Now we have the sample IDs as rows 


# Only keep columns of interest :
## Which is here the "Sample.Characteristic.genotype" which says wether each 
# sample is the wildtype or the "knock out"
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
# Look at metadata to see how the variables change with respect to each other
metadata

# Rename column
colnames(metadata) = c("genotype")
metadata

# Remove spaces in names to avoid DESeq warnings
metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
metadata$genotype[metadata$genotype == 'Snai1 knockout'] = 'knockout'
metadata

# Make sure that metadata is using "factor" type 
# Turn genotype into a factor
metadata$genotype = factor(metadata$genotype, levels=c("wildtype", "knockout"))
metadata$genotype
# We put "wildtype" as the baseline in order to compare knock out VS wild type


# ------------------------------------------------------------------------------
# Spot check expression for knockout gene SNAI1: can we see that in the row counts?
# Making sure that the data makes sense 
# ------------------------------------------------------------------------------
gene_id = genes$Gene.ID[ genes$Gene.Name == 'SNAI1' ] #consider only the SNAI1 gene
gene_counts = counts[gene_id, ]
gene_counts #we can see that some of the samples have higher expression

# Let's combine the metadata to this dataframe to clearly observe the abundance 
#of the different genotypes 
gene_data = cbind(metadata, counts=as.numeric(gene_counts))
gene_data
# The wildtype counts are in fact higher than the knockout counts : it makes sense 

# We can also plot it for better visualisation
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot()

# ------------------------------------------------------------------------------
# Run DESeq
# ------------------------------------------------------------------------------





dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~genotype)
# Run DES
dds # A DESeqDataSe

#  design=~genotype: formula that define all variables that might cause change in the analysis, variables that might cause batch effects

# Ignore genes with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Sum accross each row/gene -> counts accross all samples combined

# Run DESeq
dds <- DESeq(dds)

# Compare expression
res = results(dds, contrast=c("genotype", "knockout", "wildtype"), alpha=1e-5)
# Arg1 : name of the column in the metadata
# Adjusted p-value for multiple hypothesis testing 

res

# log2FoldChange : is something up-regulated or down-regulated

# Other examples of design formula:
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8031/Results
#     * RNA-seq on lung tumors
#     * Variables:
#         1. Location: tumor cells or normal tissue nearby
#         2. Sex of the patient
#     * Formula = ~sex + location : correct for potential batch effect due to the sex of the patients
# order does not count 
#
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10041/Results
#     * RNA-seq on mouse colon cells
#     * Variables:
#         1. Compound: cancer drug or not
#         2. Genotype: wild type mouse, and 2 other genotypes of TP53 mutations
#     * Formula = ~genotype + compound  doesn't make sense since want to compare each genotype separately (in the study: for every genotype, what happends after a change in the compound)
#     * Formula = ~group
#         metadata$group = factor(paste(metadata$genotype, metadata$compound))

# Sidenote: "~" is not a DESeq specific operator, but is used in R formula to build a model

head(iris)
model = lm(Petal.Width ~ Petal.Length, iris) # predict the width from the length from the iris dataframe
plot(iris$Petal.Length, iris$Petal.Width)
abline(model)

# ------------------------------------------------------------------------------
# Spot checks
# ------------------------------------------------------------------------------

# Compare results to GXA: 
# Do we get similar values for upregulated and down-regulated genes 
# https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Results?specific=true&geneQuery=%255B%255D&filterFactors=%257B%257D&cutoff=%257B%2522foldChange%2522%253A11%252C%2522pValue%2522%253A0.01%257D&regulation=%2522UP_DOWN%2522

# Merge gene name into data frame so can compare to GXA UI using gene names
res_df = as.data.frame(res)
head(res_df)
head(genes)
res_df = merge(res_df, genes, by='row.names')
head(res_df)
# Extract a set of genes we are interested in : 
genes_to_check = c("THY1", "SFMBT2", "PASD1", "SNAI1")
res_df[res_df$Gene.Name %in% genes_to_check, ] #Similar results

# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------

# MA plot
plotMA(res) # blue: data points that passed the threshold set of the alpha value 
# Volcano plot
# In the  value: exponent of the p vlaue -> the higher the number , the smaller the value


# Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')
# Most significant genes: in the top of the graph 

# Circos plot: more complex
# Show accross the genome: how is the gene expression: upregulation, down regulation, pvalue, where the gene names are etc..
# But gene IDs and gene names are not enough -> the software needs to know where are genes located ( chr, start and stop positions..)

BiocManager::install("biomaRt")
library(biomaRt) # Can get data from ensembl

# Find dataset name in Ensembl (Ensembl tutorial)
ensembl <- useEnsembl(biomart="genes")
datasets = listDatasets(ensembl)
head(datasets)
# Fiind the human dataset
dataset_nb = grep("human", datasets$description, ignore.case=TRUE)
dataset_nb # Human genes at row 80

dataset = datasets$dataset[dataset_nb]
dataset

# Query

ensembl <- useDataset(dataset=dataset, mart=ensembl)

# Get coordinates of all genes
attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
values <- list(ensembl_gene_id=c())
all.genes <- getBM(attributes=attributes, values=values, mart=ensembl)
head(all.genes)

# Rename column so it matches res_df in order to merge 
head(res_df)
colnames(all.genes)[1] = "Gene.ID"
head(all.genes)

# Merge the DESeq output with the Ensembl gene coordinates
merged_data <- merge(all.genes, res_df, by="Gene.ID")
head(merged_data)

# Add chr prefix to chromosome names because we have that format in the human ref genome
merged_data$chromosome_name <- paste("chr", merged_data$chromosome_name, sep = "")
head(merged_data)
# Save the dataset to a csv file 
write.csv(merged_data, "/home/nafissa/Bioinfo_Projects/dRNA-seq_DESeq/deseq.csv", row.names = FALSE)

# Same for subset: only get the rows where the gene names match the list of genes of interest
merged_data_subset = merged_data[merged_data$Gene.Name %in% genes_to_check, ]
head(merged_data_subset)
write.csv(merged_data_subset, "/home/nafissa/Bioinfo_Projects/dRNA-seq_DESeq/deseq_subset.csv", row.names = FALSE)

