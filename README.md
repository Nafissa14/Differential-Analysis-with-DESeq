# dRNA-seq_DESeq

I followed a tutorial showing how to perform a differential gene expression project with RNA-seq data.

The Data used was from the Expression Atlas database: https://www.ebi.ac.uk/gxa/home

## Steps:

The different steps of the analysis are the following:
 ### 1) Data wrangling for DESeq2:
    Which consisted in making the data ready for the differential analysis: 
    - Renaming row names of the counts and the metadata datasets
    - Removing spaces in names to avoid DESeq warnings
    - Make sure the metadata is using "factor" type

### 2) Checking the logic in the data
    Which consisted of making sure that the data makes sense: for example, the knockout version of a gene should be associated with a smaller number of counts in a given sample compared to the wild-type version

### 3) Run DESeq
    We first create a DESeqDataSet object specifying the counts, the metadata and the "model" of the analysis (any variable that could bring some change in the expression levels)


### 4) Spot Checks
    We compared our results to the GXA (Gene Expression Atlas) ones to check if we have similar results regarding the down-regulation of the up-regulation of specific genes.

### 5) Visualization
    We tried several visualization techniques for the differential analysis results: the MA plot, the Volcano plot, and the Circos Plot.
