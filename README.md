# SPT5_erythropoiesis

D. Martell
April 10, 2020

### For RNA-seq bam files start with the STAR alignment and QC file:

STAR_and_QC.txt

### For NET-seq bam files go to the following gitHub repo for the alignment and QC files:

https://github.com/churchmanlab/MiMB2019NETseq

### Next use featureCounts to get gene counts, use the appropriate annotation file (for example, NET-seq pause index requires an annotation file for the promoter region and gene body):

FeatureCounts.rtf

### Perform differential expression analysis:

DESeq2_script.Rmd

### Perform functional analysis to get enriched GO terms:

functional_analysis.Rmd

### Further NET-seq data analysis, metagene and pausing index:

PI_Jan2023.Rmd
Metagene_Plotting_Final.R
