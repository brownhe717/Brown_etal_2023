#!/usr/env/bin Rscript
library(dplyr)
library(DESeq2)
library(org.Dm.eg.db)

#Note: scriptflow: DESeq2 --> MA (Supplemental Figure 3) --> ComplexHeatmap --> Filtering for candidates --> VennDiagram --> eGO (Supplemental Figure 3)

# Read count matrix
counts <- read.delim("INPUT/counts.tsv",sep="\t",row.names="Geneid")

# Generate experimental design formula (rewritten to generate rownames in the same command)
coldata <- data.frame(row.names = colnames(counts), "tissue" = c(rep("PcKD_EAD144",3),
                                                                 rep("PcKD_EAD120",3),
                                                                 rep("PcKD_EAD108",3),
                                                                 rep("PcKD_EAD96",3),
                                                                 rep("PcKD_EAD84",3),
                                                                 rep("WT_EAD120",3),
                                                                 rep("WT_EAD108",3),
                                                                 rep("WT_EAD96",3),
                                                                 rep("WT_EAD84",3),
                                                                 rep("WT_WD120",3),
                                                                 rep("WT_WD108",3),
                                                                 rep("WT_WD96",3)))

# Input count matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~tissue)

# Filter out rows with less than 10 counts total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Perform Wald test
dds <- DESeq(dds)

########################
### WT EAD vs. WT WD ###
########################

# Perform DE testing 
res_ead_vs_wd_late <- results(dds, contrast=c("tissue","WT_WD120","WT_EAD120"))

# Add fold change classification to complete results
res_ead_vs_wd_late <- as.data.frame(res_ead_vs_wd_late) %>%
  rownames_to_column("geneid") %>%
  mutate(ead_vs_wd_change = case_when(
    log2FoldChange >= 1 ~ "enriched in WD (EAD depleted)",
    log2FoldChange <= -1 ~ "depleted in WD (EAD enriched)",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_ead_vs_wd_sig_late <- subset(res_ead_vs_wd_late, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_ead_vs_wd_sig_late, "EAD.v.WD.diff.120.txt", sep="\t", col.names=T, row.names=F, quote=F)

####################
### PcKD vs. EAD ###
####################

# Perform DE testing 
res_ead_late <- results(dds, contrast=c("tissue","PcKD_EAD144", "WT_EAD120"))

# Add fold change classification to complete results
res_ead_late <- as.data.frame(res_ead_late) %>%
  rownames_to_column("geneid") %>%
  mutate(pckd_vs_ead_change = case_when(
    log2FoldChange >= 1 ~ "enriched in absence of Pc",
    log2FoldChange <= -1 ~ "depleted in absence of Pc",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_ead_sig_late <- subset(res_ead_late, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_ead_sig_late, "PcG144.v.EAD120.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

####################
### PcKD vs. WD ###
####################

# Perform DE testing 
res_wd_late <- results(dds, contrast=c("tissue","PcKD_EAD144", "WT_WD120"))

# Add fold change classification to complete results
res_wd_late <- as.data.frame(res_wd_late) %>%
  rownames_to_column("geneid") %>%
  mutate(pckd_vs_wd_change = case_when(
    log2FoldChange >= 1 ~ "higher in Pc knockdown EAD than WD",
    log2FoldChange <= -1 ~ "higher in WD than Pc knockdown EAD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_wd_sig_late <- subset(res_wd_late, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_wd_sig_late, "PcG144.v.WD120.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

