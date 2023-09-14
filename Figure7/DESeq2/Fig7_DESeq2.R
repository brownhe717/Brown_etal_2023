#!/usr/env/bin Rscript

library(DESeq2)
library(tidyverse)

#Note: scriptflow: DESeq2 --> MA & GO (supplemental Figure 5) --> Venn Diagram & eGO (Figure 7)

# Read count matrix
counts <- read.delim("INPUT/hyp_counts.txt",sep="\t",row.names="Geneid")

# Generate experimental design formula (rewritten to generate rownames in the same command)
coldata <- data.frame(row.names = colnames(counts), "tissue" = c(rep("w1118",3),
                                                                 rep("PcKD",3),
                                                                 rep("vgLOF_PcKD",3),
                                                                 rep("sdLOF_PcKD",3)))

# Input count matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~tissue)

# Filter out rows with less than 10 counts total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Perform Wald test
dds <- DESeq(dds)

#######################
### WT EAD vs. PcKD ###
#######################

# Perform DE testing 
res_w1118_vs_PcKD <- results(dds, contrast=c("tissue","PcKD", "w1118"))

# Add fold change classification to complete results
res_w1118_vs_PcKD <- as.data.frame(res_w1118_vs_PcKD) %>%
  rownames_to_column("geneid") %>%
  mutate(w1118_vs_PcKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in absence of Pc",
    log2FoldChange <= -1 ~ "depleted in absence of Pc",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_w1118_vs_PcKD_sig <- subset(res_w1118_vs_PcKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_w1118_vs_PcKD_sig, "w1118.v.PcKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

#############################
### WT EAD vs. vgLOF_PcKD ###
#############################

# Perform DE testing 
res_w1118_vs_vgPc <- results(dds, contrast=c("tissue","vgLOF_PcKD", "w1118"))

# Add fold change classification to complete results

res_w1118_vs_vgPc <- as.data.frame(res_w1118_vs_vgPc) %>%
  rownames_to_column("geneid") %>%
  mutate(w1118_vs_vgPc_change = case_when(
    log2FoldChange >= 1 ~ "enriched in absence of Pc and vg",
    log2FoldChange <= -1 ~ "depleted in absence of Pc and vg",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_w1118_vs_vgPc_sig <- subset(res_w1118_vs_vgPc, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_w1118_vs_vgPc_sig, "w1118.v.vgLOF_PcKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

############################
### w1118 vs. sdLOF_PcKD ###
############################

# Perform DE testing 
res_w1118_vs_sdPc <- results(dds, contrast=c("tissue","sdLOF_PcKD","w1118"))

# Add fold change classification to complete results
res_w1118_vs_sdPc <- as.data.frame(res_w1118_vs_sdPc) %>%
  rownames_to_column("geneid") %>%
  mutate(w1118_vs_sdPc_change = case_when(
    log2FoldChange >= 1 ~ "enriched in absence of sd and Pc",
    log2FoldChange <= -1 ~ "depleted in absence of sd and Pc",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_w1118_vs_sdPc_sig <- subset(res_w1118_vs_sdPc, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_w1118_vs_sdPc_sig, "w1118.v.sdLOF_PcKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

###########################
### PcKD vs. vgLOF_PcKD ###
###########################

# Perform DE testing 
res_PcKD_v_vgPc <- results(dds, contrast=c("tissue","vgLOF_PcKD", "PcKD"))

# Add fold change classification to complete results
res_PcKD_v_vgPc <- as.data.frame(res_PcKD_v_vgPc) %>%
  rownames_to_column("geneid") %>%
  mutate(PcKD_v_vgPc_change = case_when(
    log2FoldChange >= 1 ~ "enriched in absence of vg",
    log2FoldChange <= -1 ~ "depleted in absence of vg",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_PcKD_v_vgPc_sig <- subset(res_PcKD_v_vgPc, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_PcKD_v_vgPc_sig, "Pc.v.vgPc.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

###########################
### PcKD vs. sdLOF_PcKD ###
###########################

# Perform DE testing 
res_PcKD_v_sdPc <- results(dds, contrast=c("tissue","sdLOF_PcKD", "PcKD"))

# Add fold change classification to complete results
res_PcKD_v_sdPc <- as.data.frame(res_PcKD_v_sdPc) %>%
  rownames_to_column("geneid") %>%
  mutate(PcKD_v_sdPc_change = case_when(
    log2FoldChange >= 1 ~ "enriched in absence of sd",
    log2FoldChange <= -1 ~ "depleted in absence of sd",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_PcKD_v_sdPc_sig <- subset(res_PcKD_v_sdPc, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_PcKD_v_sdPc_sig, "Pc.v.sdPc.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)
