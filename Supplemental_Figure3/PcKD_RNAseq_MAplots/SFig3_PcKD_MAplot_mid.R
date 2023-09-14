#!/usr/env/bin Rscript

library(ggpubr)
library(viridis)

#input = ReferenceDataTracks from figure 2 DESeq2 script

####################
##### MA plots #####
####################

# WT eye-antennal disc (EAD) v. wing disc (WD)
ggmaplot(res_ead_vs_wd_mid, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_ead_vs_wd_mid$geneid), top = 2, 
         label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))


# WT EAD v. PcKD EAD
ggmaplot(res_ead_mid, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_ead_mid$geneid), top = 2, 
         label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"), 
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))
         #where up = EAD enriched & down = PcKD enriched


# WT WD v. PcKD EAD
ggmaplot(res_wd_mid, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_ead_mid$geneid), top = 2, 
         label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))
         #where up = PcKD (EAD) enriched and down = WD enriched (lower in PCKD)

