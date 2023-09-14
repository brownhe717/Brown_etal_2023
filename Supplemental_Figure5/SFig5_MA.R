#!/usr/env/bin Rscript

library(ggpubr)
library(viridis)

#input = ReferenceDataTracks from figure 7 DESeq2 script

####################
##### MA plots #####
####################

# WT eye-antennal disc (EAD) v. sd1; ey>Pc RNAi (sdPcKD)
ggmaplot(res_w1118_vs_sdPc, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_w1118_vs_sdPc$geneid), top = 8, 
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))


# WT EAD v. vg1; DE>Pc RNAi (vgPcKD)
ggmaplot(res_w1118_vs_vgPc, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_w1118_vs_vgPc$geneid), top = 8, 
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))


# PcKD EAD v. sdPcKD
ggmaplot(res_PcKD_v_sdPc, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_PcKD_v_sdPc$geneid), top = 8, 
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# PcKD EAD v. vgPcKD
ggmaplot(res_PcKD_v_vgPc, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_PcKD_v_vgPc$geneid), top = 8, 
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))
