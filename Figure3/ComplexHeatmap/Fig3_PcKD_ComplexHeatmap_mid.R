#!/usr/env/bin Rscript
library(dplyr)
library(ComplexHeatmap)
library(circlize)

####################
##### Heatmaps #####
####################

# Get Gene IDs for significantly changed genes in the PcKD vs. EAD and WD comparisons, remove duplicates, 
# and use them to extract their corresponding log2FCs

ead_for_heatmap <- res_ead_sig_mid %>%
  dplyr::select("geneid")

wd_for_heatmap <- res_wd_sig_mid %>%
  dplyr::select("geneid")

heatmap_genes <- bind_rows(ead_for_heatmap, wd_for_heatmap) %>% 
  distinct

ead_log2fc <- filter(res_ead_mid, geneid %in% heatmap_genes$geneid) %>%
  dplyr::select("geneid", "log2FoldChange")

wd_log2fc <- filter(res_wd_mid, geneid %in% heatmap_genes$geneid) %>%
  dplyr::select("geneid", "log2FoldChange")

sig_matrix <- left_join(ead_log2fc, wd_log2fc, by = "geneid", suffix = c("EAD", "WD")) %>%
  column_to_rownames("geneid") %>%
  as.matrix

# Replace all rownames except those corresponding to genes of interest
rownames(sig_matrix)[!rownames(sig_matrix) %in% c("eya", "toy", "dac", "vg", "Ubx", "sd", "Scr", "Abd-B", "Antp", "ey")] <- ""

#Plot heatmap
Heatmap(sig_matrix, col = colorRamp2(c(-10, -5, 0, 5, 10), c("midnightblue", "cadetblue3", "snow1", "chocolate2", "firebrick")),
        cluster_rows = T, cluster_columns = T, row_gap = unit(2, "mm"), heatmap_legend_param = list(title = "log2(FC)"),row_dend_reorder = T,
        row_names_gp = gpar(fontsize = 8, fontface = 3), border = T, column_names_rot = 360, column_names_centered = T, show_column_names = T,
        column_labels = c("108hr WT eye-antennal disc", "108hr WT wing disc"), column_names_side = "top", column_title_side = "top", 
        column_title = "120hr Pc knockdown eye-antennal disc expression vs.")

