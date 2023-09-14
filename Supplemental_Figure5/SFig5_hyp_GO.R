#!/usr/env/bin Rscript

library(org.Dm.eg.db)
library(clusterProfiler)

#input = ReferenceDataTracks from figure 7 DESeq2 script

###################################
########### GO Analysis ###########
###################################

###WT EAD v. sdPcKD
w1118_vs_sdPc_go <- compareCluster(
  geneid ~ w1118_vs_sdPc_change,
  data = res_w1118_vs_sdPc_sig,
  fun = "enrichGO",
  OrgDb = "org.Dm.eg.db",
  pAdjustMethod = "fdr",
  keyType = "SYMBOL",
  ont = "BP"
)

dotplot(w1118_vs_sdPc_go, showCategory = 20, split = NULL, font.size=8, title = "Transcript XP in WT vs. sdLOF_PcKD", by="count") + scale_color_viridis_c()

###WT EAD v. vgPcKD
w1118_vs_vgPc_go <- compareCluster(
  geneid ~ w1118_vs_vgPc_change,
  data = res_w1118_vs_vgPc_sig,
  fun = "enrichGO",
  OrgDb = "org.Dm.eg.db",
  pAdjustMethod = "fdr",
  keyType = "SYMBOL",
  ont = "BP"
)

dotplot(w1118_vs_vgPc_go, showCategory = 20, split = NULL, font.size=10, title = "Transcript XP in WT vs. vgLOF,PcKD", by="count") + scale_color_viridis_c()

###PcKD v. sdPcKD
PcKD_v_sdPc_go <- compareCluster(
  geneid ~ PcKD_v_sdPc_change,
  data = res_PcKD_v_sdPc_sig,
  fun = "enrichGO",
  OrgDb = "org.Dm.eg.db",
  pAdjustMethod = "fdr",
  keyType = "SYMBOL",
  ont = "BP"
)

dotplot(PcKD_v_sdPc_go, showCategory = 15, split = NULL, font.size=8, title = "Transcript XP in PcKD vs. sdLOF_PcKD", by="count") + scale_color_viridis_c()

###PcKD v. vgPcKD
PcKD_v_vgPc_go <- compareCluster(
  geneid ~ PcKD_v_vgPc_change,
  data = res_PcKD_v_vgPc_sig,
  fun = "enrichGO",
  OrgDb = "org.Dm.eg.db",
  pAdjustMethod = "fdr",
  keyType = "SYMBOL",
  ont = "BP"
)

dotplot(PcKD_v_vgPc_go, showCategory = 15, split = NULL, font.size=8, title = "Transcript XP in PcKD vs. vgLOF_PcKD", by="count") + scale_color_viridis_c()
