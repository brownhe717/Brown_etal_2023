#GO analysis of genes in each list
#enriched comparison
library(org.Dm.eg.db)
library(clusterProfiler)

E_vgPc_sdPc_Pc = read.delim("INPUT/vgPc_sdPc_Pc_overlapE_list.txt")
#View(E_vgPc_sdPc_Pc)

E_geneList = E_vgPc_sdPc_Pc[,1]
#head(E_geneList)

E_ego = enrichGO(E_geneList, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
head(E_ego)

dotplot(E_ego, showCategory=20)
