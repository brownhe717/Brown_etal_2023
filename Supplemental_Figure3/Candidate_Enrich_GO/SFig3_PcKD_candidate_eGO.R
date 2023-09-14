#!/usr/env/bin Rscript

library(org.Dm.eg.db)
library(clusterProfiler)

###################################
###### Candidate GO Analysis ######
###################################


early_mid = read.csv("INPUT/early_mid_overlap_list.txt")
#View(early_mid)

geneList1 = early_mid[,1]
#head(geneList1)

ego1 = enrichGO(geneList1, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
#head(ego1)

dotplot(ego1, showCategory=20) 

############

early = read.csv("INPUT/early_targts_list.txt")
#View(early)

geneList4 = early[,1]
#head(geneList4)

ego4 = enrichGO(geneList4, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
#head(ego4)

dotplot(ego4, showCategory=20)
