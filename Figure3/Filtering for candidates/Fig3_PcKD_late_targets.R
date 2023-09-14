library(dplyr)

changes <- left_join(res_ead_vs_wd_late, res_ead_late, by = "geneid") %>%
  left_join(., res_wd_late, by = "geneid") %>%
  dplyr::select("geneid", "ead_vs_wd_change", "pckd_vs_ead_change", "pckd_vs_wd_change")

filtered_list <- changes %>%
  dplyr::filter(ead_vs_wd_change == "enriched in WD (EAD depleted)" & pckd_vs_ead_change == "enriched in absence of Pc")

#write table of genes
write.table(filtered_list, "target_genes_late.tsv", sep="\t", row.names=T, col.names=T, quote=F)


#separate candidate list into unknown and known genes
#get list of unknown genes
CG_targets = filtered_list %>% dplyr::filter(str_detect(geneid, "CG"))
CR_targets = filtered_list %>% dplyr::filter(str_detect(geneid, "CR"))
snoRNA_targets = filtered_list %>% dplyr::filter(str_detect(geneid, "snoRNA"))

CG_CR_targets = full_join(CG_targets, CR_targets)
unknown_targets = full_join(CG_CR_targets, snoRNA_targets)
write.table(unknown_targets, "unknown_targets_late.tsv", sep="\t", col.names=T, row.names=F, quote=F)

#get list of known genes
known_targets = anti_join(filtered_list, unknown_targets)
write.table(known_targets, "known_targets_late.tsv", sep="\t", col.names=T, row.names=F, quote=F)
