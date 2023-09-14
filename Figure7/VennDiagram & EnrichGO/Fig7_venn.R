#!/usr/env/bin Rscript
library(dplyr)
library(tidyverse)
library(VennDiagram)
library(RColorBrewer)

#read in tsv files of sig genes & store each as a variable (list#), convert to a data frame (list#_df), then store each list of genes as its own variable (sig_list#)
list1 = read.delim("~/Desktop/Brown et al 2022/code/Figure 7/Venn/INPUT/w1118.v.vgLOF_PcKD.diff.txt",sep="\t",row.names="geneid")
list1_df = as.data.frame(list1)

sig_list1 = rownames_to_column(list1_df, "geneid") %>% dplyr::select("geneid")
count(sig_list1)

sig_list1_enriched = rownames_to_column(list1_df, "geneid") %>% dplyr::select("geneid", "w1118_vs_vgPc_change") %>%
  filter(w1118_vs_vgPc_change == "enriched in absence of Pc and vg") #%>%
  head(sig_list1_enriched)

sig_list1_depleted = rownames_to_column(list1_df, "geneid") %>% dplyr::select("geneid", "w1118_vs_vgPc_change") %>%
    filter(w1118_vs_vgPc_change == "depleted in absence of Pc and vg") #%>%
  head(sig_list1_depleted)

count(sig_list1)
count(sig_list1_enriched)
count(sig_list1_depleted)  
  
##########    
list2 = read.delim("~/Desktop/Brown et al 2022/code/Figure 7/Venn/INPUT/wt.v.sdLOF_PcKD.diff.txt",sep="\t",row.names="geneid")
list2_df = as.data.frame(list2)

sig_list2 = rownames_to_column(list2_df, "geneid") %>%
  dplyr::select("geneid")

sig_list2_enriched = rownames_to_column(list2_df, "geneid") %>% dplyr::select("geneid", "w1118_vs_sdPc_change") %>%
  filter(w1118_vs_sdPc_change == "enriched in absence of sd and Pc") #%>%
  head(sig_list2_enriched)

sig_list2_depleted = rownames_to_column(list2_df, "geneid") %>% dplyr::select("geneid", "w1118_vs_sdPc_change") %>%
    filter(w1118_vs_sdPc_change == "depleted in absence of sd and Pc") #%>%
    head(sig_list2_enriched)

count(sig_list2)
count(sig_list2_enriched)
count(sig_list2_depleted)

##########
list3 = read.delim("~/Desktop/Brown et al 2022/code/Figure 7/Venn/INPUT/w1118.v.PcKD.diff.txt",sep="\t",row.names="geneid")
list3_df = as.data.frame(list3)

sig_list3 = rownames_to_column(list3_df, "geneid") %>%
  dplyr::select("geneid")

sig_list3_enriched = rownames_to_column(list3_df, "geneid") %>% dplyr::select("geneid", "w1118_vs_PcKD_change") %>%
  filter(w1118_vs_PcKD_change == "enriched in absence of Pc") #%>%
  head(sig_list3_enriched)

sig_list3_depleted = rownames_to_column(list3_df, "geneid") %>% dplyr::select("geneid", "w1118_vs_PcKD_change") %>%
    filter(w1118_vs_PcKD_change == "depleted in absence of Pc") #%>%
    head(sig_list3_depleted)

count(sig_list3)
count(sig_list3_enriched)
count(sig_list3_depleted)

#######combine the two lists [bind] and remove the duplicate gene ids [distinct]
###sdPc(2) v. vgPc(1)
sig_list_hyp <- bind_rows(sig_list1, sig_list2)
  dplyr::distinct(sig_list_hyp)
  count(sig_list_hyp)

sig_list_hyp <- bind_rows(sig_list1, sig_list2) %>%
  distinct
head(sig_list_hyp)
count(sig_list_hyp)

#view how many genes are in the final list
count(sig_list_hyp)

#write table of sig genes upregulated in both the WD and PcKD EAD
write.table(sig_list_hyp, "sdPc_and_vgPc_bind.txt", sep="\t", row.names=T, col.names=T, quote=F)

###sdPc (2) v. PcKD (3)
sig_list_sd <- bind_rows(sig_list2, sig_list3) %>%
  distinct
head(sig_list_sd)

#view how many genes are in the final list
count(sig_list_sd)

#write table of sig genes upregulated in both the WD and PcKD EAD
write.table(sig_list_sd, "sdPc_and_PcKD_bind.txt", sep="\t", row.names=T, col.names=T, quote=F)


###vgPc (1) v. PcKD (3)
sig_list_vg <- bind_rows(sig_list1, sig_list3) %>%
  distinct
head(sig_list_vg)

#view how many genes are in the final list
count(sig_list_vg)

#write table of sig genes upregulated in both the WD and PcKD EAD
write.table(sig_list_vg, "vgPc_and_PcKD_bind.txt", sep="\t", row.names=T, col.names=T, quote=F)

###vgPc (1), sdPc (2), and PcKD (3)
sig_list_all = bind_rows(sig_list1, sig_list2, sig_list3) %>%
  distinct
  count(sig_list_all)

  
###Venn Diagram attempts
myCol <- brewer.pal(3, "Pastel2")
  
SET1 <- sig_list1$geneid
SET2 <- sig_list2$geneid
SET3 <- sig_list3$geneid

  
  
venn.diagram(
    x = list(SET1, SET2, SET3),
    category.names = c("vgPcKD" , "sdPcKD" , "PcKD"),
    fill=myCol,
    filename = 'PC_sdPc_vgPC_whole_1.png',
    output=TRUE
  )
  
SET4 = sig_list1_enriched$geneid
SET5 = sig_list2_enriched$geneid
SET6 = sig_list3_enriched$geneid

venn.diagram(
  x = list(SET4, SET5, SET6),
  category.names = c("vgPcKD" , "sdPcKD" , "PcKD"),
  fill=myCol,
  filename = 'PC_sdPc_vgPC_enriched.png',
  output=TRUE
)
     
SET7 = sig_list1_depleted$geneid
SET8 = sig_list2_depleted$geneid
SET9 = sig_list3_depleted$geneid

venn.diagram(
  x = list(SET7, SET8, SET9),
  category.names = c("vgPcKD" , "sdPcKD" , "PcKD"),
  fill=myCol,
  filename = 'PC_sdPc_vgPC_depleted.png',
  output=TRUE
)

##########
#pull the genes from each section of the venn diagram
overlap = calculate.overlap(list(SET1, SET2, SET3))

#Pull appropriate rows from list to create a data frame and write that df to a tableprint(overlap[[1]])
overlap_df = data.frame(overlap[[1]]) #a5 (vgPc+sdPc+Pc = 243 genes)
write.table(overlap_df, "vgPc_sdPc_Pc_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlap[[2]])
overlap_df_2 = data.frame(overlap[[2]]) #a2 (vgPc and sdPc = 723 genes) 
write.table(overlap_df_2, "vgPc_sdPc_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlap[[3]])
overlap_df_3 = data.frame(overlap[[3]]) #a4 (vgPc and Pc = 78 genes) 
write.table(overlap_df_3, "vgPc_Pc_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlap[[4]])
overlap_df_4 = data.frame(overlap[[4]]) #a6 (sdPc and Pc = 88 genes) 
write.table(overlap_df_4, "sdPc_Pc_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlap[[5]])
overlap_df_5 = data.frame(overlap[[5]]) #a1 (vgPc = 412 genes) 
write.table(overlap_df_5, "vgPc_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlap[[6]])
overlap_df_6 = data.frame(overlap[[6]]) #a3 (sdPc = 415 genes) 
write.table(overlap_df_6, "sdPc_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlap[[7]])
overlap_df_7 = data.frame(overlap[[7]]) #a7 (Pc = 299 genes) 
write.table(overlap_df_7, "Pc_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)


##########
#pull the genes from each section of the enriched venn diagram
overlapE_list = calculate.overlap(list(SET4, SET5, SET6))
View(overlapE_list)

#Pull appropriate rows from list to create a data frame and write that df to a table
print(overlapE_list[[1]])
overlapE_df = data.frame(overlapE_list[[1]]) #a5 (vgPc+sdPc+Pc = 61 genes)
write.table(overlapE_df, "vgPc_sdPc_Pc_overlapE_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapE_list[[2]])
overlapE_df_2 = data.frame(overlapE_list[[2]]) #a2 (vgPc and sdPc = 107 genes) 
write.table(overlapE_df_2, "vgPc_sdPc_overlapE_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapE_list[[3]])
overlapE_df_3 = data.frame(overlapE_list[[3]]) #a4 (vgPc and Pc = 22 genes) 
write.table(overlapE_df_3, "vgPc_Pc_overlapE_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapE_list[[4]])
overlapE_df_4 = data.frame(overlapE_list[[4]]) #a6 (sdPc and Pc = 34 genes) 
write.table(overlapE_df_4, "sdPc_Pc_overlapE_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapE_list[[5]])
overlapE_df_5 = data.frame(overlapE_list[[5]]) #a1 (vgPc = 55 genes) 
write.table(overlapE_df_5, "vgPc_overlapE_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapE_list[[6]])
overlapE_df_6 = data.frame(overlapE_list[[6]]) #a3 (sdPc = 174 genes) 
write.table(overlapE_df_6, "sdPc_overlapE_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapE_list[[7]])
overlapE_df_7 = data.frame(overlapE_list[[7]]) #a7 (Pc = 278 genes) 
write.table(overlapE_df_7, "Pc_overlapE_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

##########
#pull the genes from each section of the depleted venn diagram
overlapD_list = calculate.overlap(list(SET7, SET8, SET9))
View(overlapD_list)

#Pull appropriate rows from list to create a data frame and write that df to a table
print(overlapD_list[[1]])
overlapD_df = data.frame(overlapD_list[[1]]) #a5 (vgPc+sdPc+Pc = 111 genes)
write.table(overlapD_df, "vgPc_sdPc_Pc_overlapD_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapD_list[[2]])
overlapD_df_2 = data.frame(overlapD_list[[2]]) #a2 (vgPc and sdPc = 676 genes) 
write.table(overlapD_df_2, "vgPc_sdPc_overlapD_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapD_list[[3]])
overlapD_df_3 = data.frame(overlapD_list[[3]]) #a4 (vgPc and Pc = 30 genes) 
write.table(overlapD_df_3, "vgPc_Pc_overlapD_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapD_list[[4]])
overlapD_df_4 = data.frame(overlapD_list[[4]]) #a6 (sdPc and Pc = 32 genes) 
write.table(overlapD_df_4, "sdPc_Pc_overlapD_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapD_list[[5]])
overlapD_df_5 = data.frame(overlapD_list[[5]]) #a1 (vgPc = 394 genes) 
write.table(overlapD_df_5, "vgPc_overlapD_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapD_list[[6]])
overlapD_df_6 = data.frame(overlapD_list[[6]]) #a3 (sdPc = 274 genes) 
write.table(overlapD_df_6, "sdPc_overlapD_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

print(overlapD_list[[7]])
overlapD_df_7 = data.frame(overlapD_list[[7]]) #a7 (Pc = 140 genes) 
write.table(overlapD_df_7, "Pc_overlapD_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

##########
#GO analysis of genes in each list
#enriched comparison
library(org.Dm.eg.db)
library(topGO)
library(GOplot)

E_vgPc_sdPc_Pc = read.csv("~/Desktop/testing/hyperproliferation/enriched_venn/vgPc_sdPc_Pc_overlapE_list.txt")
View(E_vgPc_sdPc_Pc)

E_geneList = E_vgPc_sdPc_Pc[,1]
head(E_geneList)

E_ego = enrichGO(E_geneList, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
head(E_ego)

dotplot(E_ego, showCategory=20)