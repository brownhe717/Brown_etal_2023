#!/usr/env/bin Rscript
library(dplyr)
library(tidyverse)
library(VennDiagram)
library(RColorBrewer)


###########################
###### Compile Lists ######
###########################


# Read in tsv files of ALL target genes & store each as a variable (time_list), convert to a data frame (time_list_df), 
# Then store each list of genes as its own variable (sig_time_list)


########## All Candidates ##########

early_list = read.delim("target_genes_early.tsv",sep="\t",row.names="geneid")
early_list_df = as.data.frame(early_list)

sig_early_list = rownames_to_column(early_list_df, "geneid") %>% dplyr::select("geneid")
count(sig_early_list)

##########    

mid_list = read.delim("target_genes_mid.tsv",sep="\t",row.names="geneid")
mid_list_df = as.data.frame(mid_list)

sig_mid_list = rownames_to_column(mid_list_df, "geneid") %>%
  dplyr::select("geneid")
count(sig_mid_list)

##########

late_list = read.delim("target_genes_late.tsv",sep="\t",row.names="geneid")
late_list_df = as.data.frame(late_list)

sig_late_list = rownames_to_column(late_list_df, "geneid") %>%
  dplyr::select("geneid")
count(sig_late_list)


##########################
###### Venn Diagram ######
##########################


# Venn Diagram for all targets
myCol <- brewer.pal(3, "Pastel2")

SET1 <- sig_early_list$geneid
SET2 <- sig_mid_list$geneid
SET3 <- sig_late_list$geneid

venn.diagram(
  x = list(SET1, SET2, SET3),
  category.names = c("early" , "mid" , "late"),
  fill=myCol,
  filename = 'total_targets.tiff',
  output=TRUE
)



# Venn Diagram for known targets
myCol <- brewer.pal(3, "Pastel2")

SET4 <- sig_early_list_K$geneid
SET5 <- sig_mid_list_K$geneid
SET6 <- sig_late_list_K$geneid

venn.diagram(
  x = list(SET4, SET5, SET6),
  category.names = c("early" , "mid" , "late"),
  fill=myCol,
  filename = 'known_targets.tiff',
  output=TRUE
)



# Venn Diagram for unknown targets
myCol <- brewer.pal(3, "Pastel2")

SET7 <- sig_early_list_U$geneid
SET8 <- sig_mid_list_U$geneid
SET9 <- sig_late_list_U$geneid

venn.diagram(
  x = list(SET7, SET8, SET9),
  category.names = c("early" , "mid" , "late"),
  fill=myCol,
  filename = 'unknown_targets.tiff',
  output=TRUE
)


#########################
###### Categories ######
#########################

########## All Candidates ##########
# Pull the genes from each section of the venn diagram
overlap = calculate.overlap(list(SET1, SET2, SET3))
View(overlap)

#Pull appropriate rows from list to create a data frame and write that df to a tableprint(overlap[[1]])
overlap_df = data.frame(overlap[[1]]) #a5 (early + mid + late = 13 genes)
write.table(overlap_df, "early_mid_late_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

overlap_df_2 = data.frame(overlap[[2]]) #a2 (early + mid = 22 genes) 
write.table(overlap_df_2, "early_mid_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

overlap_df_3 = data.frame(overlap[[3]]) #a4 (early + late = 37 genes) 
write.table(overlap_df_3, "early_late_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

overlap_df_4 = data.frame(overlap[[4]]) #a6 (mid + late = 23 genes) 
write.table(overlap_df_4, "mid_late_overlap_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

overlap_df_5 = data.frame(overlap[[5]]) #a1 (early = 444 genes) 
write.table(overlap_df_5, "early_targts_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

overlap_df_6 = data.frame(overlap[[6]]) #a3 (mid = 162 genes) 
write.table(overlap_df_6, "mid_targets_list.txt", sep="\t", row.names=T, col.names=T, quote=F)

overlap_df_7 = data.frame(overlap[[7]]) #a7 (late = 160 genes) 
write.table(overlap_df_7, "late_targets_list.txt", sep="\t", row.names=T, col.names=T, quote=F)