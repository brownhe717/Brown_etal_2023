#library(tsrexplorer)
library(GenomicFeatures)
library(Gviz)
library(viridis)
library(scales)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# This script imports signal tracks and displays them at the specified regions of the Drosophila 

gviz_dir = file.path("~/Desktop/Brown et al 2022/code/Figure 5")
bigwig_dir = file.path("~/Desktop/Brown et al 2022/code/Figure 5/INPUT/bigs")
rnaseq_dir <- file.path("~/Desktop/Brown et al 2022/code/Figure 5/INPUT/RNA_peaks")

# Create genomic axis track
axis.track <- GenomeAxisTrack(col = "black", scale = 0.1, col.range = "black")

options(ucscChromosomeNames = FALSE)

# Create gene annotation track
txdb <- makeTxDbFromGFF("~/Desktop/Brown et al 2022/code/Figure 3/INPUT/Dm_genome/genes.gtf")

genome.track <- GeneRegionTrack(txdb, genome = "dm6", shape = "arrow", names = "Genes", col = "black",
                                showId = FALSE, fill = "black", trancriptAnnotation = "gene_symbol", collapseTranscripts = "meta")


# Create data tracks

# Get colors from the viridis palette for the number of tracks to be plotted
# Show_col(viridis_pal()(40))



#########################
### Set y-axis EAD limits ###
#########################


########## EAD ##########
# H2AK119Ub track (replicate 2)
H2AK119Ub_EAD_WT_pos_lim <- c(0,3000)
# H2AK119Ub_WT_neg_lim <- c(0,1500)

# H3K27me3 track (58-59 10 EADs)
H3K27me3_EAD_WT_pos_lim = c(0,60000)
# H3K27me3_WT_3instar_neg_lim = c(0, )

# RNAseq track (EAD)
transcript_EAD_pos_lim = c(0,20)
transcript_EAD_neg_lim = c(0,20)

# IgG (H2AK119Ub IgG replicate 3)
IgG_EAD_K27_pos_lim = c(0,60000)
IgG_EAD_Ub_pos_lim = c(0,3000)


H2AK119Ub_EAD <- DataTrack(range = file.path(bigwig_dir, "WT_EAD_H2AK119Ub_1_SF.bigwig"), genome = "dm6", 
                           name = "H2AK119Ub", col.histogram = "#3A538BFF", fill.histogram = "#3A538BFF", ylim = H2AK119Ub_EAD_WT_pos_lim)

# H3K27me3 EAD
H3K27me3_EAD <- DataTrack(range = file.path(bigwig_dir, "WT_EAD_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                          name = "H3K27me3", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = H3K27me3_EAD_WT_pos_lim)
#try changing 120hr K27me3 to #404588FF

#RNA-seq 120hr WT ED
WT_EAD_120_1_pos <- DataTrack(range = file.path(rnaseq_dir, "HB01_EAD_120hrs_1_test.CPM.bs1.smooth25.plus.bw"), genome = "dm6", 
                              name = "EAD transcript", col.histogram = "#4CC26Cff", fill.histogram = "#4CC26Cff", ylim = transcript_EAD_pos_lim)
WT_EAD_120_1_neg <- DataTrack(range = file.path(rnaseq_dir, "HB01_EAD_120hrs_1_test.CPM.bs1.smooth25.minus.bw"), genome = "dm6", 
                              name = "EAD transcript", col.histogram = "#BADE28FF", fill.histogram = "#BADE28FF", ylim = transcript_EAD_neg_lim)

# IgG EAD
#IgG_EAD <- DataTrack(range = file.path(bigwig_dir4, "10WT_3rd_instar_IgG_Ub_2_normalized.bigwig"), genome = "dm6", 
#name = "IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_EAD_pos_lim)
IgG_EAD_K27 <- DataTrack(range = file.path(bigwig_dir, "WT_EAD_IgG_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                         name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_EAD_K27_pos_lim)
IgG_EAD_Ub <- DataTrack(range = file.path(bigwig_dir, "WT_EAD_IgG_H2AK119Ub_1_SF.bigwig"), genome = "dm6", 
                        name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_EAD_Ub_pos_lim)


#####################
### Overlay Track ###
#####################

displayPars(IgG_EAD_K27) = list(alpha.title=1,alpha=1)
displayPars(H3K27me3_EAD) = list(alpha.title=1, alpha=1)
ot_H3K27me3_EAD = OverlayTrack(trackList = list(H3K27me3_EAD, IgG_EAD_K27))

displayPars(IgG_EAD_Ub) = list(alpha.title=1,alpha=1)
displayPars(H2AK119Ub_EAD) = list(alpha.title=1, alpha=1)
ot_H2AK119Ub_EAD = OverlayTrack(trackList = list(H2AK119Ub_EAD, IgG_EAD_Ub))

displayPars(WT_EAD_120_1_pos) = list(alpha.title=1, alpha=0.8)
displayPars(WT_EAD_120_1_neg) = list(alpha.title=1, alpha=1)
otRNA_EAD = OverlayTrack(trackList = list(WT_EAD_120_1_neg, WT_EAD_120_1_pos))


##############################
### Set y-axis limits PcKD ###
##############################

#H2AK119Ub track (replicate 2)
H2AK119Ub_PcKD_pos_lim <- c(0,3000)
#H2AK119Ub_WT_neg_lim <- c(0,1500)

#H3K27me3 track (58-59 10 EADs)
H3K27me3_PcKD_pos_lim = c(0,60000)
#H3K27me3_WT_3instar_neg_lim = c(0, )

#RNAseq track (PcKD)
transcript_PcKD_pos_lim = c(0,20)
transcript_PcKD_neg_lim = c(0,20)

#IgG (H2AK119Ub IgG replicate 3)
IgG_PcKD_pos_lim = c(0,2000)
IgG_PcKD_K27_pos_lim = c(0,60000)
IgG_PcKD_Ub_pos_lim = c(0,3000)

# H2AK119Ub PcKD
H2AK119Ub_PcKD <- DataTrack(range = file.path(bigwig_dir, "PcKD_H2AK119Ub_2_SF.bigwig"), genome = "dm6", 
                            name = "H2AK119Ub", col.histogram = "#3A538BFF", fill.histogram = "#3A538BFF", ylim = H2AK119Ub_PcKD_pos_lim)

# H3K27me3 PcKD
H3K27me3_PcKD <- DataTrack(range = file.path(bigwig_dir, "PcKD_H3K27me_3_SF.bigwig"), genome = "dm6", 
                           name = "H3K27me3", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = H3K27me3_PcKD_pos_lim)
#try changing 120hr K27me3 to #404588FF

#RNA-seq 144hr PcKD EAD
PcKD_EAD_144_1_pos <- DataTrack(range = file.path(rnaseq_dir, "HB13_EAD_144hr_1_test.CPM.bs1.smooth25.plus.bw"), genome = "dm6", 
                                name = "PcKD transcript", col.histogram = "#4CC26Cff", fill.histogram = "#4CC26Cff", ylim = transcript_PcKD_pos_lim)
PcKD_EAD_144_1_neg <- DataTrack(range = file.path(rnaseq_dir, "HB13_EAD_144hr_1_test.CPM.bs1.smooth25.minus.bw"), genome = "dm6", 
                                name = "PcKD transcript", col.histogram = "#BADE28FF", fill.histogram = "#BADE28FF", ylim = transcript_PcKD_neg_lim)

# IgG PcKD
#IgG_PcKD <- DataTrack(range = file.path(bigwig_dir4, "PcKD_H3K27me_IgG_2_SF.bigwig"), genome = "dm6", 
#name = "IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_PcKD_pos_lim)
IgG_PcKD_K27 <- DataTrack(range = file.path(bigwig_dir, "PcKD_H3K27me_IgG_2_SF.bigwig"), genome = "dm6", 
                          name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_PcKD_K27_pos_lim)
IgG_PcKD_Ub <- DataTrack(range = file.path(bigwig_dir, "PcKD_H2AK119Ub_IgG_2_SF.bigwig"), genome = "dm6", 
                         name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_PcKD_Ub_pos_lim)

#####################
### Overlay Track ###
#####################

displayPars(IgG_PcKD_K27) = list(alpha.title=1,alpha=1)
displayPars(H3K27me3_PcKD) = list(alpha.title=1, alpha=1)
ot_H3K27me3_PcKD = OverlayTrack(trackList = list(H3K27me3_PcKD, IgG_PcKD_K27))

displayPars(IgG_PcKD_Ub) = list(alpha.title=1,alpha=1)
displayPars(H2AK119Ub_PcKD) = list(alpha.title=1, alpha=1)
ot_H2AK119Ub_PcKD = OverlayTrack(trackList = list(H2AK119Ub_PcKD, IgG_PcKD_Ub))

displayPars(PcKD_EAD_144_1_pos) = list(alpha.title=1, alpha=0.8)
displayPars(PcKD_EAD_144_1_neg) = list(alpha.title=1, alpha=0.8)
otRNA_PcKD = OverlayTrack(trackList = list(PcKD_EAD_144_1_neg, PcKD_EAD_144_1_pos))


############################
### Set y-axis limits WD ###
############################

#H2AK119Ub track (replicate 2)
H2AK119Ub_WT_WD_pos_lim <- c(0,3000)
#H2AK119Ub_WT_neg_lim <- c(0,1500)

#H3K27me3 track (58-59 10 EADs)
H3K27me3_WT_WD_pos_lim = c(0,60000)
#H3K27me3_WT_3instar_neg_lim = c(0, )

#RNAseq track (WD)
transcript_WD_pos_lim = c(0,15)
transcript_WD_neg_lim = c(0,15)

#IgG (H2AK119Ub IgG replicate 3)
#IgG_WD_pos_lim = c(0,2000)
IgG_WD_K27_pos_lim = c(0,6000)
IgG_WD_Ub_pos_lim = c(0,3000)


# H2AK119Ub WD
H2AK119Ub_WD <- DataTrack(range = file.path(bigwig_dir, "WT_WD_H2AK119Ub_3_SF.bigwig"), genome = "dm6", 
                          name = "H2AK119Ub", col.histogram = "#3A538BFF", fill.histogram = "#3A538BFF", ylim = H2AK119Ub_WT_WD_pos_lim)

# H3K27me3 WD
H3K27me3_WD <- DataTrack(range = file.path(bigwig_dir, "WT_WD_H3K27me_1_SF.bigwig"), genome = "dm6", 
                         name = "H3K27me3", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = H3K27me3_WT_WD_pos_lim)
#try changing 120hr K27me3 to #404588FF

#RNA-seq 120hr WT Wd
WT_WD_120_1_pos <- DataTrack(range = file.path(rnaseq_dir, "HB04_WD_120hrs_1_test.CPM.bs1.smooth25.plus.bw"), genome = "dm6", 
                             name = "WD transcript", col.histogram = "#4CC26Cff", fill.histogram = "#4CC26Cff", ylim = transcript_WD_pos_lim)
WT_WD_120_1_neg <- DataTrack(range = file.path(rnaseq_dir, "HB04_WD_120hrs_1_test.CPM.bs1.smooth25.minus.bw"), genome = "dm6", 
                             name = "WD transcript", col.histogram = "#BADE28FF", fill.histogram = "#BADE28FF", ylim = transcript_WD_neg_lim)

# IgG WD
#IgG_WD <- DataTrack(range = file.path(bigwig_dir4, "WT_WD_IgG_H2AK119Ub_3_normalized.bigwig"), genome = "dm6", 
#name = "IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_WD_pos_lim)
IgG_WD_K27 <- DataTrack(range = file.path(bigwig_dir, "WT_WD_H3K27me_IgG_1_SF.bigwig"), genome = "dm6", 
                        name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_WD_K27_pos_lim)
IgG_WD_Ub <- DataTrack(range = file.path(bigwig_dir, "WT_WD_IgG_H2AK119Ub_3_SF.bigwig"), genome = "dm6", 
                       name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_WD_Ub_pos_lim)

#####################
### Overlay Track ###
#####################

displayPars(IgG_WD_K27) = list(alpha.title=1,alpha=1)
displayPars(H3K27me3_WD) = list(alpha.title=1, alpha=1)
ot_H3K27me3_WD = OverlayTrack(trackList = list(H3K27me3_WD, IgG_WD_K27))

displayPars(IgG_WD_Ub) = list(alpha.title=1,alpha=1)
displayPars(H2AK119Ub_WD) = list(alpha.title=1, alpha=1)
ot_H2AK119Ub_WD = OverlayTrack(trackList = list(H2AK119Ub_WD, IgG_WD_Ub))

displayPars(WT_WD_120_1_pos) = list(alpha.title=1, alpha=0.8)
displayPars(WT_WD_120_1_neg) = list(alpha.title=1, alpha=0.8)
otRNA_WD = OverlayTrack(trackList = list(WT_WD_120_1_neg, WT_WD_120_1_pos))

###################
### Plot Tracks ###
###################

# EAD vg locus
plotTracks(list(axis.track, 
                ot_H3K27me3_EAD, #60,000
                ot_H2AK119Ub_EAD,#1500
                otRNA_EAD,
                genome.track), 
           chromosome = "chr2R",from = 12884201, to = 12900600, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)
#dev.off()

# PcKD vg locus
plotTracks(list(axis.track, 
                ot_H3K27me3_PcKD,#y-limit 15,000
                ot_H2AK119Ub_PcKD,
                otRNA_PcKD,
                genome.track), 
           chromosome = "chr2R",from = 12884201, to = 12900600, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)

# WD vg locus
plotTracks(list(axis.track, 
                ot_H3K27me3_WD,#y-limit 15,000
                ot_H2AK119Ub_WD,
                otRNA_WD,
                genome.track), 
           chromosome = "chr2R",from = 12884201, to = 12900600, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)

