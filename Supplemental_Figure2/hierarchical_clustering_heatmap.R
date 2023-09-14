library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)

# Read the data into R -- contains information about genes (geneid, length, # of reads aligned to the gene in each sample)
## 3 replicates for each cell type and time point
seqdata <- read.delim("~/Desktop/Brown et al 2023/scripts/SFigure 1/counts_merged.txt", stringsAsFactors = FALSE)

# Read the sample information into R
sampleinfo <- read.delim("~/Desktop/Brown et al 2023/scripts/SFigure 1/SampleInfo.txt", stringsAsFactors = TRUE)

#look ar the first 6 lines of data (head) and how many rows and columsn the data frame has (dim)
head(seqdata)
dim(seqdata)

#view sampleinfo variable in terminal -- file name, sample name, cell type, and status
sampleinfo

##For analysis, make a new matrix from the seqdata data frame containing only the counts & store gene identifiers as rownames

# Remove first two columns from seqdata -- create a new data object that contains only the counts for the 12 samples
countdata <- seqdata[,-(1:1)]
# Look at the output
head(countdata)

#Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)

#Take the column names
colnames(countdata)

#Shorten sample names to contain only the relevent info (substr)
##using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
#colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
#head(countdata)

#Note: colnames are now the same as SampleName in sampleinfo variable
table(colnames(countdata)==sampleinfo$SampleName)

#Filter to remove lowly expressed genes - filter on a minimum counts per million threshold present in 'n' (the smallest sample size for each group in the experiment)
##retain genes if they are expressed at a counts-per-milliom (CPM) above 0.5 in at least 'n' samples
###obtain cpm's with the cpm function in the edgeR library -- converting to CPMs normalizes for the different sequencing depths for each sample
myCPM = cpm(countdata)
head(myCPM)

#set a threshold for values greater than 0.5 -- corresponds to a count of 10-15 based on this library size (smaller CPM thresholds are used for larger libraries)
thresh = myCPM > 1.0
head(thresh)

#summarize how many 'trues' there are in each tow
table(rowSums(thresh))

#keep genes that have at least 2 trues in each row of thresh
keep = rowSums(thresh) >= 3

#subset the rowns of countdata to keep the more highly expressed genes
counts.keep = countdata[keep,]
summary(keep)
dim(counts.keep)

#check to see if the CPM threshold of 0.5 corresponds to ~10-15 counts
##limit the x and y axis to view smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
#add a vertical line at 0.5 CPM
abline(v=1.0)

#Convert counts to DGEList object -- object used by edgeR to store count data
y = DGEList(counts.keep)
#view y object in terminal
y

#see what slots are stored in y
names(y)
#Library size information is stored in the samples slot
y$samples

#Quality Control Plots
##library sizes (barplot) and distribution plots -- see whether there are any major discrepancies between the samples
y$samples$lib.size

#note: names argument in barplot function tells the barplot to use the sample names on the x-axis, las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Barplot of the library sizes")

#examine the distribution of raw counts by getting the log2 counts per million (since count data is not normally distributed)
logcounts = cpm(y,log=TRUE)
#check distributions of samples with boxplots -- add a blue horizontal line at the median logCPM (abline)
boxplot(logcounts,xlab="",ylab="Log2 counts per million",las=2)
abline(h=median(logcounts), col="blue")
title("Boxplots of logCPMs (unnormalized)")

#note: if a sample is far above or below the median line its quality may need to be further investigated

###Multidimensional Scaling Plots (MDSplots) -- visualize principle component analysis to determine sources of variation in the data
    #also useful to quality control and checking for outliers
#color the samples according to group information (CellType -- "basal" and "luminal"; Status -- "lactate", "pregnant" and "virgin")
#plot two plots side by side -- Cell Type first, then Status
###par(mfrow=c(1,2))
levels(sampleinfo$Tissue)
  
#code for Tissue plot
col.tissue = c("purple","orange")[sampleinfo$Tissue]
data.frame(sampleinfo$Tissue,col.tissue)
plotMDS(y,col=col.tissue)

#add a legend to the plot to tell which colors correspond to which cell types
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$Tissue))
title("Tissue")

  #code for Genotype plot
col.genotype <- c("purple", "blue", "orange")[sampleinfo$Genotype]
data.frame(sampleinfo$Genotype,col.genotype)

plotMDS(y,col=col.genotype)
legend("topleft",fill=c("blue","orange", "purple"),legend=levels(sampleinfo$Genotype),cex=0.8)
title("Genotype")

#generate an interactive MDS plot using the Glimma package (opens Chrome - can click on data points for detail)
##Glimma can be used to obtain MDS plots and mean-difference (MD) plots
###output will show the MDS plot and the amount of variation explained by each condition/'dimension' in a barplot
labels <- paste(sampleinfo$SampleName, sampleinfo$Tissue)
group <- paste(sampleinfo$Tissue,sep=".")
group <- factor(group)
glMDSPlot(y, labels=labels, groups=group, folder="mds")

#Hierarchical Clustering with Heatmaps
##heatmap.2 function from gplots package -- calculates a matrix of euclidean distances from the logcounts object (logCPM) for the 500 most variable genes
##RColorBrewer package has color schemes accesed with the brewer.pal function ("RdYiBu" or "Spectral" are common choices)
##png function will create a png file to save the plots to -- close the file with dev.off()

#estimate the variance for each row in logcounts matrix ###get different numerical output from this point forward###
var_genes = apply(logcounts, 1, var)
head(var_genes)

#get the gene names for the top 500 most variable genes
select_var = names(sort(var_genes, decreasing=TRUE))[1:300]
head(select_var)

#subset logcounts matrix  ###same output for dim, different for head####
highly_variable_lcpm = logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

#get a nicer color palette aesthetic
mypalette = brewer.pal(11,"RdYlBu")
morecols = colorRampPalette(mypalette)

#set up a color vector for celltype variable
col.genotype = c("purple", "orange", "blue")[sampleinfo$Genotype]

#plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 300 most variable genes across samples",ColSideColors=col.genotype,scale="row", margins = c(9,5))

col.cell = c("purple", "orange", "yellow", "dark green")[sampleinfo$Hr_AEL]

#plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
