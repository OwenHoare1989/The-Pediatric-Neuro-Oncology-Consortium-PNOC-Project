library(MASS) 
library(reshape2) 
library("FactoMineR")
library("factoextra")
library(heatmaply)
library(ggpubr)
library(DESeq2)
library(glmGamPoi)
library(EnhancedVolcano)
library(BiocParallel)
register(MulticoreParam(12))
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)
library("tibble")
library(reshape2)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
library(clusterProfiler)
library(enrichplot)
require(DOSE)
library(europepmc)
library(pathview)
library(survival)
library(survminer)
require("survival")
library("forestmodel")
library(devtools)
#install_github("KatrionaGoldmann/volcano3D", force = TRUE)
library(volcano3D)
#install_github("KatrionaGoldmann/volcano3Ddata")
library(volcano3Ddata)
library(volcano3D)
library(kableExtra)
library("ggrepel") 
library("plotly") 
library("dplyr")
library(knitr)
library(htmltools)
library(htmlTable)
library(htmlwidgets)
library( digest)
library( DT)
library( fs)
library(future)
library(ggrepel)
library(ggtree)
library(parallel)
library(tibble)
library(data.table)
library(speckle)
library(rstatix)
plan("multicore", workers = 12)
options(future.globals.maxSize = (1014*1024^2)*60)

# Parameters for plots
theme_set(
  theme_bw() +
    theme(
      plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size = 14)
    )
)


## @knitr PrepareData
setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/09_Bulk_RNA/Sabien/Comparison_across_strata")

# Read in the raw read counts
rawCounts <- read.csv("Combined.all_strata.exprs_matrix.csv", row.names = 1)
# NormDF <- QuantNorm(rawCounts)
# NormDF <- as.data.frame(NormDF)
# NormDF <- rownames_to_column(NormDF, var = "Gene")

#head(rawCounts)

# Read in the sample mappings
sampleData <- read.csv("Combined.all_strata.meta_data.csv")
#head(sampleData)

# Also save a copy for later
sampleData_v2 <- sampleData

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Gene
sampleIndex <- grepl("PNOC005.\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
#head(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
#head(sampleData)
rownames(sampleData) <- sampleData$Mixture
keep <- c("Strata", "Days", "Strata_timepoint")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("Strata", "Time", "Strata_Timepoint")
sampleData$Strata_Timepoint <- factor(sampleData$Strata_Timepoint)
#head(sampleData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# rename the tissue types
Strata_Timepoint <- function(x){
x <- switch(as.character(x), "A_D0" = "A_D0", "A_D14" = "A_D14", "A_D28" = "A_D28", "B_D1" = "B_D1",  "B_D4" = "B_D4",  "B_D7" = "B_D7",
                              "B_D14" ="B_D14", "B_D28" = "B_D28", "C_D28" = "C_D28", "C_D7" = "C_D7",  "C_D0" = "C_D0",  "C_D56" = "C_D56", "B_D0" = "B_D0")
    return(x)
}
sampleData$Strata_Timepoint <- unlist(lapply(sampleData$Strata_Timepoint, Strata_Timepoint))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$Strata <- factor(sampleData$Strata_Timepoint, levels=c("A_D0", "A_D14", "A_D28", "B_D1",  "B_D4",  "B_D7",  "B_D14", "B_D28", "C_D28", "C_D7",  "C_D0",  "C_D56", "B_D0"))

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ Strata_Timepoint)

dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])

# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]
# Register the number of cores to use

## @knitr Run_DESeq2
# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data, test="LRT", reduced=~1, parallel = TRUE)

# Extract differential expression results
B_C_D28 <- results(deseq2Data, contrast=c("Strata_Timepoint",  "B_D28", "C_D28"))
# Coerce to a data frame
B_C_D28 <- as.data.frame(B_C_D28)
B_C_D28 <- rownames_to_column(B_C_D28, var = "gene")
B_C_D28$significant <- ifelse(B_C_D28$padj < .05, "Significant", NA)
B_C_D28<-na.omit(B_C_D28)

# write these differentially expressed genes as a table and save to directory
write.csv(B_C_D28, file= "DEG_strata_B_vs_C_D28.csv")

## @knitr Print_Tables
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(B_C_D28, server = TRUE,
                                              class = "compact",
                                              filter="top",
                                              rownames = FALSE,
                                              colnames = c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant"),
                                              extensions = c('Buttons'),
                                              options = list(
                                                pageLength = 15,
                                                dom = 'Bfrtip',
                                                buttons = c('excel','csv','pdf','copy')
                                              ))))


## @knitr Contrast1
# Extract differential expression results
D0_v_D14_SA <- results(deseq2Data, contrast=c("Strata",  "D0", "D14"))
# Coerce to a data frame
B_C <- as.data.frame(B_C)
B_C <- rownames_to_column(B_C, var = "gene")
B_C$significant <- ifelse(B_C$padj < .05, "Significant", NA)
B_C<-na.omit(B_C)

## @knitr Print_Tables1
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(B_C, server = TRUE,
                                              class = "compact",
                                              filter="top",
                                              rownames = FALSE,
                                              colnames = c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant"),
                                              extensions = c('Buttons'),
                                              options = list(
                                                pageLength = 15,
                                                dom = 'Bfrtip',
                                                buttons = c('excel','csv','pdf','copy')
                                              ))))


## @knitr Contrast2
# Extract differential expression results
A_C <- results(deseq2Data, contrast=c("Strata",  "Strata A", "Strata C"))
# Coerce to a data frame
A_C <- as.data.frame(A_C)
A_C <- rownames_to_column(A_C, var = "gene")
A_C$significant <- ifelse(A_C$padj < .05, "Significant", NA)
A_C<-na.omit(A_C)
write.csv(A_C, file="DEG_strata_A_vs_C.csv")

## @knitr Print_Tables2
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(A_C, server = TRUE,
                                              class = "compact",
                                              filter="top",
                                              rownames = FALSE,
                                              colnames = c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant"),
                                              extensions = c('Buttons'),
                                              options = list(
                                                pageLength = 15,
                                                dom = 'Bfrtip',
                                                buttons = c('excel','csv','pdf','copy')
                                              ))))



# Extract differential expression results
## @knitr DEGSeq2_Analysis
deseq2Results <- results(deseq2Data, contrast=c("Time",  "D0", "D7"))


## @knitr DispersionPlot
plotDispEsts(deseq2Data, ylim = c(1e-6, 1e1) )

# View summary of results
#summary(deseq2Results)


# Using DEseq2 built in method
# plotMA(deseq2Results)
# Load libraries
# install.packages(c("ggplot2", "scales", "viridis"))


# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# Examine this data frame
#head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .05, "Significant", NA)

# Plot the results similar to DEseq2
MAPlot1 <- ggplot(deseq2ResDF,  aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw() + ggtitle("Control vs MV Treatment")

# Let's add some more detail
MAPlot2 <- ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2) + ggtitle("Control vs MV Treatment")

## @knitr MAPlots
MAPlot1 + MAPlot2

# Add rectangle around labels
deseq2ResDF2 <- rownames_to_column(deseq2ResDF, var = "gene")
# Add rectangle around labels

## @knitr MAPlots_Significant_Genes
ggmaplot(deseq2ResDF2, main = expression("Day 0" %->% "Day 7"),
         fdr = 0.05, fc = 1, size = 1.2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(deseq2ResDF2$gene),
         legend = "top", top = 50,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())


## @knitr VolcanoPlot
EnhancedVolcano(deseq2ResDF2,
                lab = as.character(deseq2ResDF2$gene),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                title = "Volcano Plot DEG (Day 0 vs Day 7)",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 2,
                labSize = 4,
                colAlpha = 1,
                
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 2,
                drawConnectors = TRUE,
                widthConnectors = 0.75)



# Extract counts for the gene otop2
otop2Counts <- plotCounts(deseq2Data, gene="IFI27", intgroup=c("Time", "Sample"), returnData=TRUE)

# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
## @knitr PlotCounts
ggplot(otop2Counts, aes(x=Sample, y=count, colour=Time, group=Time)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("IFI27 (Interferon Alpha Inducible Protein 27)")

deseq2ResDF["IFI27",]
rawCounts["IFI27",]
normals=row.names(sampleData[sampleData[,"Time"]=="D0",])
primaries=row.names(sampleData[sampleData[,"Time"]=="D7",])
rawCounts["IFI27",normals]
rawCounts["IFI27",primaries]

## @knitr HeatmapTopGenes
# Transform count data using the variance stablilizing transform
deseq2VST <- vst(deseq2Data)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
#head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 2,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2


# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

#head(deseq2VST_wide)
#head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
# heatmap

# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="ward.D2")
clusterSample <- hclust(distanceSample, method="ward.D2")

# Construct a dendogram for samples
# install.packages("ggdendro")

sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
#heatmap

# Combine the dendrogram and the heatmap
# install.packages("gridExtra")

#grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))


# Load in libraries necessary for modifying plots
#install.packages("gtable")


# Modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Check the widths of each grob
# sampleDendrogramGrob$widths
# heatmapGrob$widths

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
#grid.draw(finalGrob)


# Re-order the sample data to match the clustering we did
sampleData_v2$Mixture <- factor(sampleData_v2$Mixture, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52", "#00AFBB", "#E7B800", "#FC4E07", "green", "red", "navy")
sampleClinical <- ggplot(sampleData_v2, aes(x=Mixture, y=1, fill=Days)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Time", values=colours) + theme_void()

# Convert the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

# Make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# Arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
#grid.draw(finalGrob)


###############################################################################
################# Step 1: create dendrogram for genes ##########################

# we already did the clustering for genes in the tutorial, get the data to make a dendrogram with ggplot
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))

# construct the dendrogram in ggplot
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) + theme_dendro()

################################################################################
################# Step 2: Re-arrange the heatmap cells #########################

# re-factor genes for ggplot2
deseq2VST$Gene <- factor(deseq2VST$Gene, levels=clusterGene$labels[clusterGene$order])

# recreate the heatmap with this new factoring
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())

################################################################################
################# Step 3: convert to everything to grobs #######################

# note! before this step as mentioned you might need to alter the expand parameters in the plot scales for all the plots we do that here

# convert the heatmap to a grob
heatmapGrob <- ggplotGrob(heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)))

# convert the dendrogram to a grob
# note! we flipped the axis above so the x-axis is now displayed as what we would think of as the y-axis
geneDendrogramGrob <- ggplotGrob(geneDendrogram + scale_x_discrete(expand=c(0, 0)))

# we already have a sample Dendrogram, but here it is again
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0)))

# we already have our sample clinical plot but here it is again
sampleClinicalGrob <- sampleClinicalGrob

################################################################################
######### Step 4: align the gene dendrograms to match the heatmap ##############

# check that the both the heatmap and gene dendrogram have the same number of vertical elements
#length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)

# make sure every height between the two grobs is the same
maxHeight <- unit.pmax(geneDendrogramGrob$heights, heatmapGrob$heights)
geneDendrogramGrob$heights <- as.list(maxHeight)
heatmapGrob$heights <- as.list(maxHeight)

################################################################################
# Step 4b: we have a new heatmap so we need to re-align the horizontal elements #

# repeat the steps in the tutorial

# check the widths of each grob
# sampleDendrogramGrob$widths
# heatmapGrob$widths
# sampleClinicalGrob$widths

# add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

################################################################################
############### Step 5: create a blank panel ###################################

# we can use grid graphics for this
blankPanel <- grid.rect(gp=gpar(col="white"))

################################################################################
############### Step 6: Arrange the final result ###############################

# arrange all the plots together
finalGrob_v2 <- arrangeGrob(blankPanel, sampleDendrogramGrob, blankPanel, sampleClinicalGrob, geneDendrogramGrob, heatmapGrob, ncol=2, nrow=3, widths=c(1,5), heights=c(2,.8,6))

# draw the final result
grid.draw(finalGrob_v2)

################################################################################
############### Pathway and Enrichment Analysis using DEG ###############################

## @knitr GO_Pathway_Analysis_MV

# we use ggplot2 to add x axis labels (ex: ridgeplot)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2
df = read.csv("deseq2_day0_vs_day7_pathways.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(
  geneList=gene_list,
  ont = "All",
  OrgDb = organism,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea")

dotplot(gse, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

## @knitr EMapPlot_MV
x1 <- pairwise_termsim(gse)
emapplot(x1)

## @knitr cMapPlot_MV
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(x1, categorySize="pvalue", foldChange=gene_list, showCategory = 5)

## @knitr RidglePlot_MV
ridgeplot(gse) + labs(x = "enrichment distribution")

## @knitr GSEA_TopPathway_MV
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

## @knitr Published_Articles_MV
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2023, proportion=FALSE)

## @knitr KEGG_Pathways
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$Gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "human"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


## @knitr KEGG_TopPathway_MV
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

## @knitr KEGG_Pathway_Enrichment
x2 <- pairwise_termsim(kk2)
emapplot(x2)

## @knitr cnet_KEGG_Pathway_Enrichment
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(x2, categorySize="pvalue", foldChange=gene_list)

## @knitr RidgePlot_KEGG_Pathway_Enrichment
ridgeplot(kk2) + labs(x = "enrichment distribution")

## @knitr Top_KEGG_Pathway_Enrichment
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

## @knitr Measles_KEGG_Pathway
# Produce the native KEGG plot (PNG)
Measles <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05162", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
Measles <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05162", species = kegg_organism, kegg.native = F)


#########################################################################################################################################################
########################### Now let's repeat the analysis using Toca Time ##############################################################################
# Extract differential expression results
# For "Time" perform control vs treated comparison

## @knitr DEGSeq2_Analysis_Toca2
deseq2Results <- results(deseq2Data, contrast=c("Time",  "D56", "D0"))

# View summary of results
#summary(deseq2Results)


# Using DEseq2 built in method
# plotMA(deseq2Results)
# Load libraries
# install.packages(c("ggplot2", "scales", "viridis"))


# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# Examine this data frame
#head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .05, "Significant", NA)

# Plot the results similar to DEseq2
MAPlot1 <- ggplot(deseq2ResDF,  aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw() + ggtitle("Control vs MV Treatment")

# Let's add some more detail
MAPlot2 <- ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() + geom_density_2d(colour="black", size=2) + ggtitle("Control vs MV Treatment")

## @knitr MAPlots2
MAPlot1 + MAPlot2

# Add rectangle around labels
deseq2ResDF2 <- rownames_to_column(deseq2ResDF, var = "gene")
# Add rectangle around labels

## @knitr MAPlots_Significant_Genes2
ggmaplot(deseq2ResDF2, main = expression("Day 0" %->% "Day 56"),
         fdr = 0.05, fc = 1, size = 1.2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(deseq2ResDF2$gene),
         legend = "top", top = 20,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())


## @knitr VolcanoPlot2
EnhancedVolcano(deseq2ResDF2,
                lab = as.character(deseq2ResDF2$gene),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                title = "Volcano Plot DEG (Day 0 vs Day 56)",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 2,
                labSize = 4,
                colAlpha = 1,
                
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 2,
                drawConnectors = TRUE,
                widthConnectors = 0.75)



# Extract counts for the gene otop2
otop2Counts <- plotCounts(deseq2Data, gene="HOXA9", intgroup=c("Time", "Sample"), returnData=TRUE)

# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
## @knitr PlotCounts2
ggplot(otop2Counts, aes(x=Sample, y=count, colour=Time, group=Time)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("HOXA9")

deseq2ResDF["HOXA9",]
rawCounts["HOXA9",]
normals=row.names(sampleData[sampleData[,"Time"]=="D0",])
primaries=row.names(sampleData[sampleData[,"Time"]=="D56",])
rawCounts["HOXA9",normals]
rawCounts["HOXA9",primaries]

## @knitr HeatmapTopGenes2
# Transform count data using the variance stablilizing transform
deseq2VST <- vst(deseq2Data)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
#head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 2,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Convert the VST counts to long format for ggplot2


# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

#head(deseq2VST_wide)
#head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
# heatmap

# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="ward.D2")
clusterSample <- hclust(distanceSample, method="ward.D2")

# Construct a dendogram for samples
# install.packages("ggdendro")

sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
#heatmap

# Combine the dendrogram and the heatmap
# install.packages("gridExtra")

#grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))


# Load in libraries necessary for modifying plots
#install.packages("gtable")


# Modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Check the widths of each grob
# sampleDendrogramGrob$widths
# heatmapGrob$widths

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
#grid.draw(finalGrob)


# Re-order the sample data to match the clustering we did
sampleData_v2$Mixture <- factor(sampleData_v2$Mixture, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52", "#00AFBB", "#E7B800", "#FC4E07", "green", "red", "navy")
sampleClinical <- ggplot(sampleData_v2, aes(x=Mixture, y=1, fill=Days)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Time", values=colours) + theme_void()

# Convert the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

# Make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# Arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
#grid.draw(finalGrob)


###############################################################################
################# Step 1: create dendrogram for genes ##########################

# we already did the clustering for genes in the tutorial, get the data to make a dendrogram with ggplot
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))

# construct the dendrogram in ggplot
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) + theme_dendro()

################################################################################
################# Step 2: Re-arrange the heatmap cells #########################

# re-factor genes for ggplot2
deseq2VST$Gene <- factor(deseq2VST$Gene, levels=clusterGene$labels[clusterGene$order])

# recreate the heatmap with this new factoring
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())

################################################################################
################# Step 3: convert to everything to grobs #######################

# note! before this step as mentioned you might need to alter the expand parameters in the plot scales for all the plots we do that here

# convert the heatmap to a grob
heatmapGrob <- ggplotGrob(heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)))

# convert the dendrogram to a grob
# note! we flipped the axis above so the x-axis is now displayed as what we would think of as the y-axis
geneDendrogramGrob <- ggplotGrob(geneDendrogram + scale_x_discrete(expand=c(0, 0)))

# we already have a sample Dendrogram, but here it is again
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0)))

# we already have our sample clinical plot but here it is again
sampleClinicalGrob <- sampleClinicalGrob

################################################################################
######### Step 4: align the gene dendrograms to match the heatmap ##############

# check that the both the heatmap and gene dendrogram have the same number of vertical elements
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)

# make sure every height between the two grobs is the same
maxHeight <- unit.pmax(geneDendrogramGrob$heights, heatmapGrob$heights)
geneDendrogramGrob$heights <- as.list(maxHeight)
heatmapGrob$heights <- as.list(maxHeight)

################################################################################
# Step 4b: we have a new heatmap so we need to re-align the horizontal elements #

# repeat the steps in the tutorial

# check the widths of each grob
# sampleDendrogramGrob$widths
# heatmapGrob$widths
# sampleClinicalGrob$widths

# add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

################################################################################
############### Step 5: create a blank panel ###################################

# we can use grid graphics for this
blankPanel <- grid.rect(gp=gpar(col="white"))

################################################################################
############### Step 6: Arrange the final result ###############################

# arrange all the plots together
finalGrob_v2 <- arrangeGrob(blankPanel, sampleDendrogramGrob, blankPanel, sampleClinicalGrob, geneDendrogramGrob, heatmapGrob, ncol=2, nrow=3, widths=c(1,5), heights=c(2,.8,6))

# draw the final result
grid.draw(finalGrob_v2)

################################################################################
############### Pathway and Enrichment Analysis using DEG ###############################

## @knitr GO_Pathway_Analysis_MV2

# we use ggplot2 to add x axis labels (ex: ridgeplot)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2
df = read.csv("deseq2_D0_vs_D56_pathways.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(
  geneList=gene_list,
  ont = "All",
  OrgDb = organism,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea")


dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

## @knitr EMapPlot_MV2
x1 <- pairwise_termsim(gse)
emapplot(x1)

## @knitr cMapPlot_MV2
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(x1, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

## @knitr RidglePlot_MV2
ridgeplot(gse) + labs(x = "enrichment distribution")

## @knitr GSEA_TopPathway_MV2
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

## @knitr Published_Articles_MV2
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2023, proportion=FALSE)

## @knitr KEGG_Pathways2
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$Gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "human"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


## @knitr KEGG_TopPathway_MV2
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

## @knitr KEGG_Pathway_Enrichment2
x2 <- pairwise_termsim(kk2)
emapplot(x2)

## @knitr cnet_KEGG_Pathway_Enrichment2
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(x2, categorySize="pvalue", foldChange=gene_list)

## @knitr RidgePlot_KEGG_Pathway_Enrichment2
ridgeplot(kk2) + labs(x = "enrichment distribution")

## @knitr Top_KEGG_Pathway_Enrichment2
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

## @knitr Measles_KEGG_Pathway2
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05235", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05235", species = kegg_organism, kegg.native = F)



#####################################################################################################################################
############### Survival analysis comparing MV and Toca Time ######################################################################

## @knitr Survival_Analysis_MV_Infected_Time
#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/01_Input_Data/Mouse/")

# mouseMV <- read.csv(file="Mouse_OS_data_MTime.csv")
# sfit <- survfit(Surv(time, status)~1, data=mouseMV)
# sfit <- survfit(Surv(time, status)~Group_Responses, data=mouseMV)
# 
# 
# 
# ggsurvplot(
#   sfit,                     # survfit object with calculated statistics.
#   data = mouseMV,             # data used to fit survival curves.
#   risk.table = TRUE,       # show risk table.
#   pval = TRUE,             # show p-value of log-rank test.
#   conf.int = FALSE,         # show confidence intervals for 
#   # point estimates of survival curves.
#   palette = c("green", "red", "#E7B800", "#2E9FDF"),
#   title="Kaplan-Meier Curve for Mouse MV Medullablastoma Survival", 
#   xlim = c(0,150),         # present narrower X axis, but not affect
#   # survival estimates.
#   xlab = "Time in days",   # customize X axis label.
#   break.time.by = 10,     # break X axis in time intervals by 500.
#   ggtheme = theme_light(), # customize plot and risk table with a theme.
#   risk.table.y.text.col = T,# colour risk table text annotations.
#   risk.table.height = 0.25, # the height of the risk table
#   risk.table.y.text = FALSE,# show bars instead of names in text annotations
#   # in legend of risk table.
#   ncensor.plot = TRUE,      # plot the number of censored subjects at time t
#   ncensor.plot.height = 0.25,
#   conf.int.style = "step",  # customize style of confidence intervals
#   surv.median.line = "hv",  # add the median survival pointer.
#   # change legend labels.
# )
# 
# 
# 
# 
# ## @knitr ForestPlot_MV_Infected_Time
# 
# forest_model(coxph(Surv(time, status) ~ ., mouseMV)) +
#   ggtitle("Mouse Measles Time Survival Forest Plot")
# 
# 
# ##############################################################################################
# 
# ## @knitr Survival_Analysis_Toca_Infected_Time
# mouseToca <- read.csv(file="Mouse_OS_data_Toca.csv")
# sfit1 <- survfit(Surv(time, status)~1, data=mouseToca)
# sfit1 <- survfit(Surv(time, status)~Group_Responses, data=mouseToca)
# 
# 
# ggsurvplot(
#   sfit1,                     # survfit object with calculated statistics.
#   data = mouseToca,             # data used to fit survival curves.
#   risk.table = TRUE,       # show risk table.
#   pval = TRUE,             # show p-value of log-rank test.
#   conf.int = FALSE,         # show confidence intervals for 
#   # point estimates of survival curves.
#   palette = c("green", "red", "#E7B800", "#2E9FDF"),
#   title="Kaplan-Meier Curve for Mouse Toca Medullablastoma Survival", 
#   xlim = c(0,150),         # present narrower X axis, but not affect
#   # survival estimates.
#   xlab = "Time in days",   # customize X axis label.
#   break.time.by = 10,     # break X axis in time intervals by 500.
#   ggtheme = theme_light(), # customize plot and risk table with a theme.
#   risk.table.y.text.col = T,# colour risk table text annotations.
#   risk.table.height = 0.25, # the height of the risk table
#   risk.table.y.text = FALSE,# show bars instead of names in text annotations
#   # in legend of risk table.
#   ncensor.plot = TRUE,      # plot the number of censored subjects at time t
#   ncensor.plot.height = 0.25,
#   conf.int.style = "step",  # customize style of confidence intervals
#   surv.median.line = "hv",  # add the median survival pointer.
#   # change legend labels.
# )
# 
# ## @knitr ForestPlot_Toca_Infected_Time
# 
# forest_model(coxph(Surv(time, status) ~ ., mouseToca)) +
#   ggtitle("Mouse Toca Time Survival Forest Plot") 
# 


################################################################################
############### 3D Volcano Plots using DEG ###############################

## @knitr Volcano3D_MV
#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/05_Output/Mouse/SingleR")
Raw_reads <- read.csv("Strata_combined.csv", header = TRUE, row.names = 1)
Meta_Data_Volcano3D_MV <- read.csv("Meta_Data.csv", header = TRUE)


#head(Meta_Data_Volcano3D_MV)
rownames(Meta_Data_Volcano3D_MV) <- Meta_Data_Volcano3D_MV$Mixture
keep <- c("Strata", "Days")
Meta_Data_Volcano3D_MV <- Meta_Data_Volcano3D_MV[,keep]
colnames(Meta_Data_Volcano3D_MV) <- c("Strata", "Time")
Meta_Data_Volcano3D_MV$Strata <- factor(Meta_Data_Volcano3D_MV$Strata)
#head(Meta_Data_Volcano3D_MV)


# rename the tissue types
Strata <- function(x){
  x <- switch(as.character(x), "A"="Strata A", "B"="Strata B",  "C"="Strata C")
  return(x)
}

Meta_Data_Volcano3D_MV$Strata <- unlist(lapply(Meta_Data_Volcano3D_MV$Strata, Strata))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
Meta_Data_Volcano3D_MV$Strata<- factor(Meta_Data_Volcano3D_MV$Strata, levels=c("Strata A",  "Strata B", "Strata C"))

# Create the DEseq2DataSet object
deseq2DataMV <- DESeqDataSetFromMatrix(countData=round(Raw_reads), colData=Meta_Data_Volcano3D_MV, design= ~ Strata + Time)


# initial analysis run
dds_DE <- DESeq(deseq2DataMV)
# likelihood ratio test on 'Pathotype'
dds_LRT <- DESeq(deseq2DataMV, test = "LRT", reduced=~1, parallel = TRUE) 

# create 'volc3d' class object for plotting
res <- deseq_polar(dds_DE, dds_LRT, "Strata",  pcutoff = 0.05,
                   padj.method = "BH",
                   filter_pairwise = TRUE)


# plot 3d volcano plot
## @knitr Volcano3D_MV_DEGPlot

volcano3D(res, label_rows = c("AMBRA1", "MALAT1", "ARPC1B"), 
          label_size = 20,
          colour_code_labels = F,
          label_colour = "black")


## @knitr Volcano3D_MV_DEGP_RadialPlot
labs <- c('MALAT1')
radial_ggplot(res, label_rows = labs,
              marker_size = 2.3,
              legend_size = 10) +
  theme(legend.position = "right")

## @knitr Volcano3D_MV_Boxplots
plot1 <- boxplot_trio(res,
                      value = "SLC27A4",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D28"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot2 <- boxplot_trio(res,
                      value = "LINC01204",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D28"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot3 <- boxplot_trio(res,
                      value = "CORO1A",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D28"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

ggarrange(plot1, plot2, plot3, ncol=3)

plot4 <- boxplot_trio(res,
                      value = "CD8A",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D28"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot5 <- boxplot_trio(res,
                      value = "CD4",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D28"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot6 <- boxplot_trio(res,
                      value = "CD46",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D28"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

ggarrange(plot4, plot5, plot6, ncol=3)



pbmc.markers <- FindAllMarkers(Raw_reads, min.pct = 0.1, logfc.threshold = 1, test.use = "MAST", group.by = Meta_Data_Volcano3D_MV$Strata)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)



top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)


library("clusterProfiler")
library("org.Mm.eg.db")
library("AnnotationHub")

df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

#The output of length(dfsample) returns how many clusters you have
#Here there at 9 clusters (0, 1, 2, 3, 4, 5, 6, 7 and 8)
#I'm sure there's a better way but you have to make a line like below for each cluster

dfsample$Malignant = bitr(dfsample$Malignant, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$Microglia = bitr(dfsample$Microglia, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$NK_cells = bitr(dfsample$NK_cells, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`Endothelial cells` = bitr(dfsample$Endothelial, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$Monocytes = bitr(dfsample$Monocytes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$M1_Macrophages = bitr(dfsample$M1_Macrophages, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`B cells` = bitr(dfsample$`B cells`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`CD4+_T_cells` = bitr(dfsample$`CD4+_T_cells`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$Oligodendrocytes = bitr(dfsample$Oligodendrocytes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$M2_Macrophages = bitr(dfsample$M2_Macrophages, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$M0_Macrophages = bitr(dfsample$M0_Macrophages, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$Neurons = bitr(dfsample$Neurons, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$Astrocytes = bitr(dfsample$Astrocytes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`Epithelial cells` = bitr(dfsample$`Epithelial cells`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

#do the same here, a line like below for each cluster
genelist <- list("Malignant" = dfsample$Malignant$ENTREZID, 
                 "Microglia"  = dfsample$Microglia$ENTREZID,
                 "NK_cells" = dfsample$NK_cells$ENTREZID,
                 "Endothelial cells"  = dfsample$`Endothelial cells`$ENTREZID,
                 "Monocytes"  = dfsample$Monocytes$ENTREZID,
                 "M1_Macrophages" = dfsample$M1_Macrophages$ENTREZID,
                 "B cells"   = dfsample$`B cells`$ENTREZID,
                 "CD4+_T_cells" = dfsample$`CD4+_T_cells`$ENTREZID,
                 "Oligodendrocytes" = dfsample$Oligodendrocytes$ENTREZID,
                 "M2_Macrophages" = dfsample$M2_Macrophages$ENTREZID,
                 "M0_Macrophages"  = dfsample$M0_Macrophages$ENTREZID,
                 "Neurons"  = dfsample$Neurons$ENTREZID,
                 "Astrocytes"  = dfsample$Astrocytes$ENTREZID,
                 "Epithelial cells"  = dfsample$`Epithelial cells`$ENTREZID
)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot, showCategory=3)


## @knitr Volcano2D_MV_DEGPlots

DEGdf <- mutate(Volcano3D_Results, sig=ifelse(Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$D7<=-0.5 | Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$D7>=0.5, "Significant", "Not Significant")) #Will have different colors depending on significance
head(DEGdf)

volc = ggplot(DEGdf, aes(D7, -log10(pvalue))) + #log2Foldchange vs FDR
  geom_point(aes(col=sig)) + # the points are colored by significance
  scale_color_manual(values=c("blue", "red")) + 
  ggtitle("Volcano Plot Time in Days ") #e.g. 'Title description'
volc+geom_text_repel(DEGdf=head(DEGdf, 5), aes(label=gene)) #adding text for the top 5 genes


## @knitr Volcano2D_MV_GOPathways

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2


# we want the log2 fold change 
original_gene_list <- Volcano3D_Results$D7

# name the vector
names(original_gene_list) <- Volcano3D_Results$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(
  geneList=gene_list,
  ont = "All",
  OrgDb = organism,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

## @knitr Volcano2D_MV_GO_Enrichment
x1 <- pairwise_termsim(gse)
emapplot(x1)


## @knitr MeaslesDEGVolcano3D
# Extract differential expression results
# Coerce to a data frame
Fourweeks <- as.data.frame(res@df$scaled)
Fourweeks <- rownames_to_column(Fourweeks, var = "gene")
Fourweeks$significant <- ifelse(Fourweeks$pvalue < .05, "Significant", NA)
Fourweeks<-na.omit(Fourweeks)

## @knitr Print_TablesMV
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(Fourweeks, server = TRUE,
                                              class = "compact",
                                              filter="top",
                                              rownames = FALSE,
                                              colnames = c("gene", "D0", "D7", "D28", "x", "y", "r", "angle", "z", "pvalue", "col", "col", "significant"),
                                              extensions = c('Buttons'),
                                              options = list(
                                                pageLength = 15,
                                                dom = 'Bfrtip',
                                                buttons = c('excel','csv','pdf','copy')
                                              ))))


###############################################################################################################################################################
###################################### And now let's do the same for Toca Time ###############################################################################

## @knitr Volcano3D_MV2
#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/05_Output/Mouse/SingleR")
Raw_reads <- read.csv("Human_Bulk_RNA_seq_Raw_Read_Days_0_7_56.csv", header = TRUE, row.names = 1)
Meta_Data_Volcano3D_MV <- read.csv("Meta_Data_Days_0_7_56.csv", header = TRUE)

#head(Meta_Data_Volcano3D_MV)
rownames(Meta_Data_Volcano3D_MV) <- Meta_Data_Volcano3D_MV$Mixture
keep <- c("Sample_Number", "Days")
Meta_Data_Volcano3D_MV <- Meta_Data_Volcano3D_MV[,keep]
colnames(Meta_Data_Volcano3D_MV) <- c("Sample", "Time")
Meta_Data_Volcano3D_MV$Time <- factor(Meta_Data_Volcano3D_MV$Time)
#head(Meta_Data_Volcano3D_MV)


# rename the tissue types
Time_Infection <- function(x){
  x <- switch(as.character(x), "D0"="D0", "D7"="D7",
              "D56"="D56")
  return(x)
}

Meta_Data_Volcano3D_MV$Time <- unlist(lapply(Meta_Data_Volcano3D_MV$Time, Time_Infection))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
Meta_Data_Volcano3D_MV$Time<- factor(Meta_Data_Volcano3D_MV$Time, levels=c("D0",  "D7", "D56"))

# Create the DEseq2DataSet object
deseq2DataToca <- DESeqDataSetFromMatrix(countData=round(Raw_reads), colData=Meta_Data_Volcano3D_MV, design= ~ Sample + Time)


# initial analysis run
dds_DE <- DESeq(deseq2DataToca)
# likelihood ratio test on 'Pathotype'
dds_LRT <- DESeq(deseq2DataToca, test = "LRT", reduced=~1, parallel = TRUE) 

# create 'volc3d' class object for plotting
res <- deseq_polar(dds_DE, dds_LRT, "Time", pcutoff = 0.05,
                   padj.method = "none",
                   filter_pairwise = FALSE)

## @knitr Volcano3D_MV_DEGPlot2
# plot 3d volcano plot
volcano3D(res)


## @knitr Volcano3D_MV_DEGP_RadialPlot2
radial_ggplot(res,
              marker_size = 2.3,
              legend_size = 10) +
  theme(legend.position = "right")


## @knitr Volcano3D_MV_Boxplots2
plot1 <- boxplot_trio(res,
                      value = "DOCK4",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D56"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot2 <- boxplot_trio(res,
                      value = "RP11-734I18.1",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D56"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot3 <- boxplot_trio(res,
                      value = "POTEG",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D56"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

ggarrange(plot1, plot2, plot3, ncol=3)

plot4 <- boxplot_trio(res,
                      value = "CD8A",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D56"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot5 <- boxplot_trio(res,
                      value = "CD4",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D56"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot6 <- boxplot_trio(res,
                      value = "CD46",
                      test = "polar_pvalue",
                      levels_order = c("D0",  "D7", "D56"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

ggarrange(plot4, plot5, plot6, ncol=3)




Volcano3D_Results <- rownames_to_column(res@df$scaled, var = "gene")

# Syntax
Volcano3D_Results <- na.omit(Volcano3D_Results)


## @knitr Volcano2D_MV_DEGPlots2
DEGdf <- mutate(Volcano3D_Results, sig=ifelse(Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$D56<=-0.5 | Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$D56>=0.5, "Significant", "Not Significant")) #Will have different colors depending on significance
#head(DEGdf)

volc = ggplot(DEGdf, aes(D56, -log10(pvalue))) + #log2Foldchange vs FDR
  geom_point(aes(col=sig)) + # the points are colored by significance
  scale_color_manual(values=c("blue", "red")) + 
  ggtitle("Volcano Plot Post Toca Treatment ") #e.g. 'Title description'
volc+geom_text_repel(DEGdf=head(DEGdf, 5), aes(label=gene)) #adding text for the top 5 genes


## @knitr Volcano2D_MV_GOPathways2
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2


# we want the log2 fold change 
original_gene_list <- Volcano3D_Results$D56

# name the vector
names(original_gene_list) <- Volcano3D_Results$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(
  geneList=gene_list,
  ont = "All",
  OrgDb = organism,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

## @knitr Volcano2D_MV_GO_Enrichment2
x1 <- pairwise_termsim(gse)
emapplot(x1)

## @knitr TocaDEGVolcano3D
# Extract differential expression results
# Coerce to a data frame
Eightweeks <- as.data.frame(res@df$scaled)
Eightweeks <- rownames_to_column(Eightweeks, var = "gene")
Eightweeks$significant <- ifelse(Eightweeks$pvalue < .05, "Significant", NA)
Eightweeks<-na.omit(Eightweeks)

## @knitr Print_TablesToca
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(Eightweeks, server = TRUE,
                                              class = "compact",
                                              filter="top",
                                              rownames = FALSE,
                                              colnames = c("gene", "D0", "D7", "D56", "x", "y", "r", "angle", "z", "pvalue", "col", "col", "significant"),
                                              extensions = c('Buttons'),
                                              options = list(
                                                pageLength = 15,
                                                dom = 'Bfrtip',
                                                buttons = c('excel','csv','pdf','copy')
                                              ))))




venn <- read.csv("Venn_Diagrams_28_Days_Up_Regulated_Genes.csv")
library("ggVennDiagram")
library("ggvenn")
# Default plot
# Change category names
# Change the gradient fill color
venn <- as.list(venn)
#names(venn) <- c("Strata A most variable genes","Strata B DEG (2 weeks)","Strata B DEG (4 weeks)")
V1 <- ggvenn(
  venn, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

V1 + # Add title to ggplot2 venn diagram
  ggtitle("Top 2000 Genes Overlapping Across Stratas Bulk RNA")
# Remove labels background color

V2 <- ggVennDiagram(
  venn, label_alpha = 0,
  category.names = c("Strata A","Strata B","Strata C")
) +
  ggplot2::scale_fill_gradient2(low="blue", mid = "white", high = "red", midpoint = 2000,  space = "Lab",
                                na.value = "grey50",
                                guide = "colourbar",
                                aesthetics = "fill")

V2 + # Add title to ggplot2 venn diagram
  ggtitle("Up-regulated genes overlapping across strata PNOC human bulk-RNA Day 28")


de1_symbols <- read.csv("DEG_B_C_D28")
de1_symbols <- B_C_D28
keyvals <- ifelse(de1_symbols$log2FoldChange < -2, 'blue', ifelse(de1_symbols$log2FoldChange > 2, 'red', 'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'blue'] <- 'Down-Regulated'
names(keyvals)[keyvals == 'black'] <- 'NS'
names(keyvals)[keyvals == 'red'] <- 'Up-Regulated'
EnhancedVolcano(de1_symbols, lab = as.character(de1_symbols$gene), x = 'log2FoldChange', y = 'padj', selectLab = as.character(de1_symbols$gene)[which(names(keyvals) %in% c('Down-Regulated', 'Up-Regulated'))], xlab = bquote(~log[2]~ 'fold change'), title = "B vs C (Day 28)", pCutoff = 0.05, subtitle = 'Cutoff values (dashed line) at padj=0.05 log2FC=2', FCcutoff = 2, pointSize = 1.5, labSize = 4,  colCustom = keyvals, legendPosition = 'right', legendLabSize = 13, colAlpha = 1, legendIconSize = 5.0, drawConnectors = TRUE, widthConnectors = 0.75, gridlines.major = FALSE, gridlines.minor = FALSE, border = 'partial', 
                max.overlaps = 20, borderWidth = 1.5, borderColour = 'black')



pbmc.markers <- read.csv(file="Master_DEG_clusterProfiler_human_MV_timepoint_analysis_both.csv")
#top100 <- pbmc.markers[which(abs(pbmc.markers$log2FoldChange) > 2 & pbmc.markers$padj < 0.05),]


library("clusterProfiler")
library("org.Hs.eg.db")
library("AnnotationHub")

data(geneList, package="DOSE")
head(geneList)
gene <- names(geneList)[abs(geneList) > 2]
head(gene)

pbmc.markers <- data.frame(Entrez=names(df$gene), FC=df$log2FoldChange)
pbmc.markers <- pbmc.markers[abs(pbmc.markers$log2FoldChange) > 2,]
pbmc.markers$group <- "upregulated"
pbmc.markers$group[pbmc.markers$FC < 2] <- "downregulated"
pbmc.markers$othergroup <- "A"
pbmc.markers$othergroup[abs(pbmc.markers$log2FoldChange) > 2] <- "B"

formula_res <- compareCluster(pbmc.markers$gene~pbmc.markers$strata, data=pbmc.markers, fun="enrichKEGG")

head(formula_res)


df <- pbmc.markers[,1:10]
dfsample <- split(df$gene,df$cluster, df$strata)
length(dfsample)

#The output of length(dfsample) returns how many clusters you have
#Here there at 9 clusters (0, 1, 2, 3, 4, 5, 6, 7 and 8)
#I'm sure there's a better way but you have to make a line like below for each cluster

dfsample$A_D0_D14 = bitr(dfsample$A_D0_D14, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$A_D0_D28 = bitr(dfsample$A_D0_D28, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$B_D0_D1 = bitr(dfsample$B_D0_D1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$B_D0_D14 = bitr(dfsample$B_D0_D14, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$B_D0_D28 = bitr(dfsample$B_D0_D28, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$B_D0_D4 = bitr(dfsample$B_D0_D4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$B_D0_D7 = bitr(dfsample$B_D0_D7, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$C_D0_D28 = bitr(dfsample$C_D0_D28, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$B_C_D28 = bitr(dfsample$B_C_D28, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$C_D0_D7 = bitr(dfsample$C_D0_D7, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#do the same here, a line like below for each cluster
genelist <- list("A_D0_D14" = dfsample$A_D0_D14$ENTREZID, 
                 "A_D0_D28"  = dfsample$A_D0_D28$ENTREZID,
                 "B_D0_D1" = dfsample$B_D0_D1$ENTREZID,
                 "B_D0_D14" = dfsample$B_D0_D14$ENTREZID, 
                 "B_D0_D28"  = dfsample$B_D0_D28$ENTREZID,
                 "B_D0_D4" = dfsample$B_D0_D4$ENTREZID,
                 "B_D0_D7" = dfsample$B_D0_D7$ENTREZID, 
                 "C_D0_D28"  = dfsample$C_D0_D28$ENTREZID,
                 "B_C_D28" = dfsample$B_C_D28$ENTREZID,
                 "C_D0_D7" = dfsample$C_D0_D7$ENTREZID)



GOclusterplot <- compareCluster(geneCluster = genelist, ont = "All", fun="enrichGO", OrgDb="org.Hs.eg.db")

dotplot(GOclusterplot, showCategory=3)
resTimeGOTable = as.data.frame(GOclusterplot)
head(resTimeGOTable)


n_11 = resTimeGOTable$Count
n_10 = 983  # length(intersect(resTimeGO@gene, resTimeGO@universe))
n_01 = as.numeric(gsub("/.*$", "", resTimeGOTable$BgRatio))
n = 28943  # length(resTimeGO@universe)

resTimeGOTable$DE_Ratio = n_11/n_01
resTimeGOTable$GS_size = n_01  # size of gene sets
resTimeGOTable$log2_Enrichment = log( (n_11/n_10)/(n_01/n) )

hyper_mean = n_01*n_10/n

n_02 = n - n_01
n_20 = n - n_10
hyper_var = n_01*n_10/n * n_20*n_02/n/(n-1)
resTimeGOTable$zScore = (n_11 - hyper_mean)/sqrt(hyper_var)

library(ggplot2)
ggplot(resTimeGOTable[1:20, ], 
       aes(x = log2_Enrichment, y = factor(Description, levels = rev(Description)), 
           fill = DE_Ratio)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = log2_Enrichment, 
                label = sprintf("%.2e", p.adjust)), hjust = 1, col = "white") +
  ylab("")


ggplot(resTimeGOTable[1:20, ], 
       aes(x = zScore, y = factor(Description, levels = rev(Description)), 
           col = DE_Ratio, size = Count)) +
  geom_point() +
  ylab("")

ggplot(resTimeGOTable, 
       aes(x = log2_Enrichment, y = -log10(p.adjust), 
           color = DE_Ratio, size = GS_size)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.01), lty = 2, col = "#444444") +
  geom_vline(xintercept = 1.5, lty = 2, col = "#444444")


# up-regulated genes
timeDEup <- as.data.frame(subset(pbmc.markers, padj < 0.05 & log2FoldChange > log2(2)))
#timeDEupGenes <- rownames(timeDEup$gene)
resTimeGOup = enrichGO(gene = timeDEup$gene, 
                       keyType = "SYMBOL",
                       ont = "All", 
                       OrgDb = org.Hs.eg.db,
                       #universe = rownames(se),
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
barplot(resTimeGOup, showCategory=20) 
resTimeGOupTable = as.data.frame(resTimeGOup)
n_11 = resTimeGOupTable$Count
n_10 = length(intersect(resTimeGOup@gene, resTimeGOup@universe))
n_01 = as.numeric(gsub("/.*$", "", resTimeGOupTable$BgRatio))
n = length(resTimeGOup@universe)
resTimeGOupTable$log2_Enrichment = log( (n_11/n_10)/(n_01/n) )

# down-regulated genes
timeDEdown <- as.data.frame(subset(pbmc.markers, padj < 0.05 & log2FoldChange < -log2(2)))
#timeDEdownGenes <- rownames(timeDEdown$gene)
resTimeGOdown = enrichGO(gene = timeDEdown$gene, 
                         keyType = "SYMBOL",
                         ont = "All", 
                         OrgDb = org.Hs.eg.db,
                         #universe = rownames(se),
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
barplot(resTimeGOdown, showCategory=20) 
resTimeGOdownTable = as.data.frame(resTimeGOdown)
n_11 = resTimeGOdownTable$Count
n_10 = length(intersect(resTimeGOdown@gene, resTimeGOdown@universe))
n_01 = as.numeric(gsub("/.*$", "", resTimeGOdownTable$BgRatio))
n = length(resTimeGOdown@universe)
resTimeGOdownTable$log2_Enrichment = log( (n_11/n_10)/(n_01/n) )

# The name of the 3rd term is too long, we wrap it into two lines.
resTimeGOupTable[3, "Description"] = paste(strwrap(resTimeGOupTable[3, "Description"]), collapse = "\n")

direction = c(rep("up-regulated genes", 20), rep("down-regulated genes", 20))
ggplot(rbind(resTimeGOupTable[1:20, ],
             resTimeGOdownTable[1:20, ]),
       aes(x = log2_Enrichment, y =factor(Description, levels = rev(unique(Description)),  ordered=TRUE), 
           fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("up-regulated genes" = "red", "down-regulated genes" = "blue")) +
  geom_text(aes(x = log2_Enrichment, 
                label = sprintf("%.2e", p.adjust)), hjust = 1, col = "black") +
  ylab("")

GO_ID = resTimeGOTable$ID[resTimeGOTable$p.adjust < 0.05]
library(simplifyEnrichment)
set.seed(888)
mat = GO_similarity(GO_ID, ont = "MF", db = 'org.Hs.eg.db', measure = "Rel")
df = simplifyGO(mat)
compare_clustering_methods(mat, plot_type = "heatmap")



set.seed(123)
compare_clustering_methods(mat)
set.seed(123)
compare_clustering_methods(mat, plot_type = "heatmap")
clt = cmp_make_clusters(mat) # just a list of cluster labels
cmp_make_plot(mat, clt)
cmp_make_plot(mat, clt, plot_type = "heatmap")

download.file("https://jokergoo.github.io/cola_examples/TCGA_GBM/TCGA_GBM_subgroup.rds", 
              destfile = "TCGA_GBM_subgroup.rds", quiet = TRUE)
rl = readRDS("TCGA_GBM_subgroup.rds")
file.remove("TCGA_GBM_subgroup.rds")
library(cola)
res = rl["ATC:skmeans"]


set.seed(123)
get_signatures(res, k = 4, row_km = 4)



library(RColorBrewer)
library(DOSE)
library(enrichplot)
gene <- names(top100$gene)
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
de <- names(genelist)[1:50]

hm.palette <- colorRampPalette(c("blue", "white", "red"))
ego <- enrichGO(de,
                OrgDb = "org.Hs.eg.db",
                ont="ALL",
                readable=TRUE)

ego2 <- simplify(ego)

final_product_GO_BP_Over_rep <- heatplot(ego2, foldChange = geneList, showCategory = 20)+
  #ggplot2::coord_flip()+
  ggplot2::scale_fill_gradientn(colours = hm.palette(100))+
  ggplot2::ylab('Annotations of Biological Processes')+
  ggplot2::xlab('Gene Symbols')+
  ggplot2::ggtitle('GO Over-representation for Biological Processes')+
  ggplot2::theme(panel.background=element_rect(fill="white", colour="white"))+
  ggplot2::guides(fill=guide_legend(title="Fold Change"))

final_product_GO_BP_Over_rep



