###############################################################################
################# CIBERSORTx cell type quantification ##########################
# load libraries in R
## @knitr Load_libraries
setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/01_Input_Data/Mouse/")
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

## @knitr loadData
CIBERSORTx_SingleR_matrix <- read.csv("CIBERSORTx_Corrected_CellTypes_Martrix.csv", header = T, row.names = 1)
CIBERSORTx_SingleR_matrix <- t(CIBERSORTx_SingleR_matrix)
CIBERSORTx_SingleR_matrix <- as.data.frame(CIBERSORTx_SingleR_matrix)
CIBERSORTx_SingleR_matrix <- rownames_to_column(CIBERSORTx_SingleR_matrix, var = "CellType")


## @knitr CIBERSORTx_PCA
iris.pca <- PCA(CIBERSORTx_SingleR_matrix[,-1], graph = FALSE)

fviz_pca_ind(iris.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = CIBERSORTx_SingleR_matrix$CellType, 
             col.ind = "black", 
             palette = c("#0062B4FF", "#FFBF47", "#FF2700FF", "mediumseagreen", "sienna4", "violet", "violetred", "violetred4", "#CC313D", "blue", "limegreen", "#3B9AB2","deeppink", "#BF812D"), 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Cell Types") +
  ggtitle("2D PCA-plot from CIBERSORTx Signature Matrix") +
  theme(plot.title = element_text(hjust = 0.5))

CIBERSORTx_SingleR_matrix <- read.csv("CIBERSORTx_Corrected_CellTypes_Martrix.csv", header = T, row.names = 1)

## @knitr CIBERSORTx_Heatmap
data <- CIBERSORTx_SingleR_matrix
data <- as.matrix(data)
library("RColorBrewer")
cellType = c(Malignant = "#CC313D", Neurons = "#3B9AB2", Endothelial =  "mediumseagreen", Monocytes = "limegreen", Microglia = "blue", M2_Macrophages = "violetred4", M1_Macrophages = "violetred",  M0_Macrophages = "violet",  CD4 = "#FF2700FF", B = "#FFBF47", Astrocytes = "#0062B4FF", Oligodendrocytes = "#BF812D", NK_cells = "deeppink", Epithelial = "sienna4")
col<- colorRampPalette(c("cyan", "black", "red"))(256)
heatmap(data, ColSideColors = cellType, col = col)


## @knitr CIBERSORTx_Correlation_Heatmap

heatmaply_cor(
  cor(data),
  xlab = "Features",
  ylab = "Features",
  k_col = 2,
  k_row = 2
)

## @knitr CIBERSORTx_Correlation_Heatmap_pvalues
r <- cor(data)
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
p <- cor.test.p(data)

heatmaply_cor(
  r,
  node_type = "scatter",
  point_size_mat = -log10(p), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation"))


## @knitr CIBERSORTx_cell_fraction_results
CIBERSORTx_Results_SingleR <- read.csv("CIBERSORTx_Corrected_CellTypes_Results_website.csv", header = T, row.names = 1)

CIBERSORTx_Results_SingleR <- round(CIBERSORTx_Results_SingleR, digits = 2)
heatmaply(CIBERSORTx_Results_SingleR,
  cellnote = CIBERSORTx_Results_SingleR
)

theme_set(
  theme_minimal() +
    theme(legend.position = "right")
)

## @knitr CIBERSORTx_cell_fraction_barplot
CIBERSORTx_Results_SingleR <- read.csv("CIBERSORTx_Corrected_CellTypes_Results_t_test_MV.csv", header = T)

melt_data_propeller <- melt(CIBERSORTx_Results_SingleR, id=c("Mixture","Viral_Status"))

ggplot(melt_data_propeller, mapping = aes(x = Mixture, y = value, fill =  variable)) + 
  geom_bar(position= "stack", stat = "identity") +
  scale_fill_manual(labels = c("Malignant", "Neurons", "Endothelial", "Monocytes", "Microglia", "M2_Macrophages", "M1_Macrophages", "M0_Macrophages", "CD4", "B", "Astrocytes", "Oligodendrocytes", "NK_cells", "Epithelial"), values =c( "#CC313D", "#3B9AB2", "mediumseagreen", "limegreen", "blue", "violetred4",  "violetred",  "violet",  "#FF2700FF", "#FFBF47",  "#0062B4FF",  "#BF812D", "deeppink", "sienna4")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## @knitr CIBERSORTx_cell_fraction_boxplot
df <- read.csv("CIBERSORTx_Corrected_CellTypes_Results_t_test_MV.csv", header = T)
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV_NIS", "Mid_treatment_MV", "Post_treatment_MV"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV_NIS"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Mid_treatment_MV", "Post_treatment_MV"))
p1 = ggboxplot(df, x = "Viral_Status", y = c("M0_Macrophages", "M1_Macrophages", "M2_Macrophages"),  fill = "Viral_Status", legend = "none", combine = TRUE,
               bxp.errorbar = TRUE,
               bxp.errorbar.width = 0.4, palette = c("#00AFBB", "#E7B800",
                                                   "green", "red", "navy", "#009E73",
                                                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                   "#52854C", "#4E84C4", "purple", "black", "blue"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.05)  +
  ggtitle("CIBERSORTx Cell Type Comparison") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

p1


df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV_NIS", "Mid_treatment_MV", "Post_treatment_MV"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV_NIS"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Mid_treatment_MV", "Post_treatment_MV"))
p2 = ggboxplot(df, x = "Viral_Status", y = c("NK_cells", "CD4"), fill = "Viral_Status", legend = "none", combine = TRUE,
               bxp.errorbar = TRUE,
               bxp.errorbar.width = 0.4, palette = c("#00AFBB", "#E7B800",
                                                   "green", "red", "navy", "#009E73",
                                                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                   "#52854C", "#4E84C4", "purple", "black", "blue"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.05)  +
  ggtitle("CIBERSORTx Cell Type Comparison") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

p2


df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV_NIS", "Mid_treatment_MV", "Post_treatment_MV"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV_NIS"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Mid_treatment_MV", "Post_treatment_MV"))
p3 = ggboxplot(df, x = "Viral_Status", y = c("Monocytes", "Microglia", "B"),  fill = "Viral_Status", legend = "none", combine = TRUE,
               bxp.errorbar = TRUE,
               bxp.errorbar.width = 0.4, palette = c("#00AFBB", "#E7B800",
                                                   "green", "red", "navy", "#009E73",
                                                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                   "#52854C", "#4E84C4", "purple", "black", "blue"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.1)  +
  ggtitle("CIBERSORTx Cell Type Comparison") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

p3


df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV_NIS", "Mid_treatment_MV", "Post_treatment_MV"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV_NIS"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Mid_treatment_MV", "Post_treatment_MV"))
p4 = ggboxplot(df, x = "Viral_Status", y = c("Malignant", "Neurons"),fill = "Viral_Status", legend = "none", combine = TRUE,
               bxp.errorbar = TRUE,
               bxp.errorbar.width = 0.4, palette = c("#00AFBB", "#E7B800",
                                                   "green", "red", "navy", "#009E73",
                                                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                   "#52854C", "#4E84C4", "purple", "black", "blue"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1.2)  +
  ggtitle("CIBERSORTx Cell Type Comparison") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

p4


df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV_NIS", "Mid_treatment_MV", "Post_treatment_MV"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV_NIS"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Mid_treatment_MV", "Post_treatment_MV"))
p5 = ggboxplot(df, x = "Viral_Status", y = c("Endothelial", "Astrocytes", "Oligodendrocytes"), fill = "Viral_Status", legend = "none", combine = TRUE,
                bxp.errorbar = TRUE,
               bxp.errorbar.width = 0.4, palette = c("#00AFBB", "#E7B800",
                                                   "green", "red", "navy", "#009E73",
                                                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                   "#52854C", "#4E84C4", "purple", "black", "blue"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.10)  +
  ggtitle("CIBERSORTx Cell Type Comparison") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

p5

ggarrange(p1, p2, p3, p4, p5, nrow = 5, common.legend = FALSE)


# stat.test <- melt_data_propeller %>%
#   group_by(Viral_Status) %>%
#   tukey_hsd(value ~ variable) 


# propeller_test_SingleR <-propeller(clusters = melt_data_propeller$variable,
#                                    sample = melt_data_propeller$Mixture,
#                                    group = melt_data_propeller$Viral_Status)


# propeller_test_SingleR <- as.data.frame(propeller_test_SingleR)
# propeller_test_SingleR <- rownames_to_column(propeller_test_SingleR, var = "cellType")
#write.csv(stat.test, file="tukey_hsd_SingleR.csv")
########################################################################################################################################################
############################ Now let's take a look at the viral composition ############################################################################

#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/05_Output/Mouse/Virus/")
CIBERSORTx_Virus_matrix <- read.table("CIBERSORTx_Virus_Matrix_Normalized_Results.txt", header = T, sep = "\t", row.names = 1)
CIBERSORTx_Virus_matrix <- t(CIBERSORTx_Virus_matrix)
CIBERSORTx_Virus_matrix <- as.data.frame(CIBERSORTx_Virus_matrix)
CIBERSORTx_Virus_matrix <- rownames_to_column(CIBERSORTx_Virus_matrix, var = "CellType")

## @knitr CIBERSORTx_PCA2
res.pca <- prcomp(CIBERSORTx_Virus_matrix[,-1], scale = TRUE)
Viral_Status <- as.factor(CIBERSORTx_Virus_matrix$CellType)
fviz_pca_ind(res.pca,  geom.ind = "point", 
             pointsize = 5, 
             col.ind = Viral_Status, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Viral_Status",
             repel = TRUE
)

CIBERSORTx_Virus_matrix <- read.table("CIBERSORTx_Virus_Matrix_Normalized_Results.txt", header = T, sep = "\t", row.names = 1)

## @knitr CIBERSORTx_Heatmap2
data <- CIBERSORTx_Virus_matrix
data <- as.matrix(data)
heatmap(data, cexRow = 1, cexCol = 1, 
        main = "CIBERSORTx Virus Matrix")


## @knitr CIBERSORTx_Correlation_Heatmap2
heatmaply_cor(
  cor(data),
  xlab = "Features",
  ylab = "Features",
  k_col = 2,
  k_row = 2
)

## @knitr CIBERSORTx_Correlation_Heatmap_pvalues2
r <- cor(data)
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
p <- cor.test.p(data)


heatmaply_cor(
  r,
  node_type = "scatter",
  point_size_mat = -log10(p), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
)


## @knitr CIBERSORTx_cell_fraction_results2
CIBERSORTx_Results_Virus <- read.csv("CIBERSORTx_Virus_SMode_Normalized_Results_website.csv", header = T, row.names = 1)

heatmaply(
  CIBERSORTx_Results_Virus,
  cellnote = CIBERSORTx_Results_Virus
)


theme_set(
  theme_minimal() +
    theme(legend.position = "right")
)

## @knitr CIBERSORTx_cell_fraction_barplot2
CIBERSORTx_Results_Virus <- read.csv("CIBERSORTx_Virus_SMode_Normalized_Results_t_test.csv", header = T)

melt_data <- melt(CIBERSORTx_Results_Virus, id=c("Mixture","Viral_Status"))

ggplot(melt_data, mapping = aes(x = Mixture, y = value, fill =  variable)) + 
  geom_bar(position= "stack", stat = "identity") +
  scale_fill_manual(labels = c("Measles_virus_infected", "Non_Infected"), values = c( "#E7B800", "#00AFBB")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## @knitr CIBERSORTx_cell_fraction_boxplot2
df <- read.csv("CIBERSORTx_Virus_SMode_Normalized_Results_t_test.csv", header = T)
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV_NIS", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV_NIS"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Mid_treatment_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"), c("Mid_treatment_Toca", "Post_treatment_Toca"))
p6 = ggboxplot(df, x = "Viral_Status", y = c("Measles_virus_infected", "Non_Infected"),  add = "jitter", legend = "none", combine = TRUE,
               color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                   "green", "red", "navy", "#009E73",
                                                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                   "#52854C", "#4E84C4", "purple", "black", "blue"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5)  +
  ggtitle("CIBERSORTx Proportion of Infected Cells Comparison") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

p6


###############################################################################
################# DEGSeq2 analysis of viral infections ##########################

## @knitr PrepareData
#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/01_Input_Data/Mouse")

# Read in the raw read counts
rawCounts <- read.csv("Mouse_Bulk_RNA_seq_Raw_Reads.csv")
#head(rawCounts)

# Read in the sample mappings
sampleData <- read.csv("Meta_Data_SingleR.csv")
#head(sampleData)

# Also save a copy for later
sampleData_v2 <- sampleData

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$gene_name
sampleIndex <- grepl("X\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
#head(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
#head(sampleData)
rownames(sampleData) <- sampleData$Sample_ID
keep <- c("Sample_Number", "Viral_Status")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("Sample", "Virus")
sampleData$Virus <- factor(sampleData$Virus)
#head(sampleData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# rename the tissue types
Virus_Infection <- function(x){
  x <- switch(as.character(x), "Control_MV"="Control MV", "Control_Toca"="Control Toca",  "HI_MV_NIS"="Heat Inactivated MV", "HI_Toca"="Heat Inactivated Toca",  "Mid_treatment_MV"="Mid treatment MV", 
              "Mid_treatment_Toca"="Mid treatment Toca",  "Post_treatment_MV"="Post treatment MV", "Post_treatment_Toca"="Post treatment Toca")
  return(x)
}

sampleData$Virus <- unlist(lapply(sampleData$Virus, Virus_Infection))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$Virus <- factor(sampleData$Virus, levels=c("Control MV", "Control Toca",  "Heat Inactivated MV", "Heat Inactivated Toca",  "Mid treatment MV", 
                                                      "Mid treatment Toca",  "Post treatment MV", "Post treatment Toca"))

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ Sample + Virus)

dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])

# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]
# Register the number of cores to use

## @knitr Run_DESeq2
# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data, test="LRT", reduced=~1, parallel = TRUE)

# Pull out counts from DESeq object and normalize using the same method to perform differential gene expression
# dds <- estimateSizeFactors(deseq2Data)
# dds <- counts(dds, normalized=TRUE)

# Extract differential expression results
CvMT <- results(deseq2Data, contrast=c("Virus", "Post treatment MV", "Control MV"))
# Coerce to a data frame
CvMT <- as.data.frame(CvMT)
CvMT <- rownames_to_column(CvMT, var = "gene")
CvMT$significant <- ifelse(CvMT$padj < .05, "Significant", NA)
CvMT<-na.omit(CvMT)
write.csv(CvMT, file= "deseq2_control_vs_MV_PT.csv")
cat('\n', '<br>', '\n\n')

## @knitr Print_Tables
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(CvMT, server = TRUE,
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
CvMHI <- results(deseq2Data, contrast=c("Virus", "Heat Inactivated MV", "Control MV",))
# Coerce to a data frame
CvMHI <- as.data.frame(CvMHI)
CvMHI <- rownames_to_column(CvMHI, var = "gene")
CvMHI$significant <- ifelse(CvMHI$padj < .05, "Significant", NA)
CvMHI<-na.omit(CvMHI)

## @knitr Print_Tables1
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(CvMHI, server = TRUE,
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
CvM_MT <- results(deseq2Data, contrast=c("Virus",  "Mid treatment MV", "Control MV"))
# Coerce to a data frame
CvM_MT <- as.data.frame(CvM_MT)
CvM_MT <- rownames_to_column(CvM_MT, var = "gene")
CvM_MT$significant <- ifelse(CvM_MT$padj < .05, "Significant", NA)
CvM_MT<-na.omit(CvM_MT)

## @knitr Print_Tables2
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(CvM_MT, server = TRUE,
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


## @knitr Contrast3
# Extract differential expression results
CvPTT <- results(deseq2Data, contrast=c("Virus", "Post treatment Toca",  "Control Toca"))
# Coerce to a data frame
CvPTT <- as.data.frame(CvPTT)
CvPTT <- rownames_to_column(CvPTT, var = "gene")
CvPTT$significant <- ifelse(CvPTT$padj < .05, "Significant", NA)
CvPTT<-na.omit(CvPTT)

## @knitr Print_Tables3
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(CvPTT, server = TRUE,
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

## @knitr Contrast4
# Extract differential expression results
CvHIT <- results(deseq2Data, contrast=c("Virus", "Heat Inactivated Toca",  "Control Toca"))
# Coerce to a data frame
CvHIT <- as.data.frame(CvHIT)
CvHIT <- rownames_to_column(CvHIT, var = "gene")
CvHIT$significant <- ifelse(CvHIT$padj < .05, "Significant", NA)
CvHIT<-na.omit(CvHIT)

## @knitr Print_Tables4
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(CvHIT, server = TRUE,
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

## @knitr Contrast5
# Extract differential expression results
CvMTT <- results(deseq2Data, contrast=c("Virus", "Mid treatment Toca",  "Control Toca"))
# Coerce to a data frame
CvMTT <- as.data.frame(CvMTT)
CvMTT <- rownames_to_column(CvMTT, var = "gene")
CvMTT$significant <- ifelse(CvMTT$padj < .05, "Significant", NA)
CvMTT<-na.omit(CvMTT)

## @knitr Print_Tables5
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(CvMTT, server = TRUE,
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
deseq2Results <- results(deseq2Data, contrast=c("Virus", "Post treatment MV", "Control MV"))


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
# write these differentially expressed genes as a table and save to directory
write.csv(deseq2ResDF, file= "deseq2_control_vs_MV_PT.csv")

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
ggmaplot(deseq2ResDF2, main = expression("Control MV" %->% "Post treatment MV"),
         fdr = 0.05, fc = 1, size = 1.2,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(deseq2ResDF2$gene),
         legend = "top", top = 20,
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
                     title = "Volcano Plot DEG (Control and Post treatment MV)",
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
otop2Counts <- plotCounts(deseq2Data, gene="Arg1", intgroup=c("Virus", "Sample"), returnData=TRUE)

# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
## @knitr PlotCounts
ggplot(otop2Counts, aes(x=Sample, y=count, colour=Virus, group=Virus)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Arginase 1 (M2 Macrophage Marker)")

deseq2ResDF["Arg1",]
rawCounts["Arg1",]
normals=row.names(sampleData[sampleData[,"Virus"]=="Control MV",])
primaries=row.names(sampleData[sampleData[,"Virus"]=="Post treatment MV",])
rawCounts["Arg1",normals]
rawCounts["Arg1",primaries]

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
sampleData_v2$Sample_ID <- factor(sampleData_v2$Sample_ID, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52", "#00AFBB", "#E7B800", "#FC4E07", "green", "red", "navy")
sampleClinical <- ggplot(sampleData_v2, aes(x=Sample_ID, y=1, fill=Viral_Status)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tissue", values=colours) + theme_void()

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
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2
df = read.csv("deseq2_control_vs_MV_PT_Pathways.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$gene

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

## @knitr EMapPlot_MV
x1 <- pairwise_termsim(gse)
emapplot(x1)

## @knitr cMapPlot_MV
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(x1, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

## @knitr RidglePlot_MV
ridgeplot(gse) + labs(x = "enrichment distribution")

## @knitr GSEA_TopPathway_MV
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
#trace("gseaplot2", edit = TRUE)
gseaplot2(gse, geneSetID = 27, pvalue_table = TRUE,  ES_geom = "line", title = gse$Description[27])


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
df2 = df[df$gene %in% dedup_ids$SYMBOL,]

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
kegg_organism = "mouse"
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
pathview(gene.data=kegg_gene_list, pathway.id="mmu05235", species = kegg_organism, kegg.native = TRUE, match.data = FALSE, multi.state = TRUE, same.layer = TRUE)

# Produce a different plot (PDF) (not displayed here)
pathview(gene.data=kegg_gene_list, pathway.id="mmu05235", species = kegg_organism, kegg.native = TRUE, match.data = FALSE, multi.state = TRUE, same.layer = TRUE)


#########################################################################################################################################################
########################### Now let's repeat the analysis using Toca virus ##############################################################################
# Extract differential expression results
# For "Virus" perform control vs treated comparison

## @knitr DEGSeq2_Analysis_Toca2
deseq2Results <- results(deseq2Data, contrast=c("Virus",  "Control Toca", "Post treatment Toca"))

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
# write these differentially expressed genes as a table and save to directory
write.csv(deseq2ResDF, file= "deseq2_control_vs_Toca_PT.csv")

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
ggmaplot(deseq2ResDF2, main = expression("Control Toca" %->% "Post treatment Toca"),
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
                title = "Volcano Plot DEG (Control and Post treatment Toca)",
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
otop2Counts <- plotCounts(deseq2Data, gene="Jchain", intgroup=c("Virus", "Sample"), returnData=TRUE)

# Plot the data using ggplot2
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
## @knitr PlotCounts2
ggplot(otop2Counts, aes(x=Sample, y=count, colour=Virus, group=Virus)) + geom_point() + geom_line() + theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Jchain")

deseq2ResDF["Jchain",]
rawCounts["Jchain",]
normals=row.names(sampleData[sampleData[,"Virus"]=="Control Toca",])
primaries=row.names(sampleData[sampleData[,"Virus"]=="Post treatment Toca",])
rawCounts["Jchain",normals]
rawCounts["Jchain",primaries]

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
sampleData_v2$Sample_ID <- factor(sampleData_v2$Sample_ID, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52", "#00AFBB", "#E7B800", "#FC4E07", "green", "red", "navy")
sampleClinical <- ggplot(sampleData_v2, aes(x=Sample_ID, y=1, fill=Viral_Status)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tissue", values=colours) + theme_void()

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
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2
df = read.csv("DESeq2_control_vs_Toca_PT_Pathways.csv", header=TRUE)

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
kegg_organism = "mouse"
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
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu05168", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu05168", species = kegg_organism, kegg.native = F)



#####################################################################################################################################
############### Survival analysis comparing MV and Toca virus ######################################################################

## @knitr Survival_Analysis_MV_Infected_Virus
#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/01_Input_Data/Mouse/")

mouseMV <- read.csv(file="Mouse_OS_data_MVirus.csv")
sfit <- survfit(Surv(time, status)~1, data=mouseMV)
sfit <- survfit(Surv(time, status)~Group_Responses, data=mouseMV)



ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  data = mouseMV,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("green", "red", "#E7B800", "#2E9FDF"),
  title="Kaplan-Meier Curve for Mouse MV Medullablastoma Survival", 
  xlim = c(0,150),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  # change legend labels.
)




## @knitr ForestPlot_MV_Infected_Virus

forest_model(coxph(Surv(time, status) ~ ., mouseMV)) +
  ggtitle("Mouse Measles Virus Survival Forest Plot")


##############################################################################################

## @knitr Survival_Analysis_Toca_Infected_Virus
mouseToca <- read.csv(file="Mouse_OS_data_Toca.csv")
sfit1 <- survfit(Surv(time, status)~1, data=mouseToca)
sfit1 <- survfit(Surv(time, status)~Group_Responses, data=mouseToca)


ggsurvplot(
  sfit1,                     # survfit object with calculated statistics.
  data = mouseToca,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = c("green", "red", "#E7B800", "#2E9FDF"),
  title="Kaplan-Meier Curve for Mouse Toca Medullablastoma Survival", 
  xlim = c(0,150),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.height = 0.25, # the height of the risk table
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  # change legend labels.
)

## @knitr ForestPlot_Toca_Infected_Virus

forest_model(coxph(Surv(time, status) ~ ., mouseToca)) +
  ggtitle("Mouse Toca Virus Survival Forest Plot") 



################################################################################
############### 3D Volcano Plots using DEG ###############################

## @knitr Volcano3D_MV
#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/05_Output/Mouse/SingleR")
Raw_reads <- read.csv("Mouse_Bulk_RNA_seq_Raw_Reads_Volcano3D_MV.csv", header = TRUE, row.names = 1)
Meta_Data_Volcano3D_MV <- read.csv("Meta_Data_Volcano3D_MV.csv", header = TRUE)


#head(Meta_Data_Volcano3D_MV)
rownames(Meta_Data_Volcano3D_MV) <- Meta_Data_Volcano3D_MV$Sample_ID
keep <- c("Sample_Number", "Viral_Status")
Meta_Data_Volcano3D_MV <- Meta_Data_Volcano3D_MV[,keep]
colnames(Meta_Data_Volcano3D_MV) <- c("Sample", "Virus")
Meta_Data_Volcano3D_MV$Virus <- factor(Meta_Data_Volcano3D_MV$Virus)
#head(Meta_Data_Volcano3D_MV)


# rename the tissue types
Virus_Infection <- function(x){
  x <- switch(as.character(x), "Control_MV"="Control_MV", "Mid_treatment_MV"="Mid_treatment_MV",
              "Post_treatment_MV"="Post_treatment_MV")
  return(x)
}

Meta_Data_Volcano3D_MV$Virus <- unlist(lapply(Meta_Data_Volcano3D_MV$Virus, Virus_Infection))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
Meta_Data_Volcano3D_MV$Virus<- factor(Meta_Data_Volcano3D_MV$Virus, levels=c("Control_MV",  "Mid_treatment_MV", "Post_treatment_MV"))

# Create the DEseq2DataSet object
deseq2DataMV <- DESeqDataSetFromMatrix(countData=round(Raw_reads), colData=Meta_Data_Volcano3D_MV, design= ~ Sample + Virus)


# initial analysis run
dds_DE <- DESeq(deseq2DataMV)
# likelihood ratio test on 'Pathotype'
dds_LRT <- DESeq(deseq2DataMV, test = "LRT", reduced=~1, parallel = TRUE) 

# create 'volc3d' class object for plotting
res <- deseq_polar(dds_DE, dds_LRT, "Virus", pcutoff = 0.05,
                   padj.method = "none",
                   filter_pairwise = FALSE)


# plot 3d volcano plot
## @knitr Volcano3D_MV_DEGPlot
volcano3D(res)

## @knitr Volcano3D_MV_DEGP_RadialPlot
radial_ggplot(res,
              marker_size = 2.3,
              legend_size = 10) +
  theme(legend.position = "right")

## @knitr Volcano3D_MV_Boxplots
plot1 <- boxplot_trio(res,
                      value = "Jchain",
                      test = "polar_pvalue",
                      levels_order = c("Control_MV",  "Mid_treatment_MV", "Post_treatment_MV"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot2 <- boxplot_trio(res,
                      value = "Igkc",
                      test = "polar_pvalue",
                      levels_order = c("Control_MV",  "Mid_treatment_MV", "Post_treatment_MV"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot3 <- boxplot_trio(res,
                      value = "Arg1",
                      test = "polar_pvalue",
                      levels_order = c("Control_MV",  "Mid_treatment_MV", "Post_treatment_MV"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

ggarrange(plot1, plot2, plot3, ncol=3)

plot4 <- boxplot_trio(res,
                      value = "C1qa",
                      test = "polar_pvalue",
                      levels_order = c("Control_MV",  "Mid_treatment_MV", "Post_treatment_MV"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot5 <- boxplot_trio(res,
                      value = "Grn",
                      test = "polar_pvalue",
                      levels_order = c("Control_MV",  "Mid_treatment_MV", "Post_treatment_MV"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot6 <- boxplot_trio(res,
                      value = "Ctsb",
                      test = "polar_pvalue",
                      levels_order = c("Control_MV",  "Mid_treatment_MV", "Post_treatment_MV"),
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


## @knitr Volcano2D_MV_DEGPlots

DEGdf <- mutate(Volcano3D_Results, sig=ifelse(Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$Post_treatment_MV<=-0.5 | Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$Post_treatment_MV>=0.5, "Significant", "Not Significant")) #Will have different colors depending on significance
head(DEGdf)

volc = ggplot(DEGdf, aes(Post_treatment_MV, -log10(pvalue))) + #log2Foldchange vs FDR
  geom_point(aes(col=sig)) + # the points are colored by significance
  scale_color_manual(values=c("blue", "red")) + 
  ggtitle("Volcano Plot Post MV Treatment ") #e.g. 'Title description'
volc+geom_text_repel(DEGdf=head(DEGdf, 5), aes(label=gene)) #adding text for the top 5 genes


## @knitr Volcano2D_MV_GOPathways

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2


# we want the log2 fold change 
original_gene_list <- Volcano3D_Results$Post_treatment_MV

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
measles <- as.data.frame(res@df$scaled)
measles <- rownames_to_column(measles, var = "gene")
measles$significant <- ifelse(measles$pvalue < .05, "Significant", NA)
measles<-na.omit(measles)

## @knitr Print_TablesMV
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(measles, server = TRUE,
                                              class = "compact",
                                              filter="top",
                                              rownames = FALSE,
                                              colnames = c("gene", "Control_MV", "Mid_treatment_MV", "Post_treatment_MV", "x", "y", "r", "angle", "z", "pvalue", "col", "col", "significant"),
                                              extensions = c('Buttons'),
                                              options = list(
                                                pageLength = 15,
                                                dom = 'Bfrtip',
                                                buttons = c('excel','csv','pdf','copy')
                                              ))))


###############################################################################################################################################################
###################################### And now let's do the same for Toca virus ###############################################################################

## @knitr Volcano3D_MV2
#setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Bulk_RNA/05_Output/Mouse/SingleR")
Raw_reads <- read.csv("Mouse_Bulk_RNA_seq_Raw_Reads_Volcano3D_Toca.csv", header = TRUE, row.names = 1)
Meta_Data_Volcano3D_MV <- read.csv("Meta_Data_Toca_3DVolcano.csv", header = TRUE)
#head(Meta_Data_Volcano3D_MV)
rownames(Meta_Data_Volcano3D_MV) <- Meta_Data_Volcano3D_MV$Sample_ID
keep <- c("Sample_Number", "Viral_Status")
Meta_Data_Volcano3D_MV <- Meta_Data_Volcano3D_MV[,keep]
colnames(Meta_Data_Volcano3D_MV) <- c("Sample", "Virus")
Meta_Data_Volcano3D_MV$Virus <- factor(Meta_Data_Volcano3D_MV$Virus)
#head(Meta_Data_Volcano3D_MV)


# rename the tissue types
Virus_Infection <- function(x){
  x <- switch(as.character(x), "Control_Toca"="Control_Toca", "Mid_treatment_Toca"="Mid_treatment_Toca",
              "Post_treatment_Toca"="Post_treatment_Toca")
  return(x)
}

Meta_Data_Volcano3D_MV$Virus <- unlist(lapply(Meta_Data_Volcano3D_MV$Virus, Virus_Infection))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
Meta_Data_Volcano3D_MV$Virus<- factor(Meta_Data_Volcano3D_MV$Virus, levels=c("Control_Toca",  "Mid_treatment_Toca", "Post_treatment_Toca"))

# Create the DEseq2DataSet object
deseq2DataToca <- DESeqDataSetFromMatrix(countData=round(Raw_reads), colData=Meta_Data_Volcano3D_MV, design= ~ Sample + Virus)


# initial analysis run
dds_DE <- DESeq(deseq2DataToca)
# likelihood ratio test on 'Pathotype'
dds_LRT <- DESeq(deseq2DataToca, test = "LRT", reduced=~1, parallel = TRUE) 

# create 'volc3d' class object for plotting
res <- deseq_polar(dds_DE, dds_LRT, "Virus", pcutoff = 0.05,
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
                      value = "BC147527",
                      test = "polar_pvalue",
                      levels_order = c("Control_Toca",  "Mid_treatment_Toca", "Post_treatment_Toca"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot2 <- boxplot_trio(res,
                      value = "Gm50076",
                      test = "polar_pvalue",
                      levels_order = c("Control_Toca",  "Mid_treatment_Toca", "Post_treatment_Toca"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot3 <- boxplot_trio(res,
                      value = "Fxyd3",
                      test = "polar_pvalue",
                      levels_order = c("Control_Toca",  "Mid_treatment_Toca", "Post_treatment_Toca"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

ggarrange(plot1, plot2, plot3, ncol=3)

plot4 <- boxplot_trio(res,
                      value = "C1qa",
                      test = "polar_pvalue",
                      levels_order = c("Control_Toca",  "Mid_treatment_Toca", "Post_treatment_Toca"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot")

plot5 <- boxplot_trio(res,
                      value = "Grn",
                      test = "polar_pvalue",
                      levels_order = c("Control_Toca",  "Mid_treatment_Toca", "Post_treatment_Toca"),
                      box_colours = c("blue", "red", "green3"),
                      text_size = 12,
                      stat_colour = "black",
                      stat_size = 3,
                      step_increase = 0.05,
                      plot_method = "ggplot") 

plot6 <- boxplot_trio(res,
                      value = "Ctsb",
                      test = "polar_pvalue",
                      levels_order = c("Control_Toca",  "Mid_treatment_Toca", "Post_treatment_Toca"),
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
DEGdf <- mutate(Volcano3D_Results, sig=ifelse(Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$Post_treatment_Toca<=-0.5 | Volcano3D_Results$pvalue<=0.05 & Volcano3D_Results$Post_treatment_Toca>=0.5, "Significant", "Not Significant")) #Will have different colors depending on significance
#head(DEGdf)

volc = ggplot(DEGdf, aes(Post_treatment_Toca, -log10(pvalue))) + #log2Foldchange vs FDR
  geom_point(aes(col=sig)) + # the points are colored by significance
  scale_color_manual(values=c("blue", "red")) + 
  ggtitle("Volcano Plot Post Toca Treatment ") #e.g. 'Title description'
volc+geom_text_repel(DEGdf=head(DEGdf, 5), aes(label=gene)) #adding text for the top 5 genes


## @knitr Volcano2D_MV_GOPathways2
# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2


# we want the log2 fold change 
original_gene_list <- Volcano3D_Results$Post_treatment_Toca

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
Toca <- as.data.frame(res@df$scaled)
Toca <- rownames_to_column(Toca, var = "gene")
Toca$significant <- ifelse(Toca$pvalue < .05, "Significant", NA)
Toca<-na.omit(Toca)

## @knitr Print_TablesToca
# Print the datatable with all markers
print( htmltools::tagList(DT::renderDataTable(Toca, server = TRUE,
                                              class = "compact",
                                              filter="top",
                                              rownames = FALSE,
                                              colnames = c("gene", "Control_Toca", "Mid_treatment_Toca", "Post_treatment_Toca", "x", "y", "r", "angle", "z", "pvalue", "col", "col", "significant"),
                                              extensions = c('Buttons'),
                                              options = list(
                                                pageLength = 15,
                                                dom = 'Bfrtip',
                                                buttons = c('excel','csv','pdf','copy')
                                              ))))




## @knitr Interrogating_low_viral_uptake
Merged <- read.csv(file="Mouse_Bulk_RNA_seq_Normalized_Meta_Data.csv", header = TRUE)
df = subset(Merged, select = -c(X))
colnames(df)[9] <- "T_cells"
colnames(df)[10] <- "B_cells"

#######################################################################################################################################################
################################ Markers of Phagocytes ##############################################################################################
## @knitr Interrogating_phagocyte_markers_correlation

s1 <- ggscatter(df, x = "T_cells", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#F0E442",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = T_cells), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s2 <- ggscatter(df, x = "B_cells", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#0072B2",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = B_cells), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s3 <- ggscatter(df, x = "Macrophages", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#009E73",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Macrophages), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s4 <- ggscatter(df, x = "Microglia", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "navy",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Microglia), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s5 <- ggscatter(df, x = "NK", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#52854C", 
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = NK), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s6 <- ggscatter(df, x = "Monocytes", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "red",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Monocytes), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s7 <- ggscatter(df, x = "Malignant", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#00AFBB",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Malignant), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s8 <- ggscatter(df, x = "Neurons", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#E7B800",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Neurons), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s9 <- ggscatter(df, x = "Endothelial", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "green",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Endothelial), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s10 <- ggscatter(df, x = "Astrocytes", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                 combine = TRUE, ylab = "Gene Expression", color = "#D55E00",
                 add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Astrocytes), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s11 <- ggscatter(df, x = "Oligodendrocytes", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                 combine = TRUE, ylab = "Gene Expression", color = "#CC79A7",
                 add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Oligodendrocytes), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s12 <- ggscatter(df, x = "Epithelial", y = c("C1qa", "Grn", "Ctsb", "Cd93", "Ctse", "Cd68"), size = 1,
                 combine = TRUE, ylab = "Gene Expression", color = "#4E84C4",
                 add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Epithelial), method = "spearman")  + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, ncol = 4, nrow = 3) 


## @knitr Interrogating_phagocyte_markers_boxplots

# Perorm pairwise comparisons for C1qa gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p1 <- ggboxplot(df, x = "Viral_Status", y = "C1qa",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$C1qa), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 20)  +
  ggtitle("Phagocytosis Gene Expression C1qa") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

# Perorm pairwise comparisons for Grn gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p2 <- ggboxplot(df, x = "Viral_Status", y = "Grn",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Grn), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 20)  +
  ggtitle("Phagocytosis Gene Expression Grn") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )


# Perorm pairwise comparisons for Ctsb gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p3 <- ggboxplot(df, x = "Viral_Status", y = "Ctsb",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Ctsb), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 18)  +
  ggtitle("Phagocytosis Gene Expression Ctsb") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )


# Perorm pairwise comparisons for Cd93 gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p4 <- ggboxplot(df, x = "Viral_Status", y = "Cd93",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Cd93), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 18) +
  ggtitle("Phagocytosis Gene Expression Cd93") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 



# Perorm pairwise comparisons for Ctse gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p5 <- ggboxplot(df, x = "Viral_Status", y = "Ctse",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Ctse), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 20)   +
  ggtitle("Phagocytosis Gene Expression Ctse") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


# Perorm pairwise comparisons for Cd68 gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p6 <- ggboxplot(df, x = "Viral_Status", y = "Cd68",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Cd68), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 15)  +
  ggtitle("Phagocytosis Gene Expression Cd68") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3, common.legend = FALSE)



#######################################################################################################################################################
################################ Markers of Immunoglobulin and Apoptotic Genes ######################################################################################

## @knitr Interrogating_immunoglobulin_apoptotic_markers_correlation
s1 <- ggscatter(df, x = "T_cells",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#F0E442",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = T_cells), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s2 <- ggscatter(df, x = "B_cells",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#0072B2",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = B_cells), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s3 <- ggscatter(df, x = "Macrophages",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#009E73",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Macrophages), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s4 <- ggscatter(df, x = "Microglia",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "navy",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Microglia), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s5 <- ggscatter(df, x = "NK",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#52854C", 
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = NK), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s6 <- ggscatter(df, x = "Monocytes",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "red",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Monocytes), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s7 <- ggscatter(df, x = "Malignant", y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#00AFBB",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Malignant), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s8 <- ggscatter(df, x = "Neurons",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "#E7B800",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Neurons), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s9 <- ggscatter(df, x = "Endothelial",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                combine = TRUE, ylab = "Gene Expression", color = "green",
                add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Endothelial), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s10 <- ggscatter(df, x = "Astrocytes",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                 combine = TRUE, ylab = "Gene Expression", color = "#D55E00",
                 add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Astrocytes), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s11 <- ggscatter(df, x = "Oligodendrocytes", y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                 combine = TRUE, ylab = "Gene Expression", color = "#CC79A7",
                 add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Oligodendrocytes), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  ) 


s12 <- ggscatter(df, x = "Epithelial",  y = c("Igkc", "Iglc1", "Ighg2c", "Bcl2", "Casp4", "Casp7"), size = 1,
                 combine = TRUE, ylab = "Gene Expression", color = "#4E84C4",
                 add = "reg.line", conf.int = TRUE, ggtheme = theme_bw()) +
  stat_cor(aes(color = Epithelial), method = "spearman") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )  


ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, ncol = 4, nrow = 3) 


## @knitr Interrogating_immunoglobulin_apoptotic_markers_boxplots
# Perorm pairwise comparisons for Igkc gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p1 <- ggboxplot(df, x = "Viral_Status", y = "Igkc",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Igkc), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 25)  +
  ggtitle("Immunolglobulin Gene Expression Igkc") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

# Perorm pairwise comparisons for Iglc1 gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p2 <- ggboxplot(df, x = "Viral_Status", y = "Iglc1",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Iglc1), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 15)  +
  ggtitle("Immunolglobulin Gene Expression Iglc1") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )


# Perorm pairwise comparisons for Ighg2c gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p3 <- ggboxplot(df, x = "Viral_Status", y = "Ighg2c",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Ighg2c), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 20)  +
  ggtitle("Immunolglobulin Gene Expression Ighg2c") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )


# Perorm pairwise comparisons for Bcl2 gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p4 <- ggboxplot(df, x = "Viral_Status", y = "Bcl2",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Bcl2), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 15)  +
  ggtitle("Apoptosis Gene Expression Bcl2") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )


# Perorm pairwise comparisons for Casp4 gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p5 <- ggboxplot(df, x = "Viral_Status", y = "Casp4",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Casp4), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 20)  +
  ggtitle("Apoptosis Gene Expression Casp4") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )



# Perorm pairwise comparisons for Casp7 gene
df$Viral_Status <- factor(df$Viral_Status, levels = c("Control_MV", "HI_MV", "Mid_treatment_MV", "Post_treatment_MV", "Control_Toca", "HI_Toca", "Mid_treatment_Toca", "Post_treatment_Toca"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control_MV", "HI_MV"),  c("Control_MV", "Mid_treatment_MV"), c("Control_MV", "Post_treatment_MV"), c("Control_Toca", "HI_Toca"), c("Control_Toca", "Mid_treatment_Toca"), c("Control_Toca", "Post_treatment_Toca"))
p6 <- ggboxplot(df, x = "Viral_Status", y = "Casp7",  add = "jitter", legend = "none",
                color = "Viral_Status", palette = c("#00AFBB", "#E7B800",
                                                    "green", "red", "navy", "#009E73",
                                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                    "#52854C", "#4E84C4"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$Casp7), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 12)  +
  ggtitle("Apoptosis Gene Expression Casp7") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3, common.legend = FALSE)

## @knitr KEGG_Pathway_Phagocytosis

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# reading in data from deseq2
df = read.csv("deseq2_control_vs_MV_PT_Pathways.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

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
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "mouse"

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
phagocytosis <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04666", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
phagocytosis <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04666", species = kegg_organism, kegg.native = F)

## @knitr KEGG_Pathway_Apoptosis
apoptosis <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04210", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
apoptosis <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04210", species = kegg_organism, kegg.native = F)
