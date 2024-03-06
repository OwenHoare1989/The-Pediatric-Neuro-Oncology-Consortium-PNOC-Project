setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/07_Investigating_Low_Viral_Uptake/")
library(Seurat)
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
library(dittoSeq)
library(SCpubr)
library("Nebulosa")
library(MetaVolcanoR)
plan("multicore", workers = 12)
options(future.globals.maxSize = (1014*1024^2)*60)
seurat.data <- readRDS("sc10x.rna.seurat_sorted_cellTypes.rds")
seurat.data <- subset(seurat.data, subset = CellType_Phenotypes != "NA")

seurat.data <- readRDS("Three_way_neoplastic_viral_cell_analysis.rds")
# seurat.data <- subset(seurat.data, subset = Sample_ID != "MV_PD1_604")
seurat.data <- subset(seurat.data, subset=CellType_Phenotypes %in% c("Microglia", "Monocytes",         "M1_Macrophages",            
                                                                      "M2_Macrophages",    "M0_Macrophages" ))

# # Use the RNA assay to find variable genes
DefaultAssay(seurat.data)  = "RNA"
seurat.data = FindVariableFeatures(seurat.data, selection.method = "vst", nfeatures = 2000)
var_genes_from_rna = seurat.data@assays$RNA@var.features

# # Use the anchors features
# DefaultAssay(seurat.data)  = "integrated"
# features = seurat.data@assays$integrated@var.features

# # # Use the integrated assay (with all genes integrated) to find variable genes 
# DefaultAssay(seurat.data)  = "integrated"
# seurat.data = FindVariableFeatures(seurat.data, selection.method = "vst", nfeatures = 2000)
# var_genes_from_int = seurat.data@assays$integrated@var.features

# Plot Variable genes
# Identify the 20 most highly variable genes
top20 <- head(var_genes_from_rna, 50)
# 
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.data)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
print(plot2)



SIGNATURES<-list("infection_genes"= c("Slamf1", "Slamf6", "Slamf7", "Slamf8", "Slamf9", "Cd46"))

for(i in seq_along(SIGNATURES)){
  seurat.data<- AddModuleScorePerso(seurat.data, assay = "RNA", features = list(SIGNATURES[[i]]), name = names(SIGNATURES)[i])
}



dittoSeq::dittoDimPlot(seurat.data, "CellType_Phenotypes",  reduction.use = "umap",  color.panel = c("#0062B4FF", "#FFBF47", "#FF2700FF", "mediumseagreen",  "sienna4", "violet", "violetred", "violetred4", "#CC313D", "blue", "limegreen", "#3B9AB2","deeppink", "#BF812D"),  size = 0.7)
D2 <- dittoSeq::dittoBarPlot(seurat.data, var = "CellType_Phenotypes", group.by = "Sample_ID", split.by = "Treatment", color.panel = c("#0062B4FF", "#FFBF47", "#FF2700FF", "mediumseagreen",  "sienna4", "violet", "violetred", "violetred4", "#CC313D", "blue", "limegreen", "#3B9AB2","deeppink", "#BF812D"))


Phagocytic_markers <- c("Ctsb", "Grn", "C1qa") 

dittoBoxPlot(sample, "C1qa", group.by = "Treatment")
dittoBoxPlot(seurat.data, "Slamf8",  group.by = "Virus_Infection_true", split.by = "CellType_Phenotypes")
dittoBoxPlot(seurat.data, "C1qa",  group.by = "Virus_Infection_true", split.by = "CellType_Phenotypes")


sample <- seurat.data
# Use custom colors.
p <- SCpubr::do_AlluvialPlot(sample = sample,
                             first_group = "Virus_Infection_true",
                             middle_groups = "CellType_Phenotypes",
                             last_group = "Treatment", fill.by = "CellType_Phenotypes",
                              repel = TRUE,  use_labels = TRUE,
                             stratum.fill.conditional = FALSE,
                             use_geom_flow = FALSE,
                             alluvium.color = "white",
                             flow.color = "white",
                             flip = FALSE,
                             label.color = "black",
                             curve_type = "sigmoid",
                             viridis_color_map = "G",
                             viridis_direction = -1,
                             plot.grid = TRUE,
                             grid.color = "grey75",
                             grid.type = "dashed",
                             na.value = "white",
                             legend.position = "right",
                             legend.title = "Cell Types")

p

p <- SCpubr::do_ChordDiagramPlot(sample = sample,
                                 from = "CellType_Phenotypes",
                                 to = "Virus_Infection_true",
                                 link.arr.type = "triangle", repel = TRUE)

p




p1 <- SCpubr::do_FeaturePlot(sample = sample, 
                             features = c("C1qa", "C1qb", "C1qc", "C2", "Fcgr2b", "Fcgr3")) 

p2 <- SCpubr::do_NebulosaPlot(sample = sample, 
                              features = c("Ctsb", "Grn", "C1qa", "Fcgr2b", "Ighg2c", "Igkc"))
p <- p1 | p2
p


p <- SCpubr::do_BoxPlot(sample = sample,
                                  feature = c(),
                                  use_test = TRUE,
                                  group.by = "CellType_Phenotypes",
                                  comparisons = list(c("Malignant", "Microglia"),
                                                     c("Malignant", "M0_Macrophages")),
                                                     map_signif_level = FALSE)
p


# Retrieve all dimensionality reductions
getReductions(seurat.data)
phagPlot1 <- dittoDimPlot(seurat.data, "Ctsb", split.by = "Virus_Infection_true")
phagPlot2 <- dittoDimPlot(seurat.data, "Grn", split.by = "Virus_Infection_true")
phagPlot3 <- dittoDimPlot(seurat.data, "C1qa", split.by = "Virus_Infection_true")

phagPlot1 + phagPlot2 + phagPlot3


# Annotating and ordering cells by some meaningful feature(s):
dittoHeatmap(seurat.data, Phagocytic_markers,
             scaled.to.max = FALSE,
             show_colnames = FALSE,
             show_rownames = TRUE,
             complex = TRUE,
             annot.by = c("Virus_Infection_true", "CellType_Phenotypes"),
             order.by = "CellType_Phenotypes")

dittoDotPlot(seurat.data, vars = Phagocytic_markers, group.by = "CellType_Phenotypes", split.by = "Virus_Infection_true")

dittoPlotVarsAcrossGroups(seurat.data, Phagocytic_markers, group.by = "Virus_Infection_true", split.by = "CellType_Phenotypes",
                          main = "Markers of Mouse Phagocytosis")

dittoDimPlot(seurat.data, infection_genes,
             split.by = c("Virus_Infection_true"))

dittoRidgePlot(seurat.data, "C1qa", group.by = "CellType_Phenotypes")
dittoRidgePlot(seurat.data, Phagocytic_markers, group.by = "CellType_Phenotypes", split.by = "Virus_Infection_true",
               split.ncol = 1)

# Increase cell size.
p <- SCpubr::do_CorrelationPlot(sample = seurat.data, column_title = "Cell Type Annotation",
                                group.by = "Treatment",
                                column_names_rot = 45,
                                legend.position = "bottom",
                                cell_size = 12)
p

p1 <- SCpubr::do_CorrelationPlot(sample = sample, column_title = "MV Infection",
                                group.by = "Virus_Infection",
                                column_names_rot = 45,
                                legend.position = "bottom",
                                cell_size = 60)
p1


genes <- c("Cd274", "Pdcd1", "Arg1", "Ighg2b", "Iglc2", "Igkc", "Ighg2c", "Cd4", "Cd8a", "Cd3e",
           "Jchain", "Igha", "Cxcl9", "Cxcl10", "Ccl2", "Ccl7", "Ccl8", "Ccr1", "Ccr5",
           "C1qa", "C1qb", "C1qc", "C2", "Fcgr2b", "Fcgr3", "H2-Aa", "H2-Q4", "H2-Ab1", "H2-Eb1",
           "Lag3", "Cxcl3", "Ctla4", "Klrd1", "Klrc1", "Xcl1", "Cx3cr1", "Ncr1", "Lyz1", "Cd19", "Cd14",
           "Csf1", "Il1b", "Fos", "Jun", "Junb", "Rela", "Relb", "Ptgs2", "Stat3")

# Default parameters.
p <- SCpubr::do_ExpressionHeatmap(sample = seurat.data,
                                features = genes, flip = FALSE,
                                group.by = c("CellType_Phenotypes", "Virus_Infection_true", "Treatment"),
                                viridis_direction = -1,
                                viridis_color_map = "C")
p



# Seurat sample.
sample <- seurat.data

# Set the identities correctly.
Seurat::Idents(sample) <- sample$Treatment

# Compute DE genes and transform to a tibble.
de_genes <- tibble::tibble(Seurat::FindAllMarkers(object = sample,  logfc.threshold = 0.1, return.thresh = 0.01, only.pos = TRUE))
de_genes <- as.data.frame(de_genes)
de_genes$significant <- ifelse(de_genes$p_val_adj < .05, "Significant", NA)
de_genes<-na.omit(de_genes)
# Add more layers of mean expression with group.by.
p <- SCpubr::do_GroupwiseDEPlot(sample = sample,
                                de_genes = de_genes,
                                top_genes = 20,
                                cell_size = 6,
                                column_title = "Differential expression T cells control vs treated",
                                row_title_p_values = "P-Values",
                                row_title_logfc = "FC",
                                viridis_map_expression = "C",
                                 group.by = c("CellType_Phenotypes",
                                             "Virus_Infection", 
                                             "Treatment"),
                                row_title_expression = c("Cell Type",
                                                         "MV",
                                                         "Treatment"))
p

my_genes <- FindMarkers(sample, ident.1 = "Treated", ident.2 = "Control", logfc.threshold = 0.1, test.use = "MAST", group.by = "Treatment")
my_genes <- rownames_to_column(my_genes, var = "gene")
EnhancedVolcano(my_genes,
                lab = as.character(my_genes$gene),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "DEG Control vs Treated Innate Immune Cells",
                pCutoff = 0.05,
                FCcutoff = 0,
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

ggplot(data=my_genes, aes(x=gene, y=avg_log2FC, fill = p_val_adj)) +
  geom_bar(stat="identity") + 
  coord_flip()

head(meta_degs_vote@metaresult, 3)
meta_degs_vote@degfreq
write.csv(de_genes, file="My_de_genes.csv")
de_genes2 <- read.csv("de_genes2.csv")
library(decoupleR)
library(OmnipathR)
library(clusterProfiler)
library(org.Mm.eg.db)

genes.use <- c("Ccl4",     "Cx3cr1" ,  "C1qc" ,    "Lgmn" ,    "Jun"  ,    "Atf3"  ,   "Fos" ,     "Ccl3"  ,   "Trem2"  ,  "C1qb"  ,   "Ly86"  ,   "Junb"  ,   "Ctss" ,   
           "C1qa",     "Tyrobp",   "Mafb"  ,   "Ctsl" ,    "Zfp36" ,   "Fcer1g"  , "Csf1r" ,   "Cd83" ,    "Sat1"  ,   "Ctsz" ,    "Dusp1" ,   "Plek"  ,   "Ctsd"  ,  
           "Btg2",     "Egr1" ,    "Nfkbia" ,  "Cd9"  ,    "Ccrl2" ,   "Selplg" ,  "Cebpb" ,   "P2ry12"  , "Rhob"  ,   "Gpr34"  ,  "Rgs10"  ,  "Ctsh" ,    "Fcgr3"  , 
           "Ier5" ,    "Laptm5" ,  "Cyba"  ,   "Tmem119"  ,"Tgfbr1" ,  "Unc93b1" , "Rgs1" ,    "Grn"   ,   "Lyn"   ,   "Cd300c2" , "Cd81"   ,  "Ctsc"  ,   "Rnase4" , 
           "Btg1",     "Klf6",     "Serinc3" , "Ppp1r15a" ,"Ier2"   ,  "Il1a")

# Compute the grouped GO terms.
# out <- SCpubr::do_GroupedGOTermPlot(genes = genes.use,
#                                     org.db = org.Mm.eg.db)
# 
# # Plot the output combined heatmap.
# out$Plots$BP$Combined
# 
# # Compute the grouped GO terms.
# out <- SCpubr::do_GroupedGOTermPlot(genes = genes.use,
#                                     org.db = org.Mm.eg.db,
#                                     cluster_rows = FALSE,
#                                     cluster_cols = FALSE)
# 
# # Plot the output combined heatmap.
# out$Plots$BP$Combined
# 
# 
# # Compute the grouped GO terms.
# out <- SCpubr::do_GroupedGOTermPlot(genes = genes.use,
#                                     org.db = org.Mm.eg.db,
#                                     min.overlap = 2,
#                                     flip = FALSE,
#                                     levels.use = c(1, 2, 3, 4),
#                                     colors.use = c("lightblue", "navyblue"))
# 
# # Plot the output combined heatmap.
# out$Plots$BP$Combined



# Compute the grouped GO terms.
out <- SCpubr::do_FunctionalAnnotationPlot(genes = de_genes$gene,
                                           org.db = org.Mm.eg.db,
                                           min.overlap = 10,
                                           p.adjust.cutoff = 0.00001,
                                           pAdjustMethod = "BH",
                                           minGSSize = 10,
                                           maxGSSize = 500)

# Retrieve the heatmap.
out$Heatmap

# Retrieve the Bar and Dot plot.
out$BarPlot | out$DotPlot

# Retrieve the Tree plot.
out$TreePlot

library(enrichR)
# Set necessary enrichR global options. This is copied from EnrichR code to avoid having to load the package.
suppressMessages({
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr"))
  
  # Set the search to Human genes.
  enrichR::setEnrichrSite(site = "Enrichr")
  
  websiteLive <- TRUE
  dbs <- enrichR::listEnrichrDbs()
  # Get all the possible databases to query.
  dbs <- sort(dbs$libraryName)
})

# Choose the dataset to query against.
dbs_use <- c("GO_Biological_Process_2021", 
             "GO_Cellular_Component_2021", 
             "Azimuth_Cell_Types_2021")

# List of genes to use as input.
genes <- c("Ctsb", "Grn", "C1qa")

# Retrieve the enriched terms.
enriched_terms <- enrichR::enrichr(genes.use, dbs_use)

# Default plot.
p <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
p

# Modify font size of the terms.
p1 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms)
p2 <- SCpubr::do_TermEnrichmentPlot(enriched_terms = enriched_terms,
                                    text_labels_size = 6)

p1
p2


library(SeuratWrappers)

# Define your sample and assay.
sample <- seurat.data
assay <- sample@assays$RNA@scale.data



# Retrieve prior knowledge network.
network <- decoupleR::get_progeny(organism = "mouse", top=500)

# Run weighted means algorithm.
activities <- decoupleR::run_wmean(mat = as.matrix(assay),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)

# General heatmap.
library(ggdist)
out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                      group.by = "CellType_Phenotypes",
                                      activities = activities,
                                      split.by = "Virus_Infection_true",
                                      plot_FeaturePlots = TRUE,
                                      plot_GeyserPlots = TRUE)
p <- out$heatmaps$average_scores
p


# Retrieve feature plots.
out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                      activities = activities,
                                      group.by = "CellType_Phenotypes",
                                      plot_FeaturePlots = TRUE)
p1 <- SCpubr::do_DimPlot(sample)
p2 <- out$feature_plots$EGFR
p <- p1 | p2
p

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)
my_genes <- as.data.frame(my_genes)
my_genes$significant <- ifelse(my_genes$p_val_adj < .05, "Significant", NA)
my_genes<-na.omit(my_genes)

# we want the log2 fold change 
original_gene_list <- my_genes$avg_log2FC

# name the vector
names(original_gene_list) <- my_genes$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign)

gseaplot2(gse, geneSetID = 1639, pvalue_table = TRUE,  ES_geom = "line", title = gse$Description[1639])

# Seurat sample.
sample <- seurat.data

# Set the identities correctly.
Seurat::Idents(sample) <- sample$CellType_Phenotypes

pbmc.markers <- FindAllMarkers(sample, min.pct = 0.1, logfc.threshold = 1, test.use = "MAST", group.by = "CellType_Phenotypes")
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

##################################################################
#Subsetting top 100 markers with adjusted p values lower than .05#
##################################################################
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
dotplot(GOclusterplot)


ggplot(top100pval, aes(x = gene, y = avg_log2FC, fill = factor(cluster))) +
  geom_col(alpha = 0.8) +
  scale_fill_brewer(palette = "Dark2", name = "Differential Expression Clusters") +
  theme_bw() +
  theme(
    legend.position = c(0.8, 0.2),
    legend.background = element_blank(),
    legend.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA)
  ) +
  labs(x = "Subject", y = "Change from baseline",
       title = "Week 4 outcome by treatment")


# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "SYMBOL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = my_genes[my_genes$gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$avg_log2FC

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

dotplot(kk2, showCategory = 15, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)




my_genes1 <- my_genes[which(abs(my_genes$avg_log2FC) > 0.5 & my_genes$p_val_adj < 0.05),]

genes <- my_genes1$gene


mat<- sample[["RNA"]]@data[genes, ] %>% as.matrix()

## scale the rows
mat<- t(scale(t(mat)))

cluster_anno<- sample@meta.data$Treatment



# what's the value range in the matrix
quantile(mat, c(0.1, 0.95))

Seurat::PurpleAndYellow()


## make the black color map to 0. the yellow map to highest and the purle map to the lowest
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
# ha = HeatmapAnnotation(type = sample@meta.data$CellType_Phenotypes,
#                        col = list(type = c("Astrocytes" =  "#0062B4FF", "B cells" = "#FFBF47", "CD4+_T_cells" = "#FF2700FF", "Endothelial cells" =  "mediumseagreen", "Epithelial cells" = "sienna4",
#                                         "M0_Macrophages" = "violet", "M1_Macrophages"  = "violetred", "M2_Macrophages"  = "violetred4", "Malignant" = "#CC313D", 
#                                            "Microglia" = "blue", "Monocytes" = "limegreen", "Neurons" = "#3B9AB2", "NK_cells" = "deeppink", "Oligodendrocytes" = "#BF812D")))

ha = HeatmapAnnotation(type = sample@meta.data$Treatment,
                       col = list(type = c("Treated" =  "#0062B4FF", "Control" = "#FFBF47")))


MyCoolHeatmap1 <- Heatmap(mat, name = "Expression",  
                          column_split = factor(cluster_anno),
                          top_annotation = ha,
                          cluster_columns = TRUE,
                          show_column_dend = FALSE,
                          cluster_column_slices = TRUE,
                          column_title_gp = gpar(fontsize = 20),
                          column_gap = unit(0.5, "mm"),
                          cluster_rows = TRUE,
                          show_row_dend = FALSE,
                          col = col_fun,
                          row_names_gp = gpar(fontsize = 4),
                          clustering_distance_rows = "euclidean",
                          clustering_method_rows = "ward.D2",
                          column_title_rot = 0,
                          #top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                          show_column_names = FALSE,
                          show_heatmap_legend = TRUE,
                          use_raster = TRUE,
                          raster_resize_mat = TRUE,
                          raster_quality = 4)


print(MyCoolHeatmap1)


cellType = c("Malignant" = "#CC313D", "Neurons" = "#3B9AB2", "Endothelial cells" =  "mediumseagreen", "Monocytes" = "limegreen", "Microglia" = "blue", "M2_Macrophages" = "violetred4", "M1_Macrophages" = "violetred",  "M0_Macrophages" = "violet",  "CD4+_T_cells" = "#FF2700FF", "B cells" = "#FFBF47",  "Astrocytes" = "#0062B4FF", "Oligodendrocytes" = "#BF812D", "NK_cells" = "deeppink", "Epithelial cells" = "sienna4")
VlnPlot(seurat.data, features = c("Igkc", "H2-Eb1", "Slamf7", "Cd274", "Fcgr2b", "Cd46"), group.by = "CellType_Phenotypes", cols = cellType)
p = VlnPlot(seurat.data, features = c("Igkc", "H2-Eb1", "Slamf7", "Cd274", "Fcgr2b", "Cd46"), group.by = "Treatment", cols = c("#0062B4FF", "#FFBF47"))
p + theme(legend.position = 'right')

p1 = VlnPlot(seurat.data, features = c("Cd274", "Fcer1g"), group.by = "CellType_Phenotypes")
p1 + theme(legend.position = 'right')
