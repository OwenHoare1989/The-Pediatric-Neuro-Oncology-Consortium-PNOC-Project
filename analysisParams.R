###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "02_Global_Heterogeneity_Analysis"
LITERAL_TITLE = "Cell heterogeneity exploration"

# Maximum number of variable features to keep for PCA analysis
SELECTION_METHOD     = "vst"
NVARIABLE_FEATURES   = 2000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE
DATA_SCALE        = FALSE

# PCA parameters
PCA_NPC              = 50  # Default number of dimensions to use for PCA (see Seurat::RunPCA())

# UMAP/TSNE parameters (see Seurat::RunUMAP() and Seurat::RunTSNE())
REDUCTION = "pca"
UMAP_METHOD = "uwot"
METRIC = "correlation"
MAX_VIZDIMS = 20

# Clustering parameters (see Seurat::FindNeighbors() and Seurat::FindClusters())
MAX_CLUSDIMS = 20
FORCE = TRUE

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.75;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10;       # Number of marker genes to show in report and tables (NULL for all)


# Signatures used for scoring (see Seurat::AddModuleScore())
SIGNATURES<-list("Oligodendrocyte"=c("Syt1", "Mobp","Fyn", "Ptgds", "Olig1", "Cldn11", "Apod", "Mal"),
                 "Myeloid_signature" = c("Lyz1", "Lgals2","Cd14"),
                 "Astrocyte_signature"= c("Agt","Aqp4", "Gfap", "Gja1"),
                 "Endothelial_signature" = c("Cldn5", "Adgrf5"),
                 "Fibroblast_signature" = c("Col1a1", "Col1a2", "Col5a1", "Loxl1", "Lum", "Fbln1", "Fbln2", "Cd34", "Pdgfra", "Fgf2", "Fgfr2"),
                 "Macrophage_microglia_signature" = c("Jchain", "Cx3cr1", "C1qa", "Csf1r", "Tmem119", "P2ry12"),
                 "Neuron_signature"  =c("Gad1", "Slc17a6",  "Snap25", "Syt1"),
                 "OPC_signature"  =c("Pdgfra", "Cspg4",  "Olig1"),
                 "VLMC_signature"  =c("Igfbp2", "Tbx18",  "Slc6a13", "Pdgfra", "Apod", "Ptgds", "Gja1"),
                 "B_signature"= c("Ms4a1", "Cd79a","Cd79b","Cd19"),
                 "T_signature"= c("Cd3e","Cd3d"),
                 "NK_signature"=c("Klrd1", "Klrc1",  "Xcl1", "Cx3cr1", "Ncr1"))
                 
                
  




























