## @knitr test
Sys.setenv(RETICULATE_MINICONDA_ENABLED=FALSE)
######## R Libraries
library(infercnv)
library(knitr)
library(dplyr)
library(future)
plan("multicore", workers = 12)
options(future.globals.maxSize = 8000 * 1024^2)
# library(Seurat)
# setwd("/mnt/DOSI/BNSRLAB/BIOINFO/Project/BIOFLIRT/03_BIOINFO_ANALYSIS/04_BCR_Diversity/05_Output/CimlSamples_BCR_Diversity_Analysis_Merged/05_B_cell_malignancy_assessment/GG67")
# GG67 <- readRDS(file="GG67_seurat_malignancy_clones_attributed.rds")
# setwd("/mnt/DOSI/BNSRLAB/BIOINFO/Project/BIOFLIRT/03_BIOINFO_ANALYSIS/03_Global_Heterogeneity_Analysis/05_Output/210212_BIOFLIRT_EXP6BL")
# seurat <- readRDS(file="sc10x.rna.seurat.RDS")
# GG67.T.cells <- SplitObject(seurat, split.by = "RNA_snn_res.1.B.TME.Assignment")

# 
# TME.cells = GG67.T.cells$TME_cells
# T_cells = subset(TME.cells, T_signature > 1)
# GG67_B.T.cells = merge(x = GG67, y = T_cells)
# 
# GG67_B.T.cells.meta.data = GG67_B.T.cells@meta.data
# setwd("/mnt/DOSI/BNSRLAB/BIOINFO/Project/BIOFLIRT/03_BIOINFO_ANALYSIS/05_InferCNV/GG67_inferCNV")
# write.csv(GG67_B.T.cells.meta.data, file="GG67_B.T.cells.meta.data.csv")
# 
# 
# expr_raw <- GetAssayData(object = GG67_B.T.cells, assay = "RNA", slot = "counts")
# expr <- as(Class = 'matrix', object = expr_raw)
# 
# write.csv(expr, file="GG67_inferCNV_malignancy_assessment_raw.counts.csv", row.names = TRUE)
# 
# 

###########################
#### Tidy position data ###
###########################
library(tidyverse)
position_df <- read.table(file="Gene_Chromosome_position_file.txt")
# 
# # Add colnames to position dataframe
colnames(position_df) <- c("gene_name","chromosome","start","end")
# 
# # Relevel chromosome data as from raw 23 to 38 it correspond to clones
position_df$chromosome <- fct_relevel(position_df$chromosome, c(as.character(1:22),"X","Y",levels(position_df$chromosome)[23:38]))
# 
# # Order genes first according to ascending position in chromosome then according to ascending start position
position_df <- position_df %>% dplyr::arrange(chromosome, start)

# Remove Colnames
colnames(position_df) <- NULL
# 
# # Save the data
Reference_meta.data <- read.csv(file="Reference_meta.data.csv")

write.table(Reference_meta.data, file = "Reference_meta.data.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

## @knitr loadData
setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/03_inferCNV/01_Input_Data/")
out_dir = "/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/03_inferCNV/05_Output/"

infercnv_obj <- c()

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="raw_matrix.txt",
                                    annotations_file="Reference_meta.data.txt",
                                    delim="\t",
                                    gene_order_file="your_gen_pos.txt",
                                    ref_group_names=c("Immune_cells (non-malignant)"))


# Run inferCNV parameters
CUTOFF = 0.1  # 0.1 for 10x-genomics
ANALYSIS_MODE='subclusters'
CLUSTER_BY_GROUP_PARAM = TRUE
NUMBER_OF_HIDDEN_CNV_CLUSTERS = 8
SMOOTH_METHOD = 'pyramidinal' # based on genes-center base pair distance
DENOISE_PARAM = TRUE
HMM_PARAM = TRUE
HMM_TYPE ='i6'
HMM_PVAL = 0.3
SUBCLUSTER_METHOD ="leiden"
TUMOR_SUBCLUSTER_PVAL = 0.05
MODE=FALSE
PLOT_STEPS=FALSE
SD_AMPLIFIER=3  # sets midpoint for logistic
NOISE_LOGISTIC=TRUE # turns gradient filtering on
NBR_THREAD = 12
diagnostics = TRUE
HMM_report_by = "subcluster"
min_cells_per_gene = 3


## @knitr Run_inferCNV_analysis
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             analysis_mode = ANALYSIS_MODE,
                             cutoff=CUTOFF,
                             cluster_by_groups=CLUSTER_BY_GROUP_PARAM,
                             k_obs_groups=NUMBER_OF_HIDDEN_CNV_CLUSTERS,
                             out_dir=out_dir,
                             smooth_method = SMOOTH_METHOD,
                             denoise=DENOISE_PARAM,
                             plot_steps=PLOT_STEPS,
                             sd_amplifier=SD_AMPLIFIER,
                             noise_logistic=NOISE_LOGISTIC,
                             HMM=HMM_PARAM,
                             HMM_type = HMM_TYPE,
                             HMM_report_by = "subcluster",
                             min_cells_per_gene = 3,
                             tumor_subcluster_partition_method=SUBCLUSTER_METHOD,
                             tumor_subcluster_pval = TUMOR_SUBCLUSTER_PVAL,
                             BayesMaxPNormal = HMM_PVAL,
                             resume_mode=MODE,
                             num_threads=NBR_THREAD
)


