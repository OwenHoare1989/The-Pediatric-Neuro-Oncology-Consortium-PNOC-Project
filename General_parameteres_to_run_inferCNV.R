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
Threshold = 0.1  # 0.1 for 10x-genomics
my_mode='subclusters'
Cluster_settings = TRUE
CNV_cluster_to_hide = 8
Smooth = 'pyramidinal' # based on genes-center base pair distance
Denoise = TRUE
Hmm_settings = TRUE
Hmm_mode ='i6'
Hmm_p_val = 0.3
Subclusters ="leiden"
Subcluster_p_val = 0.05
my_mode=FALSE
Plotting_steps=FALSE
amplifier=3  # sets midpoint for logistic
Logistic=TRUE # turns gradient filtering on
Number_of_threads = 12
diagnostics = TRUE
HMM_report_by = "subcluster"
min_cells_per_gene = 3


## @knitr Run_inferCNV_analysis
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             analysis_mode = my_mode,
                             cutoff=Threshold,
                             cluster_by_groups=Cluster_settings,
                             k_obs_groups=CNV_cluster_to_hide,
                             out_dir=out_dir,
                             smooth_method = Smooth,
                             denoise=Denoise,
                             plot_steps=Plotting_steps,
                             sd_amplifier=amplifier,
                             noise_logistic=Logistic,
                             HMM=Hmm_settings,
                             HMM_type = Hmm_mode,
                             HMM_report_by = "subcluster",
                             min_cells_per_gene = 3,
                             tumor_subcluster_partition_method=Subclusters,
                             tumor_subcluster_pval = Subcluster_p_val,
                             BayesMaxPNormal = Hmm_p_val,
                             resume_mode=my_mode,
                             num_threads=Number_of_threads
)


