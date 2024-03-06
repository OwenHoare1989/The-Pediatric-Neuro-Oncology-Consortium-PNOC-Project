library(speckle)
setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/08_TME_Analysis/05_Output/Immune_Cell_Heterogeneity_SingleR/")
Seurat.viral.analysis_SingeR_true <- readRDS("Seurat.viral.analysis_SingeR_true.rds")
Seurat.viral.analysis_SingeR_true = subset( Seurat.viral.analysis_SingeR_true, subset = (SingleR_true_annot == "none"), invert =TRUE)
saveRDS(Seurat.viral.analysis_SingeR_true, file="Seurat.viral.analysis_SingeR_true.rds")

D1 <- dittoSeq::dittoDimPlot(Seurat.viral.analysis_SingeR_true, reduction.use = "umap", "SingleR_true_annot", color.panel = c("hotpink", "#E7B800", "#FC4E07", "skyblue","navy","green","#CC313D","#0F8554", "#CC3A8E", "#00AFBB","cyan1","brown", "blue"),  size = 0.7)
D2 <- dittoSeq::dittoBarPlot(Seurat.viral.analysis_SingeR_true, var = "Virus_Infection_true", split.by = "Treatment", group.by = "Sample_ID",  color.panel = c("hotpink", "#E7B800", "#FC4E07", "skyblue","navy","green","#CC313D","#0F8554", "#CC3A8E", "#00AFBB","cyan1","brown", "blue"))

D1 + D2


dittoSeq::dittoDimPlot(Seurat.viral.analysis_SingeR_true, reduction.use = "umap", "Virus_Infection_true", split.by = "Virus_Infection_true",   color.panel = c( "orange", "dodgerblue"),  size = 0.7)
Seurat.viral.analysis_SingeR_true[["Virus_Infection_true"]] <- paste0(Seurat.viral.analysis_SingeR_true$Sample_ID, "_", Seurat.viral.analysis_SingeR_true$Virus_Infection)
Seurat.viral.analysis_SingeR_true@meta.data$Virus_Infection_true <- plyr::mapvalues(Seurat.viral.analysis_SingeR_true@meta.data$Virus_Infection_true, 
                                                  from = c("MV_302_Non_Infected", "MV_302_Measles_virus_infected",       "MV_310_Non_Infected",   "MV_310_Measles_virus_infected",      
                                                           "MV_311_Non_Infected",   "MV_311_Measles_virus_infected",       "MV_312_Non_Infected", "MV_312_Measles_virus_infected",      
                                                           "MV_313_Non_Infected",  "MV_313_Measles_virus_infected",       "MV_PD1_601_Measles_virus_infected",   "MV_PD1_601_Non_Infected",            
                                                           "MV_PD1_602_1_Non_Infected",           "MV_PD1_602_1_Measles_virus_infected", "MV_PD1_602_2_Measles_virus_infected", "MV_PD1_602_2_Non_Infected",         
                                                           "MV_PD1_604_Non_Infected",             "MV_PD1_604_Measles_virus_infected"  ),
                                                  to = c("Non_Infected", "Measles_virus_infected",       "Non_Infected",   "Measles_virus_infected",      
                                                         "Non_Infected",   "Measles_virus_infected",       "Non_Infected", "Measles_virus_infected",      
                                                         "Non_Infected",  "Measles_virus_infected",       "Measles_virus_infected",   "Non_Infected",            
                                                         "Non_Infected",           "Measles_virus_infected", "Measles_virus_infected", "Non_Infected",         
                                                         "Non_Infected",             "Non_Infected"  ))

D3 <- dittoSeq::dittoDimPlot(Seurat.viral.analysis_SingeR_true, "Virus_Infection", split.by = "SingleR_true_annot",   size = 0.7)
D4 <- dittoSeq::dittoBarPlot(Seurat.viral.analysis_SingeR_true, var = "Virus_Infection", split.by = "SingleR_true_annot", group.by = "Sample_ID")

D3 + D4


D5 <- dittoSeq::dittoDimPlot(Seurat.viral.analysis_SingeR_true, "reassigned.malignancy", split.by = "Virus_Infection",   size = 0.7)
D6 <- dittoSeq::dittoBarPlot(Seurat.viral.analysis_SingeR_true, var = "reassigned.malignancy", split.by = "Virus_Infection", group.by = "Sample_ID")

D5 + D6

setwd("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/Figures_PNOC_Paper/")

# Run propeller testing for cell type proportion differences between the two groups
propeller_test <- propeller(clusters = Seurat.viral.analysis_SingeR_true$SingleR_true_annot, sample = Seurat.viral.analysis_SingeR_true$Sample_ID,
                            group = Seurat.viral.analysis_SingeR_true$Treatment)

write.csv(propeller_test, file= "propeller_test.csv")
