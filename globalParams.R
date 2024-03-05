###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#
# Parameters for individual samples (e.g. SAMPLE_NAME, PATH_DATA) must be set
# in the file 'sampleParams.R' (stored in the sample folder). Any value defined 
# there will override values defined in this file.
#

#### General

PROJECT_NAME = "UCSF Aaron Diaz Lab - Mouse Bulk RNA";
PROJECT_SHORT_NAME = "Viral_Infection"

PATH_PROJECT = "/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/09_Bulk_RNA/01_Input_Data/Mouse/"

#### Input / Output

# RawData folder name of the experiment
#PATH_RAWDATA = file.path("/media/owen/Backup Plus/UCSF_Project/01_RAW_DATA/PNOC005/cellranger/Mouse_Virus_Experiments/")

# Processed data for all experiments
PATH_PROCESSED_DATA = file.path("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/09_Bulk_RNA/01_Input_Data/Mouse/")

#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



