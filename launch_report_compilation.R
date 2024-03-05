# ####################################################################
# This script launch the compilation of both reports (one per sample)
# and rename them accordingly with the sample name
# ####################################################################

library( knitr)
library( rmarkdown)

### Define working folder (contains R/Rmd file for current sample, parent contains global project files)
WORKING_DIR   = getwd();

EXPERIMENT_NAME_SET = c("Mouse_Bulk_Heterogeneity")

for( EXPERIMENT_NAME in EXPERIMENT_NAME_SET){

	### Load parameters
	# Define an environment that will contain parameters
	paramsEnv = new.env();

	# Load file defining global parameters
	globalParamsFilePath = file.path( WORKING_DIR, "../globalParams.R");
	if(file.exists(globalParamsFilePath)) {
	  source( globalParamsFilePath, local = paramsEnv);
	} else {
	  warning("The file 'globalParamsFilePath.R' containing global parameters is missing.");
	}

	# Load file defining analysis parameters
	analysisParamsFilePath = file.path( WORKING_DIR, "analysisParams.R");
	if(file.exists(analysisParamsFilePath)) {
	  source( analysisParamsFilePath, local = paramsEnv);
	} else {
	  warning("The file 'analysisParamsFilePath.R' containing analysis-specific parameters is missing.");
	}

	# # # Load file defining sample parameters
	# sampleParamsFilePath = file.path( WORKING_DIR, "../04_Config", paste0( EXPERIMENT_NAME, "_sampleParams.R"));
	# if(file.exists(sampleParamsFilePath)) {
	#   source( sampleParamsFilePath, local = paramsEnv);
	# } else {
	#   warning("The file 'sampleParamsFilePath.R' containing sample parameters is missing.");
	# }

	# Assign loaded values to current environment (Not .GlobalEnv as it may not be the current one depending on where rendering is started from)
	invisible( lapply( ls( paramsEnv, all.names = TRUE), function(x, envir)
	{ 
	  assign( x = x, 
		  value = get( x, pos = paramsEnv), 
		  pos = envir)
	}, 
	environment()));

	# Loop over the tissue to make one report per tissue

	PATH_EXPERIMENT_ANALYSIS_OUTPUT = file.path ("/media/owen/Backup Plus/UCSF_Project/02_PROCESSED_DATA/PNOC/ANALYSIS/09_Bulk_RNA/05_Output/Mouse/")

	rmarkdown::render( input = file.path( WORKING_DIR, "Bulk_RNA_seq_Report.Rmd"),
		             output_dir = PATH_EXPERIMENT_ANALYSIS_OUTPUT,
		             #output_file  = paste0( PROJECT_SHORT_NAME, "_", ANALYSIS_STEP_NAME, "_", EXPERIMENT_NAME, ".html"),
		             quiet = FALSE)
	
	rm(list = ls( paramsEnv))
}

