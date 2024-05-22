rm(list=ls())
set.seed(821)

# environment ====
library(metaboprep)

# VARS ====
## record todays date
today = Sys.Date()
today = gsub("-","_",today)

# arguments ====
## single argument = the path to the paramater_file.txt
args <- commandArgs(trailingOnly=TRUE)
## Check that you passed an argument to the script
if(length(args) != 1){
  stop( 
    paste0("You must provide a single argument (i.e., the full path to the paramater file) when you call the script."),
    call.=FALSE)
} 

# paramaters ====
VAR_pfile <- read.table(args[1], header = FALSE, sep = "=", as.is = TRUE)
## project name
VAR_project <- as.character(VAR_pfile[1,2] )
## path	
VAR_data_dir <- as.character(VAR_pfile[2,2] )
if (Sys.info()[['sysname']] == 'Windows') {
  if (substring(VAR_data_dir, nchar(VAR_data_dir)) != "\\") {
    VAR_data_dir <- paste0(VAR_data_dir, "\\")
  }
} else {
  if (substring(VAR_data_dir, nchar(VAR_data_dir)) != "/") {
    VAR_data_dir <- paste0(VAR_data_dir, "/")
  }
}
## identify file type
VAR_METABO_file2process <- as.character(VAR_pfile[3,2] )
VAR_xl_ftype <- c(".xls",".xlsx", ".XLS", ".XLSX")
VAR_isexcel <- sum(unlist(sapply(VAR_xl_ftype, function(x){ grep(x, VAR_METABO_file2process) } ) ) )
VAR_flat_ftype <- c(".csv",".txt", ".tsv", ".TXT", ".CSV", ".TSV")
VAR_isflat <- sum(unlist(sapply(VAR_flat_ftype, function(x){ grep(x, VAR_METABO_file2process) } ) ) )
## annotation file
VAR_FeatureAnno_file2process <- as.character(VAR_pfile[4,2] )
VAR_SampleAnno_file2process <- as.character(VAR_pfile[5,2] )
## platform
VAR_platform <- as.character(VAR_pfile[6,2] )
# QC variables
VAR_feature_missingness <- as.numeric(VAR_pfile[7,2] )
VAR_sample_missingness <- as.numeric(VAR_pfile[8,2] )
VAR_total_peak_area_SD <- as.numeric(VAR_pfile[9,2] )
VAR_outlier_udist <- as.numeric(VAR_pfile[10,2] )
VAR_outlier_treatment <- as.character(VAR_pfile[11,2] )
VAR_tree_cut_height <- as.numeric(VAR_pfile[12,2] )
VAR_PC_outlier_SD <- as.numeric(VAR_pfile[13,2] )
## data
VAR_derived_var_exclusion <- VAR_pfile[14,2] 
VAR_feat_anno_run_mode_col <- as.character(VAR_pfile[15,2] )
VAR_plot_feature_distributions <- VAR_pfile[16,2] 
VAR_data <- paste0(VAR_data_dir, VAR_METABO_file2process)

# setup ====
cat(paste0("1. Setting up your pipeline\n"))
cat(paste0("\t - Creating directory\n"))

## directories ====
### main directory
subdir <- paste0(VAR_data_dir, "metaboprep_release_", today)
if (Sys.info()[['sysname']] == 'Windows') {
  subdir <- gsub("/", "\\\\", subdir)  # Ensure backslashes on Windows
}
if (!dir.exists(subdir)) {
  dir.create(subdir, recursive = TRUE)
}
### sumstats directory
sumstats_dir <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats")
if (Sys.info()[['sysname']] == 'Windows') {
  sumstats_dir <- gsub("/", "\\\\", sumstats_dir)  # Ensure backslashes on Windows
}
if (!dir.exists(sumstats_dir)) {
  dir.create(sumstats_dir, recursive = TRUE)
}
### filtered directory
filtered_data_dir <- paste0(VAR_data_dir, "metaboprep_release_", today, "/filtered_data")
if (Sys.info()[['sysname']] == 'Windows') {
  filtered_data_dir <- gsub("/", "\\\\", filtered_data_dir)  # Ensure backslashes on Windows
}
if (!dir.exists(filtered_data_dir)) {
  dir.create(filtered_data_dir, recursive = TRUE)
}
### filtered dataset
filtered_dataset_dir <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/filtered_dataset")
if (Sys.info()[['sysname']] == 'Windows') {
  filtered_dataset_dir <- gsub("/", "\\\\", filtered_dataset_dir)  # Ensure backslashes on Windows
}
if (!dir.exists(filtered_dataset_dir)) {
  dir.create(filtered_dataset_dir, recursive = TRUE)
}
### raw directory
raw_dataset_dir <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/raw_dataset")
if (Sys.info()[['sysname']] == 'Windows') {
  raw_dataset_dir <- gsub("/", "\\\\", raw_dataset_dir)  # Ensure backslashes on Windows
}
if (!dir.exists(raw_dataset_dir)) {
  dir.create(raw_dataset_dir, recursive = TRUE)
}

## log file ====
cat(paste0("\t - Creating log\n"))
logfilename = paste0(VAR_data_dir, "metaboprep_release_", today, "/", VAR_project, "_", today,  "_logfile.txt")
sink(file = logfilename , split = TRUE  )

## cat() ====
cat(paste0("Thresholds:"))
cat(paste0("\t- Feature extreme missingness: 0.8 \n"))
cat(paste0("\t- Sample extreme missingness: 0.8 \n"))
cat(paste0("\t- Feature missingness: ", VAR_feature_missingness, "\n"))
cat(paste0("\t- Sample missingness: ", VAR_sample_missingness, "\n"))
cat(paste0("\t- Total peak area: ", VAR_total_peak_area_SD, "SD \n"))
cat(paste0("\t- Interquartile range distance from the median of each feature: ", VAR_outlier_udist, "SD \n"))
cat(paste0("\t- Outliers for PCA: ", VAR_outlier_treatment, "\n"))
cat(paste0("\t- Metabolite independence, Spearmans rho: ", VAR_tree_cut_height, "\n"))
cat(paste0("\t- PCA exclusion: ", VAR_PC_outlier_SD, " SD\n"))

cat(paste0("Data:"))
cat(paste0("\t- Platform: ", VAR_platform, "\n"))
cat(paste0("\t- Directory: ", VAR_data_dir, "\n"))
cat(paste0("\t- Data: ", VAR_METABO_file2process, "\n"))

cat(paste0("Arguments:"))
## file type
if(VAR_isexcel > 0){
  cat(paste0("\t- Feature file is excel and will be processed as a commercial source file \n"))  
} else {
  if(VAR_isflat > 0){
    cat(paste0("\t- Feature file is flat text file.\n"))    
  } else {
    stop(paste0("\t- Feature file type not detected, should be xls, xlsx, txt, or csv \n"), call.=FALSE )    
  } 
}
## annotation file
if( !is.na(VAR_FeatureAnno_file2process)){
  cat(paste0("\t- Feature annotation file: ", VAR_FeatureAnno_file2process, "\n"))  
} else{
  cat(paste0("\t- Feature annotation file: not provided \n"))  
}
## sample/batch annotation file
if( !is.na(VAR_SampleAnno_file2process)){
  cat(paste0("\t- Sample/batch annotation file: ", VAR_SampleAnno_file2process, "\n"))  
} else{
  cat(paste0("\t- Sample/batch annotation file: not provided \n"))  
}
# annotation files
## Nightingale derived variable exclusions when evaluting SAMPLE quality for QC
if(VAR_platform == "Nightingale"){
  if(VAR_derived_var_exclusion == TRUE){
    cat(paste0("\t- Nightingale derived variables: will be excluded from data filtering steps \n"))
  }
}
## Mass Spec Run Mode for each metabolite; This variable defines the column name, in the feature_annotation_file, indexing the run mode string(s). There should in turn be column name(s) in the sample_annotation_file that match the run mode string(s) and hold the batch IDs for each sample, in that run mode. 
if( !is.na(VAR_feat_anno_run_mode_col)){
  cat(paste0("\t- Run mode column name: ", VAR_feat_anno_run_mode_col, "\n"))  
} 
## visual report
if (VAR_plot_feature_distributions == TRUE) {
  cat(paste0("\t- Create visual report: yes\n"))
} else {
  cat(paste0("\t- Create visual report: no\n"))
}

# data ====
cat(paste0("2. Data \n"))
if(VAR_isflat > 0){
  metabolitedata <- read_data(VAR_data = VAR_data, VAR_isflat = VAR_isflat, VAR_platform = VAR_platform)
}
if(!is.na(VAR_FeatureAnno_file2process)){
  featuredata <- feature_annotation(VAR_data_dir = VAR_data_dir, VAR_FeatureAnno_file2process = VAR_FeatureAnno_file2process, VAR_platform = VAR_platform, metabolitedata = metabolitedata)
}
if(!is.na(VAR_FeatureAnno_file2process)){
  sampledata <- sample_annotation(VAR_data_dir = VAR_data_dir, VAR_SampleAnno_file2process = VAR_SampleAnno_file2process, metabolitedata = metabolitedata)
}

# generate working data ====
if(VAR_isflat > 0){
  mydata <- generate_sample_and_feature_data(VAR_isflat = VAR_isflat, metabolitedata = metabolitedata, ng_anno = ng_anno, VAR_platform = VAR_platform)
}
if(VAR_isexcel > 0){
  mydata <- generate_sample_and_feature_data_excel(VAR_isexcel = VAR_isexcel, VAR_platform = VAR_platform, VAR_METABO_file2process = VAR_METABO_file2process, VAR_data_dir = VAR_data_dir, VAR_project = VAR_project)
}

# normalize data ====
if(VAR_platform == "Metabolon"){
  mydata <- normalize_ms_data(VAR_platform = VAR_platform, mydata = mydata, VAR_feat_anno_run_mode_col = VAR_feat_anno_run_mode_col)
}

# summary stats ====
result <- generate_summary_statistics(mydata = mydata, VAR_data_dir = VAR_data_dir, VAR_project = VAR_project, VAR_outlier_udist = VAR_outlier_udist, VAR_derived_var_exclusion = VAR_derived_var_exclusion, VAR_tree_cut_height = VAR_tree_cut_height)
featuresumstats <- result$featuresumstats
samplesumstats <- result$samplesumstats
varexp <- result$varexp
feature_tree <- result$feature_tree

raw_data = list(metabolite_data = mydata$metabolitedata,
                sample_data = samplesumstats,
                feature_data = featuresumstats$table,
                feature_tree = feature_tree,
                varexp = varexp)

# perform QC ====
qcdata <- perform_qc(
  mydata = mydata, 
  VAR_data_dir = VAR_data_dir, 
  VAR_project = VAR_project, 
  VAR_derived_var_exclusion = VAR_derived_var_exclusion, 
  VAR_feature_missingness = VAR_feature_missingness, 
  VAR_sample_missingness = VAR_sample_missingness, 
  VAR_total_peak_area_SD = VAR_total_peak_area_SD, 
  VAR_outlier_treatment = VAR_outlier_treatment, 
  VAR_outlier_udist = VAR_outlier_udist, 
  VAR_PC_outlier_SD = VAR_PC_outlier_SD, 
  VAR_tree_cut_height = VAR_tree_cut_height
)
qcing_data <- qcdata$qcing_data
qcdata <- qcdata$qcdata

# estimate summary stats ====
qc_data <- estimate_summary_statistics(qcdata = qcdata, VAR_data_dir = VAR_data_dir, VAR_project = VAR_project, VAR_derived_var_exclusion = VAR_derived_var_exclusion, VAR_outlier_udist = VAR_outlier_udist, VAR_tree_cut_height = VAR_tree_cut_height)

# log ====
## stop writing to log file.
sink()

# report ====
output_dir_path = paste0(VAR_data_dir, "metaboprep_release_", today, "/")
rdfile = paste0(output_dir_path, "ReportData.Rdata")
generate_report(full_path_2_Rdatafile = rdfile, 
                dir_4_report = output_dir_path)

# feature plots ====
if(VAR_plot_feature_distributions == TRUE){
  cat( paste0("VIII. Plot feature distributions and summary statistics to pdf.\n") )
  
  feature_plots(wdata = raw_data$metabolite_data, outdir = output_dir_path, nsd = VAR_PC_outlier_SD)
  
  f = paste0(output_dir_path, "feature_distribution.pdf")
  cat( paste0("\t- plots for each figure written to the pdf ", f,".\n") )
}
