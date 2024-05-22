#' Generate Summary Statistics
#'
#' This function generates summary statistics for metabolite data including sample and feature statistics.
#'
#' @param mydata List containing metabolite, sample, and feature data.
#' @param VAR_data_dir Directory where summary statistics will be saved.
#' @param VAR_project Name of the project.
#' @param VAR_outlier_udist Threshold for defining outliers.
#' @param VAR_derived_var_exclusion Whether to exclude derived variables.
#' @param VAR_tree_cut_height Height for hierarchical clustering cut.
#'
#' @return List containing sample and feature summary statistics.
#'
#' @export
generate_summary_statistics <- function(mydata, VAR_data_dir, VAR_project, VAR_outlier_udist, VAR_derived_var_exclusion, VAR_tree_cut_height) {
  cat(paste0("3. Summary statistics \n"))
  
  ## Sample summary statistics
  cat(paste0("\ta. Sample summary statistics \n"))
  
  if (length(mydata$featuredata$SUPER_PATHWAY) > 0) {
    w <- which(mydata$featuredata$SUPER_PATHWAY %in% c("xenobiotics", "Xenobiotics"))
    xeno_names <- mydata$featuredata$feature_names[w]
    samplesumstats <- sample.sum.stats(wdata = mydata$metabolitedata, feature_names_2_exclude = xeno_names, outlier_udist = VAR_outlier_udist)
  } else {
    if (length(mydata$featuredata$derived_features) > 0 & VAR_derived_var_exclusion == "TRUE") {
      w <- which(mydata$featuredata$derived_features == "yes")
      derivedfeature_names <- as.character(mydata$featuredata$feature_names[w])
      samplesumstats <- sample.sum.stats(wdata = mydata$metabolitedata, feature_names_2_exclude = derivedfeature_names, outlier_udist = VAR_outlier_udist)
    } else {
      samplesumstats <- sample.sum.stats(wdata = mydata$metabolitedata, outlier_udist = VAR_outlier_udist)
    }
  }
  
  ### Write sample sum stats to file
  cat(paste0("\t\t- Writing sample summary statistics to file.\n"))
  
  ## Make a sum stats directory in VAR_data_dir
  dd <- gsub(" ", "\\\\ ", VAR_data_dir)
  system(paste0("mkdir -p ", dd, "metaboprep_release_", today, "/", "sumstats"))
  ## Make a raw_dataset directory inside the sumstats folder
  system(paste0("mkdir -p ", dd, "metaboprep_release_", today, "/", "sumstats/raw_dataset"))
  
  ### SAMPLES
  if ("sampledata" %in% names(mydata)) {
    samplesumstats <- cbind(mydata$sampledata, samplesumstats[,-1])
  }
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/raw_dataset/", VAR_project, "_", today, "_sample_anno_sumstats.txt")
  write.table(samplesumstats, file = n,
              row.names = FALSE, col.names = TRUE, 
              sep = "\t", quote = FALSE)
  
  ## Feature summary statistics
  cat(paste0("\tb. Feature summary statistics \n"))
  
  ### Sample missingness
  if (length(samplesumstats$VAR_sample_missingness_w_exclusions) > 0) {
    sammis <- samplesumstats$VAR_sample_missingness_w_exclusions
  } else {
    sammis <- samplesumstats$VAR_sample_missingness
  }
  
  ### Feature names to exclude
  if (length(mydata$featuredata$derived_features) > 0 & VAR_derived_var_exclusion == "TRUE") {
    w <- which(mydata$featuredata$derived_features == "yes") 
    fn2e <- as.character(mydata$featuredata$feature_names[w])
  } else {
    fn2e <- NA
  }
  
  ### Run feature summary stats function
  featuresumstats <- feature.sum.stats(wdata = mydata$metabolitedata,
                                       sammis = sammis, 
                                       tree_cut_height = VAR_tree_cut_height,
                                       outlier_udist = VAR_outlier_udist,
                                       feature_names_2_exclude = fn2e)
  
  ### Write feature sum stats to file
  cat(paste0("\t\t- Writing feature summary statistics to file.\n"))
  
  if ("featuredata" %in% names(mydata)) {
    featuresumstats$table <- cbind(mydata$featuredata, featuresumstats$table[,-1])
  }
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/raw_dataset/", VAR_project, "_", today, "_feature_anno_sumstats.txt")
  write.table(featuresumstats$table, file = n,
              row.names = FALSE, col.names = TRUE, 
              sep = "\t", quote = TRUE)
  
  ## PC outliers
  cat(paste0("\tc. Performing principle component analysis and identifying outliers.\n"))
  
  ## Identify independent feature names as reported in featuresumstats
  w <- which(featuresumstats$table$independent_features_binary == 1)
  indf <- as.character(featuresumstats$table[w, 1])
  if (sum(indf %in% colnames(mydata$metabolitedata)) == 0) {
    indf <- as.character(featuresumstats$table$feature_names[w])
  }
  PCs_outliers <- pc.and.outliers(metabolitedata = mydata$metabolitedata, indfeature_names = indf)
  
  ### Re-write sample sum stats to file to include PCs
  cat(paste0("\t\t- Re-Writing sample summary statistics to file to include PCs.\n"))
  
  samplesumstats <- cbind(samplesumstats, PCs_outliers[[1]])
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/raw_dataset/", VAR_project, "_", today, "_sample_anno_sumstats.txt")
  write.table(samplesumstats, file = n,
              row.names = FALSE, col.names = TRUE, 
              sep = "\t", quote = FALSE)
  
  ### Write the variance explained by PCs out to file
  cat(paste0("\t\t- Writing PC statistics to file.\n\n"))
  
  varexp <- data.frame(VarExp = PCs_outliers[[2]])
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/raw_dataset/", VAR_project, "_", today, "_pc_varexp.txt")
  write.table(varexp, file = n,
              row.names = TRUE, col.names = TRUE, 
              sep = "\t", quote = FALSE)
  
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/raw_dataset/", VAR_project, "_", today, "_featuretree.Rdata")
  feature_tree <- featuresumstats[[2]]
  save(feature_tree, file = n)
  
  return(list(featuresumstats = featuresumstats, samplesumstats = samplesumstats, varexp = varexp, feature_tree = feature_tree))
}

#' Perform Quality Control on Metabolite Data
#'
#' This function performs quality control on metabolite data, including filtering and the exclusion of specified variables. It returns two lists: `qcdata` and `qcing_data`.
#'
#' @param mydata A list containing the raw metabolite data, sample data, and feature data.
#' @param VAR_data_dir A string specifying the directory where output files will be saved.
#' @param VAR_project A string specifying the project name used in output file names.
#' @param VAR_derived_var_exclusion A string indicating whether derived variables should be excluded ("TRUE" or "FALSE").
#' @param VAR_feature_missingness A numeric value specifying the threshold for feature missingness.
#' @param VAR_sample_missingness A numeric value specifying the threshold for sample missingness.
#' @param VAR_total_peak_area_SD A numeric value specifying the standard deviation threshold for total peak area.
#' @param VAR_outlier_treatment A string indicating the method of outlier treatment ("leave_be", "turn_NA", "winsorize").
#' @param VAR_outlier_udist A numeric value specifying the threshold for outlier distance.
#' @param VAR_PC_outlier_SD A numeric value specifying the standard deviation threshold for principal component outliers.
#' @param VAR_tree_cut_height A numeric value specifying the height for tree cutting in hierarchical clustering.
#' @return A list containing two lists: `qcdata` and `qcing_data`.
#' @examples
#' \dontrun{
#' result <- perform_qc(
#'   mydata = mydata, 
#'   VAR_data_dir = "output_directory/", 
#'   VAR_project = "my_project", 
#'   VAR_derived_var_exclusion = "TRUE", 
#'   VAR_feature_missingness = 0.1, 
#'   VAR_sample_missingness = 0.1, 
#'   VAR_total_peak_area_SD = 3, 
#'   VAR_outlier_treatment = "turn_NA", 
#'   VAR_outlier_udist = 3, 
#'   VAR_PC_outlier_SD = 2, 
#'   VAR_tree_cut_height = 0.5
#' )
#' qcdata <- result$qcdata
#' qcing_data <- result$qcing_data
#' }
#' @export
perform_qc <- function(mydata, VAR_data_dir, VAR_project, VAR_derived_var_exclusion, VAR_feature_missingness, VAR_sample_missingness, VAR_total_peak_area_SD, VAR_outlier_treatment, VAR_outlier_udist, VAR_PC_outlier_SD, VAR_tree_cut_height) {
  ### xenobiotics to exclude
  w <- which(mydata$featuredata$SUPER_PATHWAY %in% c("xenobiotics", "Xenobiotics"))
  xeno_names <- mydata$featuredata$feature_names[w]
  if (length(xeno_names) == 0) { xeno_names <- NA }
  
  ### derived variables to exclude
  w <- which(mydata$featuredata$derived_features == "yes")
  derived_names <- as.character(mydata$featuredata$feature_names[w])
  if (length(derived_names) == 0) { derived_names <- NA }
  if (VAR_derived_var_exclusion != "TRUE") { derived_names <- NA }
  
  ### execute super function
  cat(paste0("\ta. Performing data filtering.\n"))
  
  dataQC <- perform.metabolite.qc(
    wdata = mydata$metabolitedata,
    fmis = VAR_feature_missingness,
    smis = VAR_sample_missingness,
    tpa_out_SD = VAR_total_peak_area_SD,
    outlier_treatment = VAR_outlier_treatment, ## options are "leave_be", "turn_NA", "winsorize"
    winsorize_quantile = 1,                ## winsorize to what quantile of remaining (not outlier) values
    outlier_udist = VAR_outlier_udist,
    PC_out_SD = VAR_PC_outlier_SD,
    tree_cut_height = VAR_tree_cut_height,
    feature_colnames_2_exclude = xeno_names,
    derived_colnames_2_exclude = derived_names
  )
  
  #################################
  ## B. write QC data to file
  #################################
  cat(paste0("\tb. Writing QC data to file.\n\n"))
  
  #############################
  ##
  ## B.1. Make a QCing data set object
  ##
  #############################
  qcing_data <- list(
    metabolite_data = dataQC$wdata,
    exclusion_data = dataQC$exclusion_data,
    feature_sumstats = dataQC$featuresumstats$table,
    feature_tree = dataQC$featuresumstats$tree,
    pcs = dataQC$pca$pcs,
    varexp = dataQC$pca$varexp,
    accelerationfactor = dataQC$pca$accelerationfactor,
    nparallel = dataQC$pca$nsig_parrallel
  )
  
  ##################################
  ## B.2. Add sample and feature data to qcdata
  ##################################
  temp <- dataQC$wdata
  m <- match(rownames(temp), mydata$sampledata[, 1])
  n <- match(colnames(temp), mydata$featuredata[, 1])
  if (length(n) == sum(is.na(n))) {
    n <- match(colnames(temp), mydata$featuredata$feature_names)
  }
  qcdata <- list(
    metabolitedata = temp,
    sampledata = as.data.frame(mydata$sampledata[m, ]),
    featuredata = as.data.frame(mydata$featuredata[n, ])
  )
  rm(temp)
  
  if (colnames(qcdata$sampledata)[1] == "mydata$sampledata[m, ]") { colnames(qcdata$sampledata)[1] = "SampleID" }
  if (colnames(qcdata$featuredata)[1] == "mydata$featuredata[n, ]") { colnames(qcdata$featuredata)[1] = "feature_names" }
  
  ##################################
  ## B.3. Write to file
  ##################################
  ## Create the directory for filtered data if it doesn't exist
  dir.create(paste0(VAR_data_dir, "metaboprep_release_", today, "/filtered_data"), recursive = TRUE, showWarnings = FALSE)
  
  ## qc metabolite data
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/filtered_data/", VAR_project, "_", today, "_Filtered_metabolite_data.txt")
  write.table(qcdata$metabolitedata, file = n,
              row.names = TRUE, col.names = TRUE, 
              sep = "\t", quote = FALSE)
  
  ## qc sample data
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/filtered_data/", VAR_project, "_", today, "_Filtered_sample_data.txt")
  write.table(qcdata$sampledata, file = n,
              row.names = FALSE, col.names = TRUE, 
              sep = "\t", quote = FALSE)
  
  ## qc feature data
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/filtered_data/", VAR_project, "_", today, "_Filtered_feature_data.txt")
  write.table(qcdata$featuredata, file = n,
              row.names = FALSE, col.names = TRUE, 
              sep = "\t", quote = FALSE)
  
  return(list(qcdata = qcdata, qcing_data = qcing_data))
}

#' Estimate Summary Statistics on the QC'd Data Set
#'
#' This function estimates summary statistics on the filtered QC'd data set, including summary statistics for samples and features, and performs principal component analysis.
#'
#' @param qcdata A list containing the QC'd metabolite data, sample data, and feature data.
#' @param VAR_data_dir A string specifying the directory where output files will be saved.
#' @param VAR_project A string specifying the project name used in output file names.
#' @param VAR_derived_var_exclusion A string indicating whether derived variables should be excluded ("TRUE" or "FALSE").
#' @param VAR_outlier_udist A numeric value specifying the threshold for outlier distance.
#' @param VAR_tree_cut_height A numeric value specifying the height for tree cutting in hierarchical clustering.
#' @return A list containing the final filtered data set, including metabolite data, sample data, feature data, feature tree, variance explained by principal components, and other statistics.
#' @examples
#' \dontrun{
#' qc_data <- estimate_summary_statistics(
#'   qcdata = qcdata, 
#'   VAR_data_dir = "output_directory/", 
#'   VAR_project = "my_project", 
#'   VAR_derived_var_exclusion = "TRUE", 
#'   VAR_outlier_udist = 3, 
#'   VAR_tree_cut_height = 0.5
#' )
#' }
#' @export
estimate_summary_statistics <- function(qcdata, VAR_data_dir, VAR_project, VAR_derived_var_exclusion, VAR_outlier_udist, VAR_tree_cut_height) {
## Estimate  Summary Statistics on the QC'd Data Set
  cat(paste0("VI. Estimating Summary Statistics on Filtered Data Set.\n"))
  ## A. Summary Statistics for samples
  cat(paste0("\ta. Estimating summary statistics for Filtered samples\n"))
  
  ## A.1. Estimate sum stats
  ## Is this metabolon data?? 
  ##  -- is the column SUPER_PATHWAY present in the feature data
  ##  -- if yes, exclude Xenobiotics from one of the missingness estimate
  if (length(qcdata$featuredata$SUPER_PATHWAY) > 0) {
    w <- which(qcdata$featuredata$SUPER_PATHWAY %in% c("xenobiotics", "Xenobiotics"))
    xeno_names <- rownames(qcdata$featuredata)[w]
    samplesumstats <- sample.sum.stats(wdata = qcdata$metabolitedata, feature_names_2_exclude = xeno_names, outlier_udist = VAR_outlier_udist)
  } else {
    ## Is this Nightingale data?? 
    ##  -- is the column derived_features present in the feature data
    ##  -- if yes, exclude derived variables from one of the missingness estimate
    if (length(qcdata$featuredata$derived_features) > 0 & VAR_derived_var_exclusion == "TRUE") {
      w <- which(qcdata$featuredata$derived_features == "yes")
      derivedfeature_names <- as.character(qcdata$featuredata$feature_names[w])
      samplesumstats <- sample.sum.stats(wdata = qcdata$metabolitedata, feature_names_2_exclude = derivedfeature_names, VAR_outlier_udist = VAR_outlier_udist)
    } else {
      samplesumstats <- sample.sum.stats(wdata = qcdata$metabolitedata, outlier_udist = VAR_outlier_udist)
    }
  }
  
  ### A.2. Write sample sum stats to file
  cat(paste0("\tb. Writing filtered sample summary statistics to file.\n"))
  
  ### WRITE
  if ("sampledata" %in% names(qcdata)) {
    ## add sample stats to the sample annotation file
    samplesumstats <- cbind(qcdata$sampledata, samplesumstats[,-1])
  }
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/filtered_dataset/", VAR_project, "_", today, "_sample_anno_sumstats.txt")
  write.table(samplesumstats, file = n, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  ##################################
  ## B. Summary Statistics for features
  ##################################
  cat(paste0("\tc. Estimating summary statistics for filtered features.\n"))
  
  ### sample missingness
  if (length(samplesumstats$VAR_sample_missingness_w_exclusions) > 0) {
    sammis <- samplesumstats$VAR_sample_missingness_w_exclusions
  } else {
    sammis <- samplesumstats$VAR_sample_missingness
  }
  
  ### feature names to exclude
  if (length(qcdata$featuredata$derived_features) > 0 & VAR_derived_var_exclusion == "TRUE") {
    w <- which(qcdata$featuredata$derived_features == "yes")
    fn2e <- as.character(qcdata$featuredata$feature_names[w])
  } else {
    fn2e <- NA
  }
  
  ### RUN feature summary stats function
  featuresumstats <- feature.sum.stats(wdata = qcdata$metabolitedata, sammis = sammis, tree_cut_height = VAR_tree_cut_height, outlier_udist = VAR_outlier_udist, feature_names_2_exclude = fn2e)
  
  ## count of independent features
  icount <- sum(featuresumstats$table$independent_features_binary)
  cat(paste0("\t\t\t- A total of ", icount ," independent features were identified in the total filtered data set.\n"))
  
  ### write feature sum stats to file
  cat(paste0("\td. Writing feature summary statistics to file.\n"))
  
  if ("featuredata" %in% names(qcdata)) {
    ## add feature stats to the feature annotation file
    featuresumstats$table <- cbind(qcdata$featuredata, featuresumstats$table[, -1])
  }
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/filtered_dataset/", VAR_project, "_", today, "_feature_anno_sumstats.txt")
  write.table(featuresumstats$table, file = n, row.names = FALSE, col.names = TRUE, sep = "\t", quote = TRUE)
  
  ##################################
  ## C. Generation of PCs
  ##################################
  cat(paste0("\te. Performing principal component analysis on final filtered data set.\n"))
  
  ## identify independent feature names as reported in featuresumstats
  w <- which(featuresumstats$table$independent_features_binary == 1)
  indf <- featuresumstats$table[w, 1]
  if (sum(indf %in% colnames(mydata$metabolitedata)) == 0) {
    indf <- as.character(featuresumstats$table$feature_names[w])
  }
  PCs_outliers <- pc.and.outliers(metabolitedata = qcdata$metabolitedata, indfeature_names = indf)
  
  cat(paste0("\t\t The number of informative principal components:\n"))
  cat(paste0("\t\t\t 1. Cattel's Scree Test : acceleration factor = ", PCs_outliers$accelerationfactor, "\n"))
  cat(paste0("\t\t\t 2. Parallel Analysis = ", PCs_outliers$nsig_parrallel, "\n"))
  
  ### write sample sum stats to file
  cat(paste0("\tf. Re-Writing filtered sample summary statistics to file to include PCs.\n"))
  
  ### SAMPLES
  samplesumstats <- cbind(samplesumstats, PCs_outliers[[1]][, 1:10])
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/filtered_dataset/", VAR_project, "_", today, "_sample_anno_sumstats.txt")
  write.table(samplesumstats, file = n, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  ### write the variance explained by PCs out to file
  cat(paste0("\tg. Writing PC statistics to file.\n\n"))
  
  varexp <- data.frame(VarExp = PCs_outliers[[2]])
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/filtered_dataset/", VAR_project, "_", today, "_pc_varexp.txt")
  write.table(varexp, file = n, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/sumstats/filtered_dataset/", VAR_project, "_", today, "_featuretree.Rdata")
  feature_tree <- featuresumstats[[2]]
  save(feature_tree, file = n)
  
  #############################
  ##
  ## Make a FILTERED data set object
  ##
  #############################
  qc_data <- list(
    metabolite_data = qcdata$metabolitedata,
    sample_data = samplesumstats,
    feature_data = featuresumstats$table,
    feature_tree = feature_tree,
    varexp = varexp,
    accelerationfactor = PCs_outliers$accelerationfactor,
    nparallel = PCs_outliers$nsig_parrallel
  )
  
  #########################
  ##
  ## (IX) Save Data
  ##
  #########################
  cat(paste0("VII. Generate Data Description html report.\n"))
  
  ############
  n <- paste0(VAR_data_dir, "metaboprep_release_", today, "/ReportData.Rdata")
  save(raw_data, qcing_data, qc_data, VAR_project, VAR_platform, VAR_data_dir, 
       VAR_feature_missingness, VAR_sample_missingness, VAR_total_peak_area_SD, VAR_PC_outlier_SD, VAR_tree_cut_height, file = n)
  ############
  
  return(qc_data)
}
