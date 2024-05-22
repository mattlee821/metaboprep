#' Normalize MS Data
#'
#' This function performs normalization on Metabolon or other MS data based on the platform.
#'
#' @param VAR_platform Platform of the data (e.g., "Metabolon").
#' @param mydata List containing metabolite, sample, and feature data.
#' @param VAR_feat_anno_run_mode_col Name of the feature annotation column for normalization.
#'
#' @return Updated data list with normalized metabolite data.
#'
#' @export
normalize_ms_data <- function(VAR_platform, mydata, VAR_feat_anno_run_mode_col) {
  if(VAR_platform == "Metabolon") {
    cat(paste0("Normalization. Performing normalization on Metabolon Data.\n\n"))
    if(!is.na(VAR_feat_anno_run_mode_col)) {
      cat(paste0("\t- Performing normalization with parameter file provided feature annotation column '",VAR_feat_anno_run_mode_col,"'.\n"))
      
      norm_metabolite_data <- batch_normalization( 
        wdata = mydata$metabolitedata, 
        feature_data_sheet =  mydata$featuredata, 
        sample_data_sheet = mydata$sampledata, 
        feature_runmode_col = VAR_feat_anno_run_mode_col, 
        batch_ids = NULL  
      )
      
      mydata$raw_metabolitedata <- mydata$metabolitedata
      mydata$metabolitedata <- norm_metabolite_data
      rm(norm_metabolite_data)
      cat( paste0("\t- Normalization completed.\n\n") )
    } else {
      cat( paste0("\t- Looking for run mode information automatically given Metabolon standard data release formatting.\n") )
      fanno_col_number <- which(colnames(mydata$featuredata) %in% c("VAR_platform","VAR_platform"))
      if(length(fanno_col_number) == 1) {
        runmode <- unlist( mydata$featuredata[,fanno_col_number] )
        runmode <- gsub("LC\\/MS\\ ","",runmode)
        runmode <- gsub(" ","",runmode)
        runmode <- tolower(runmode)
        
        mydata$featuredata[,fanno_col_number] <- runmode
        batchrunmodes <- unique(runmode)
      } else {
        cat(paste0("\t- NOTE: Unable to identify a column header called 'VAR_platform' or 'VAR_platform' in the feature data sheet. This is necessary to perform batch normalization\n\n") )
      }
      
      n <- tolower(colnames(mydata$sampledata))
      n <- gsub(" ","",n)
      n <- gsub("_","",n)
      n <- gsub("\\.","",n)
      
      k <- which(n %in% runmode)
      colnames(mydata$sampledata)[k] <- n[k]
      
      if(length(k) == length(batchrunmodes)) {
        cat( paste0("\t- Performing normalization with automatically identified Metabolon standard data release formatting column '",colnames(mydata$featuredata)[fanno_col_number],"'.\n") )
        
        norm_metabolite_data <- batch_normalization( 
          wdata = mydata$metabolitedata, 
          feature_data_sheet =  mydata$featuredata, 
          sample_data_sheet = mydata$sampledata, 
          feature_runmode_col = fanno_col_number, 
          batch_ids = NULL  
        )
        
        mydata$raw_metabolitedata <- mydata$metabolitedata
        mydata$metabolitedata <- norm_metabolite_data
        rm(norm_metabolite_data)
        cat( paste0("\t- Normalization completed.\n\n") )
      } else {
        cat(paste0("\t- NOTE: Unable to identify a column headers in sample sheet that match the VAR_platform run modes found in the feature data sheet. This should be something like neg, polar, pos early, pos late.\n") )
        
        cat( paste0("\t- NOTE: We will take the ScaledImpData data and remove the imputed data to extract the normalized data.\n\n") )
        
        scaled_imputed_data_tab <- grep("ScaledImp", names(mydata))
        if(length(scaled_imputed_data_tab) == 1) {
          ndata <- mydata[[scaled_imputed_data_tab]]
          
          if(sum(is.na(mydata$metabolitedata)) > 0) {
            ndata[is.na(mydata$metabolitedata)] <- NA  
          }
          mydata$raw_metabolitedata <- mydata$metabolitedata
          mydata$metabolitedata <- ndata 
          rm(ndata) 
          cat( paste0("\t- Normalization completed.\n\n") )
        } else {
          cat( paste0("\t- NOTE: Unable to identify a 'ScaledImp' data tab in the excel file. No normalization carried out.\n\n") )
        }    
      }
    }
  }
  
  if(!is.na(VAR_feat_anno_run_mode_col) & VAR_platform == "Other") {
    cat( paste0("Normalization. Performing normalization parameter file provided feature annotation column '",VAR_feat_anno_run_mode_col,"'.\n") )
    
    norm_metabolite_data <- batch_normalization( 
      wdata = mydata$metabolitedata, 
      feature_data_sheet =  mydata$featuredata, 
      sample_data_sheet = mydata$sampledata, 
      feature_runmode_col = VAR_feat_anno_run_mode_col, 
      batch_ids = NULL  
    )
    
    mydata$raw_metabolitedata <- mydata$metabolitedata
    mydata$metabolitedata <- norm_metabolite_data
    rm(norm_metabolite_data)
    cat( paste0("        - Normalization completed.\n\n") )
  }
  
  return(mydata)
}
