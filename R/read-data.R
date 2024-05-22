#' Read and Process Metabolite Data
#'
#' This function reads and processes metabolite data from a file, performs various formatting and cleaning steps,
#' and returns the processed data as a data frame. The function supports CSV, TXT, and TSV file formats.
#'
#' @param VAR_data A character string specifying the path to the metabolite data file.
#' @param VAR_isflat A numeric value indicating whether the data is in a flat file format (1 for flat, 0 otherwise).
#' @param VAR_platform A character string specifying the platform. If "Nightingale", specific name formatting is applied.
#' @return A data frame containing the processed metabolite data.
#' @examples
#' \dontrun{
#'   processed_data <- read_data("path/to/metabolite_data.csv", 1, "Nightingale")
#' }
#' @export
read_data <- function(VAR_data, VAR_isflat, VAR_platform) {
  cat(paste0("II. Data \n"))
  
  if (VAR_isflat > 0) {
    if (length(grep(".csv", VAR_data)) > 0) {
      cat(paste0("\t- Reading in your csv metabolite file\n"))
      metabolitedata <- read.csv(VAR_data, header = TRUE, as.is = TRUE, na.strings = c("NA", "NDEF", "TAG", -1, -9), row.names = 1)
    } else if (length(c(grep(".txt", VAR_data), grep(".tsv", VAR_data))) > 0) {
      cat(paste0("\t- Reading in your txt|tsv metabolite file\n"))
      metabolitedata <- read.table(VAR_data, header = TRUE, as.is = TRUE, sep = "\t", na.strings = c("NA", "NDEF", "TAG", -1, -9), row.names = 1)
    }
    
    # Remove any possible commas
    metabolitedata <- apply(metabolitedata, 2, function(x) {
      sapply(x, function(y) {
        gsub(",", "", y)
      })
    })
    metabolitedata <- as.data.frame(metabolitedata)
    
    # Format metabolite data:
    # Check if row names are just numerics. If yes, redefine rownames and column 1 values
    editrownames <- sum(rownames(metabolitedata) == 1:nrow(metabolitedata)) / nrow(metabolitedata)
    if (editrownames == 1) {
      cat(paste0("\t\t- Assuming sample IDs are in column 1 and redefining rownames\n"))
      rownames(metabolitedata) <- as.character(metabolitedata[, 1])
      metabolitedata <- metabolitedata[, -1]
    }
    
    # Check if column names are just numerics. R will add X to numeric column names
    editcolnames <- sum(substring(colnames(metabolitedata), 1, 1) == "X", na.rm = TRUE) / ncol(metabolitedata)
    if (editcolnames == 1) {
      cat(paste0("\t\t- Column names are numeric. Adding 'featID_' prefix to each.\n"))
      numeric_ids <- substring(colnames(metabolitedata), 2, nchar(colnames(metabolitedata)))
      new_col_id <- paste0("featID_", as.character(numeric_ids))
      colnames(metabolitedata) <- new_col_id
    }
    
    # If VAR_platform is Nightingale, edit metabolite names
    if (VAR_platform == "Nightingale") {
      cat(paste0("\t\t- Your defined platform is Nightingale, so editing metabolite names in an attempt to match the metaboprep annotation file.\n"))
      colnames(metabolitedata) <- gsub("_.", "pct", colnames(metabolitedata))
      colnames(metabolitedata) <- gsub("%", "pct", colnames(metabolitedata))
      colnames(metabolitedata) <- gsub("by", "", colnames(metabolitedata))
      colnames(metabolitedata) <- gsub("/", "_", colnames(metabolitedata))
      colnames(metabolitedata) <- gsub("\\.", "", colnames(metabolitedata))
      colnames(metabolitedata) <- gsub("-", "", colnames(metabolitedata))
      colnames(metabolitedata) <- gsub("_", "", colnames(metabolitedata))
      colnames(metabolitedata) <- tolower(colnames(metabolitedata))
    }
    
    # Ensure everything is numeric
    ids <- rownames(metabolitedata)
    metabolitedata <- apply(metabolitedata, 2, function(x) {
      as.numeric(as.character(x))
    })
    rownames(metabolitedata) <- ids
  }
  
  return(metabolitedata)
}

#' Read and Process Feature Annotation File
#'
#' This function reads and processes a feature annotation file, performs various formatting and cleaning steps,
#' and returns the processed data as a data frame. The function supports CSV, TXT, and TSV file formats.
#'
#' @param VAR_data_dir A character string specifying the directory path where the feature annotation file is located.
#' @param VAR_FeatureAnno_file2process A character string specifying the name of the feature annotation file.
#' @param VAR_platform A character string specifying the platform. If "Nightingale", specific name formatting is applied.
#' @param metabolitedata A data frame containing the metabolite data.
#' @return A data frame containing the processed feature annotation data.
#' @examples
#' \dontrun{
#'   processed_feature_data <- process_feature_annotation("path/to/data_dir", "feature_annotation.csv", "Nightingale", metabolitedata)
#' }
#' @export
feature_annotation <- function(VAR_data_dir, VAR_FeatureAnno_file2process, VAR_platform, metabolitedata) {
  if (!is.na(VAR_FeatureAnno_file2process)) {
    # Full path to feature annotation file
    n <- paste0(VAR_data_dir, VAR_FeatureAnno_file2process)
    
    if (length(grep(".csv", n)) > 0) {
      cat(paste0("\t- Reading in your csv feature annotation file\n"))
      featuredata <- read.csv(n, header = TRUE, as.is = TRUE, quote = "", fill = TRUE)
    } else if (length(c(grep(".txt", n), grep(".tsv", n))) > 0) {
      cat(paste0("\t- Reading in your txt feature annotation file\n"))
      featuredata <- readr::read_delim(n, show_col_types = FALSE)
      featuredata <- as.data.frame(featuredata)
    }
    
    # Format featuredata
    editrownames <- sum(rownames(featuredata) == 1:nrow(featuredata)) / nrow(featuredata)
    if (editrownames == 1) {
      # Redefine only if the number of unique strings is the same as the number of rows
      if (length(unique(featuredata[, 1])) == nrow(featuredata)) {
        cat(paste0("\t\t- Assuming feature IDs are in column 1 and redefining rownames\n"))
        # Check if column names of metabolite data file were numeric
        editcolnames <- sum(substring(colnames(metabolitedata), 1, 1) == "X", na.rm = TRUE) / ncol(metabolitedata)
        if (editcolnames == 1) {
          cat(paste0("\t\t- Column were numeric so will also add 'featID_' as a prefix to each row name here.\n"))
          rownames(featuredata) <- paste0("featID_", as.character(featuredata[, 1]))
        } else {
          rownames(featuredata) <- as.character(featuredata[, 1])
        }
      }
    }
    
    # Make sure there is a "feature_names" column and that it has the same values as rownames
    featuredata$feature_names <- rownames(featuredata)
    
    # If VAR_platform is Nightingale, edit metabolite names
    if (VAR_platform == "Nightingale") {
      cat(paste0("\t\t- Your defined VAR_platform is Nightingale,\n\t\t  so editing metabolite names in an attempt to match the metaboprep annotation file.\n"))
      rownames(featuredata) <- gsub("_.", "pct", rownames(featuredata))
      rownames(featuredata) <- gsub("%", "pct", rownames(featuredata))
      rownames(featuredata) <- gsub("/", "_", rownames(featuredata))
      rownames(featuredata) <- gsub("\\.", "", rownames(featuredata))
      rownames(featuredata) <- gsub("-", "", rownames(featuredata))
      rownames(featuredata) <- gsub("_", "", rownames(featuredata))
      rownames(featuredata) <- tolower(rownames(featuredata))
    }
    
    # Verify that the first column of feature names matches the column names in metabolite data
    idsmatch <- sum(rownames(featuredata) %in% colnames(metabolitedata)) == nrow(featuredata)
    if (!idsmatch) {
      stop(paste0("The feature or metabolite IDs in the first column of your feature annotation file
        \tdo not match the feature (column) ids in the metabolite data file.
        \tPlease ensure that these ids match and come back and try again."),
           call. = FALSE)
    }
    
    # If the IDs do match, make sure they are ordered properly
    if (idsmatch) {
      m <- match(colnames(metabolitedata), rownames(featuredata))
      featuredata <- featuredata[m, ]
    }
  }
  
  return(featuredata)
}

#' Read and Process Sample Annotation File
#'
#' This function reads and processes a sample annotation file, performs various formatting and cleaning steps,
#' and returns the processed data as a data frame. The function supports CSV, TXT, and TSV file formats.
#'
#' @param VAR_data_dir A character string specifying the directory path where the sample annotation file is located.
#' @param VAR_SampleAnno_file2process A character string specifying the name of the sample annotation file.
#' @param metabolitedata A data frame containing the metabolite data.
#' @return A data frame containing the processed sample annotation data. Returns NULL if VAR_SampleAnno_file2process is NA.
#' @examples
#' \dontrun{
#'   processed_sample_data <- process_sample_annotation("path/to/data_dir", "sample_annotation.csv", metabolitedata)
#' }
#' @export
sample_annotation <- function(VAR_data_dir, VAR_SampleAnno_file2process, metabolitedata) {
  if (!is.na(VAR_SampleAnno_file2process)) {
    # Full path to sample annotation file
    n <- paste0(VAR_data_dir, VAR_SampleAnno_file2process)
    
    if (length(grep(".csv", n)) > 0) {
      cat(paste0("\t- Reading in your csv sample annotation file\n"))
      sampledata <- read.csv(n, header = TRUE, quote = "", as.is = TRUE)
    } else if (length(c(grep(".txt", n), grep(".tsv", n))) > 0) {
      cat(paste0("\t- Reading in your txt sample annotation file\n"))
      sampledata <- read.table(n, header = TRUE, as.is = TRUE, sep = "\t")
    }
    
    # Format sampledata
    editrownames <- sum(rownames(sampledata) == 1:nrow(sampledata)) / nrow(sampledata)
    if (editrownames == 1) {
      cat(paste0("\t\t- Assuming sample IDs are in column 1 and redefining rownames\n"))
      rownames(sampledata) <- as.character(sampledata[, 1])
    }
    
    # Verify that the first column of sample IDs matches the rownames in metabolite data
    idsmatch <- sum(rownames(sampledata) %in% rownames(metabolitedata)) == nrow(sampledata)
    if (!idsmatch) {
      stop(paste0("The sample IDs in the first column of your sample annotation file
        \tdo not match the sample (row) IDs in the metabolite data file.
        \tPlease ensure that these IDs match and come back and try again."),
           call. = FALSE)
    }
    
    # If the IDs do match, make sure they are ordered properly
    if (idsmatch) {
      m <- match(rownames(metabolitedata), rownames(sampledata))
      sampledata <- sampledata[m, ]
    }
  }
  
  return(sampledata)
}

#' Generate Sample and Feature Data
#'
#' This function generates sample and feature data based on the metabolite data and platform.
#'
#' @param VAR_isflat A numeric value indicating whether the data is in a flat file format (1 for flat, 0 otherwise).
#' @param metabolitedata A data frame containing the metabolite data.
#' @param ng_anno A data frame containing the Nightingale annotation data.
#' @param VAR_platform A character string specifying the platform. If "Nightingale", specific name formatting is applied.
#' @return A list containing the sample and feature data.
#' @examples
#' \dontrun{
#'   sample_and_feature_data <- generate_sample_and_feature_data(VAR_isflat = 1, metabolitedata = my_metabolite_data, ng_anno = my_annotation_data, VAR_platform = "Nightingale")
#' }
generate_sample_and_feature_data <- function(VAR_isflat, metabolitedata, ng_anno, VAR_platform) {
  
  if (VAR_isflat > 0) {
    if (!exists("sampledata")) {
      sampledata <- data.frame(SampleID = rownames(metabolitedata))
    }
    
    if (!exists("featuredata")) {
      featuredata <- data.frame(feature_names = colnames(metabolitedata))
      rownames(featuredata) <- as.character(featuredata[, 1])
    }
    
    if (VAR_platform == "Nightingale") {
      m <- match(rownames(featuredata), ng_anno$metabolite)
      featuredata <- cbind(featuredata, ng_anno[m, -1])
      
      w <- which(is.na(m))
      if (length(w) > 0) {
        ids_i_could_not_match <- rownames(featuredata)[w]
        featuredata[w, c("feature_names", "raw.label", "class", "subclass", "label", "label.no.units", "derived_features")] <- "unknown"
        cat(paste0("\n\n
                   \t- !!!! WARNING !!!!
                   \t\tPlease be aware that there is|are ",  length(w) ," metabolite IDs in your Nightingale data set that metaboprep could not annotate.
                   \t\tAt present Nightingale Health data annotation is performed with the metaboprep data frame ng_anno. 
                   \t\tYou can see it by opening an R session and typing
                   \n\t\t\t> metaboprep::ng_anno
                   \n\t\t\tThe annotation is necessary for proper treatment of derived features in the Nightingale data sets 
                   \t\t(i.e. ratios or features derived from two or more features already present in the data).
                   \t\tThis error may be the product of a new spelling Nightingale Health has chosen or the addition of a new feature.
                   \t\tThis warning is ONLY an issue IF the new or unmapped metabolite is a derived feature and you are excluding derived features
                   \t\tfrom the independent feature identification step, that also feeds into the PCA.
                   \t\tIf it is spelling mismatch you can edit your files and try again.
                   \t\tIf it is a new feature, we are aware of this issue, and are working on a solution.\n"))
        cat(paste0("\n\t\t- The features metaboprep could not match are:\n"))
        cat(paste0(ids_i_could_not_match, "\n\n\n------\n"))
      }
    }
    
    mydata <- list(metabolitedata = metabolitedata, sampledata = sampledata, featuredata = featuredata)
  }
  
  return(mydata)
}

#' Process Excel Data
#'
#' This function processes Excel data based on the platform.
#'
#' @param VAR_isexcel A numeric value indicating whether the data is in Excel format (1 for Excel, 0 otherwise).
#' @param VAR_platform A character string specifying the platform. If "Nightingale" or "Metabolon", specific processing is applied.
#' @param VAR_METABO_file2process A character string specifying the name of the Excel file to process.
#' @param VAR_data_dir A character string specifying the directory path where the data file is located.
#' @param VAR_project A character string specifying the project name.
#' @return A list containing the processed data.
#' @examples
#' \dontrun{
#'   processed_data <- process_excel_data(VAR_isexcel = 1, VAR_platform = "Nightingale", VAR_METABO_file2process = "data.xlsx", VAR_data_dir = "data/", VAR_project = "Project")
#' }
generate_sample_and_feature_data_excel <- function(VAR_isexcel, VAR_platform, VAR_METABO_file2process, VAR_data_dir, VAR_project) {
  mydata <- NULL
  
  if (VAR_isexcel > 0) {
    if (VAR_platform == "Nightingale") {
      cat(paste0("II. Processing your Nightingale data.\n"))
      if (!is.na(VAR_METABO_file2process[3, 2])) {
        mydata <- read.in.nightingale(file2process = VAR_METABO_file2process, VAR_data_dir = VAR_data_dir, projectname = VAR_project)
        cat(paste0("\t- Your raw Nightingale data has been read in and converted to working tab delimited text files.\n\n"))
      } 
    }
    
    if (VAR_platform == "Metabolon") {
      cat(paste0("II. Processing your Metabolon data.\n"))
      if (!is.na(VAR_METABO_file2process[3, 2])) {
        mydata <- read.in.metabolon(file2process = VAR_METABO_file2process, VAR_data_dir = VAR_data_dir, projectname = VAR_project)
        cat(paste0("\t- Your raw Metabolon data has been read in and converted to working tab delimited text files.\n\n"))
      } 
    }
  }
  
  return(mydata)
}
