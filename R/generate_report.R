#' generate metaboprep summary html report
#'
#' This function generates the html report.
#'
#' @param full_path_2_Rdatafile full path to the Rdatafile 
#' @param dir_4_report directory to place the report
#' @param path_2_Rmd_template full path to the html report template
#'
#' @keywords knit html report
#'
#' @return Writes a html report to file
#'
#' @importFrom knitr knit2html
#' 
#' @return an html file written to file
#'
#' @export
#'
#' @examples
#' 
#'
generate_report <- function(full_path_2_Rdatafile, dir_4_report, path_2_Rmd_template = file.path( system.file("rmarkdown", package="metaboprep"), "metaboprep_Report_v0.Rmd")) {
  # Load the R data file
  cat("Loading the R data file\n")
  load(full_path_2_Rdatafile)
  
  # Define the output file path
  output_file <- file.path(dir_4_report, paste0("metaboprep_report_", Sys.Date(), ".html"))
  
  # Ensure the directory for the report exists
  if (!dir.exists(dir_4_report)) {
    dir.create(dir_4_report, recursive = TRUE)
  }
  
  # Render the R Markdown document
  cat("Knit the report to html\n")
  rmarkdown::render(input = path_2_Rmd_template, 
                    output_file = output_file,
                    params = list(full_path_2_Rdatafile = full_path_2_Rdatafile))
  
  cat("Report generated at:", output_file, "\n")
}



