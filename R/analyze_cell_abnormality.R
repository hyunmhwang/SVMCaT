#' Analyze cell abnormality
#'
#' @export
#' @keywords
#' @examples
#' analyze_cell_abnormality()

analyze_cell_abnormality <- function(fluo4_data_test){
  cell_status_df = fluo4_data_test[, c(2,3)]
  names(cell_status_df)[1] <- paste("cell_id")
  names(cell_status_df)[2] <- paste("cell_status")
  cell_status_df[,"cell_status"] <- NA
  cell_status_df$cell_status = "normal"
  for(i in 1:nrow(cell_status_df)){
    cell_id = cell_status_df$cell_id[i]
    peak_status_analytical_algorithm = peak_df$peak_status_analytical_algorithm[peak_df$cell_id == cell_id]
    if(length(peak_status_analytical_algorithm) == 0){
      cell_status_df$cell_status[i] = "abnormal"
    }
    else if("abnormal" %in% peak_status_analytical_algorithm){
      cell_status_df$cell_status[i] = "abnormal"
    }
  }
  cell_status_df <- filter(cell_status_df, !(cell_id %in% c(single_peak_cells$cell_id)))
  return(cell_status_df)
}
