#' Create peak_df
#'
#' @export
#' @keywords
#' @examples
#' create_peak_df(fluo4_data_test)

create_peak_df <- function(fluo4_data_test){
  n_cells = nrow(fluo4_data_test)
  peak_df = data.frame(well = NULL, cell_id = NULL, peak_id = NULL, left = NULL, max = NULL, right = NULL)
  for(i in 1:n_cells){
    well = unlist(fluo4_data_test$Well[i])
    cell_id = unlist( fluo4_data_test$Cell_ID[i] )
    y = unlist( fluo4_data_test[i, -(1:2)] )
    peaks = detect_peak(y = y)
    n_peaks = nrow(peaks)
    peak_df = rbind(peak_df, data.frame(well = rep(well, n_peaks), cell_id = rep(cell_id, n_peaks), peak_id = paste0("peak", 1:n_peaks), peaks))
  }
  return(peak_df)
}
