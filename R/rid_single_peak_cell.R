#' Single-peak cell elimination
#'
#' @export
#' @keywords
#' @examples
#' rid_single_peak_cell()

rid_single_peak_cell <- function(fluo4_data_test, single_peak_cells){
  fluo4_data_test <- dplyr::filter(fluo4_data_test, !(Cell_ID %in% unlist(single_peak_cells$cell_id)))
  return( fluo4_data_test)
}
