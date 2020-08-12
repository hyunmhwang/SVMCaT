#' # Single-peak signal elimination
#'
#' @export
#' @keywords
#' @examples
#' rid_single_peak_signal()

# Single-peak signal elimination
rid_single_peak_signal <- function(peak_df, single_peak_cells){
  peak_df <- dplyr::filter(peak_df, !(cell_id %in% unlist(single_peak_cells$cell_id)))
  return(peak_df)
}
