#' Update peak status
#'
#' @export
#' @keywords
#' @examples
#' update_peak_status()

# Update peak status
update_peak_status <- function(peak_info_df){
  peak_info_df$cell_peak_id = paste(peak_info_df$cell_id, peak_info_df$peak_id, sep = ":")
  peak_info_df = peak_info_df %>% distinct(cell_peak_id, .keep_all = TRUE)
  return(peak_info_df)
}
