#' Single-peak cells
#'
#' @export
#' @keywords
#' @examples
#' single_peak_signal()

single_peak_signal <- function(peak_df){
  unique_cell_count <- peak_df %>% group_by(cell_id) %>% count(cell_id)
  single_peak_cells <- unique_cell_count[ which(unique_cell_count$n == 1) , 1]
  return(single_peak_cells)
}
