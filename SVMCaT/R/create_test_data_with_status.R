#' Create test data with cell status from analytical method
#'
#' @export
#' @keywords
#' @examples
#' create_test_data_with_status()

create_test_data_with_status <- function(fluo4_data_test, cell_status_df){
  names(fluo4_data_test)[2] <- paste("cell_id")
  data_with_status = merge(cell_status_df, fluo4_data_test, by = c("cell_id"))
  data_with_status = data_with_status %>% distinct(cell_id, .keep_all = TRUE)
  return(data_with_status)
}
