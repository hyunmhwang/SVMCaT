#' Function to calculate peak area
#'
#' @export
#' @keywords
#' @examples
#' get_peak_area()

get_peak_area <- function(x, y){
  R = 0
  p = length(x)
  for(j in 2:p){
    R = R + (y[j] + y[j-1]) * (x[j] - x[j-1]) * 0.5
  }
  R = R - (y[p] + y[1]) * (x[p] - x[1]) * 0.5
  return(R)
}
