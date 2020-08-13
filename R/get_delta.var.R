#' get delta.var
#'
#' @export
#' @keywords
#' @examples
#' get_delta.var()

get_delta.var <- function(delta) {
  if(length(delta)>2)
  {
    return(var(delta[-1]))
  }
  else{
    return(0)
  }
}
