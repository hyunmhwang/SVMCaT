#' First derivative
#'
#' @export
#' @keywords
#' @examples
#' finite.differences(x, y)


# Function to obtain first dirivative
finite.differences <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }
  n <- length(x)
  fdx <- vector(length = n)
  for (i in 2:n) {
    fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])
  }
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  return(fdx)
}
