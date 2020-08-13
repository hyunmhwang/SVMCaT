#' get cell-well amp_l, amp_m, and amp_r
#'
#' @export
#' @keywords
#' @examples
#' get_peak_amp()

get_peak_amp <- function(x, y, peaks){
  names(y) = NULL
  n_peak = nrow(peaks)
  peak_amp = data.frame(cell_ID = NULL, peak_ID = NULL,
                        x_l = NULL, x_m = NULL, x_r = NULL,
                        y_l = NULL, y_max = NULL, y_r = NULL, y_min = NULL,
                        s_l = NULL, s_r = NULL, s_h = NULL)
  for(k in 1:n_peak){
    x_l = peaks$left[k];
    x_m = peaks$max[k];
    x_r = peaks$right[k];
    y_l = y[x_l]; # left amplitutde
    y_m = y[x_m]; # max amplitude
    y_r = y[x_r]; # right amplitude
    y_min = min(y_l, y_r) # section minimum
    s_l = x_m - x_l; # section left
    s_r = x_r - x_m; # section right
    s_h = y_m - (y_l + y_r) / 2; # section height
    peak_amp = rbind(peak_amp,
                     data.frame(x_l = x_l, x_m = x_m, x_r = x_r,
                                y_l = y_l, y_m = y_m, y_r = y_r, y_min = y_min,
                                s_l = s_l, s_r = s_r, s_h = s_h))
  }
  return(peak_amp)
}
