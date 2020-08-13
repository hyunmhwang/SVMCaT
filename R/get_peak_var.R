#' Function to obtain peak variables
#'
#' @export
#' @keywords
#' @examples
#' get_peak_var()

get_peak_var <- function(x=1:60, y, peaks){
  names(y) = NULL
  n_peak = nrow(peaks)
  dy = finite.differences(x, y) # derivative of y
  d2y = finite.differences(x, dy) # second derivative of y
  peak_var = data.frame(A_l = NULL, A_r = NULL, A_d = NULL, D_l = NULL, D_r = NULL,
                        Dy_max = NULL, Dy_min = NULL, D2y_max = NULL,
                        D2y_min = NULL,
                        R = NULL, delta = NULL,
                        delta_l2Dymax = NULL, delta_m2Dymin = NULL)
  for(k in 1:n_peak){
    x_l = peaks$left[k]; x_m = peaks$max[k]; x_r = peaks$right[k];
    A_l = y[x_m] - y[x_l]; # left amplitutde
    A_r = y[x_m] - y[x_r]; # right amplitude
    A_d = A_r - A_l; # difference between right and left amplitude
    D_l = x_m - x_l; # left duration
    D_r = x_r - x_m; # right duration
    Dy_max = max(dy[x_l:x_m]); # max of left side first derivative
    x_l_Dy = (x_l:x_m)[which.max(dy[x_l:x_m])];
    delta_l2Dymax = x_l_Dy - x_l; # duration from peak beginning to Dy_max
    Dy_min = abs(min(dy[x_m:x_r])) # absolute min of right side first derivative
    x_r_Dy = (x_m:x_r)[which.min(dy[x_m:x_r])]
    delta_m2Dymin = x_r_Dy - x_m; # duration from peak maximum to Dy_min
    D2y_max = max(d2y[x_l:x_m]); # maximum of right side second derivative
    D2y_min = abs(min(d2y[x_m:x_r])) # absolute minimum of right side second derivative
    R = get_peak_area(x[x_l:x_r], y[x_l:x_r]) # peak area
    if(k == 1){
      delta = x_m - 1
    }else{
      delta = x_m - peaks$max[k-1]
    }
    peak_var = rbind(peak_var,
                     data.frame(
                       A_l = A_l, A_r = A_r, A_d = A_d,
                       D_l = D_l, D_r = D_r, Dy_max = Dy_max,
                       Dy_min = Dy_min, D2y_max = D2y_max, D2y_min = D2y_min,
                       R = R, delta = delta,
                       delta_l2Dymax = delta_l2Dymax,
                       delta_m2Dymin = delta_m2Dymin))
  }
  return(peak_var)
}
