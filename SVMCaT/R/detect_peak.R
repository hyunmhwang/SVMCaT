#' Peak detection
#'
#' @export
#' @keywords
#' @examples
#' detect_peak(x, y, tup, rtup)

detect_peak <- function(x=1:60, y, tup=30, rtup=2){
  if (length(x) != length(y)) {
    stop('x and dy vectors must have equal length')
  }
  # remove slope by linear regression
  # y = unlist( fluo4_data_160[i, -(1:3)] )
  fit = lm(y ~ x)
  y = y - coefficients(fit)[2] * x
  dy = finite.differences(x, y) # obtain first derivative
  d2y = finite.differences(x, dy) # obtain second derivative
  n = length(x)
  peak = data.frame(left = NULL, max = NULL, right = NULL)
  peak_start_list = which(dy > tup)
  if(length(peak_start_list) == 0){
    return(NULL)
  }else if(peak_start_list[1] > 1){
    peak_l = peak_start_list[1]
  }else if(length(peak_start_list) > 1){
    peak_l = peak_start_list[2]
  }
  count = 0
  while(peak_l < n){
    # print(paste("peak left is ", peak_l))
    for(j in (peak_l+1):n){
      if(dy[j] * dy[j-1] <= 0){
        # print(paste("peak max is ", j))
        peak_m = j
        for(k in (peak_m + 1):n){
          if( dy[k]*dy[k-1] <= 0 ){
            for(k_r in (k:n) ){
              if(dy[k_r] > rtup){
                peak_r = k_r-1
                #  print(paste("peak right is ", peak_r))
                peak = rbind(peak, data.frame(left = peak_l, max = peak_m, right = peak_r))
                count = count + 1
                break
              }else if(k_r == n){
                peak_r = n
                # print(paste("peak right is ", peak_r))
                peak = rbind(peak, data.frame(left = peak_l, max = peak_m, right = peak_r))
                count = count + 1
                break
              }
            }
            if(sum(peak_start_list > peak_r) >0){
              peak_l = peak_start_list[peak_start_list > peak_r - 1][1]
            }else{
              peak_l = n
            }
            break
          }else if(k == n){
            # print(paste("peak right is ", k))
            #peak_r = k
            #peak = rbind(peak, data.frame(left = peak_l, max = peak_m, right = peak_r))
            #count = count + 1
            peak_l = n
          }
        }
        break
      }
      else if(j == n){
        # print("End of peak searching.")
        peak_l = n
      }
    }
  }
  # Function to remove noises and first peak with s < 5 and A_l < 0.5 * A_r
  # Remove noises and first peak with s < 5 and A_l < 0.5 * A_r
  if(! is.null(peak)){
    n_peaks = nrow(peak)
    # Average amplitude
    A_peaks = y[peak$max] - (y[peak$left] + y[peak$right])/2
    noise = rep(FALSE, n_peaks)
    for(k in 1:n_peaks){
      if(k == 1){
        A_l = y[peak$max[k]] - y[peak$left[k]]
        A_r = y[peak$max[k]] - y[peak$right[k]]
        if( (peak$left[k] < 5) & (A_l < 0.5 * A_r) ){
          noise[k] = TRUE
        }
      }
      if(A_peaks[k] < max(A_peaks) * 0.15 ){
        noise[k] = TRUE
      }
    }
    peak = peak[!noise, ]
    n_peaks = nrow(peak)
    return(peak)
  }
}
