

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


# Function to detect signal peaks
detect_peak <- function(x, y, tup = 30, rtup = 10){
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
  # print(peak)
  return(peak)
}


# Function to calculate peak area
get_peak_area <- function(x, y){
  R = 0
  p = length(x)
  for(j in 2:p){
    R = R + (y[j] + y[j-1]) * (x[j] - x[j-1]) * 0.5
  }
  R = R - (y[p] + y[1]) * (x[p] - x[1]) * 0.5
  return(R)
}


# Function to obtain peak variables
get_peak_var <- function(x, y, peaks){
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
    x_l = peaks$left[k]; x_m = peaks$max[k]; x_r = peaks$right[k]
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


# get cell-well amp_l, amp_m, and amp_r
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


