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



# Signal peak detection
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



# Create peak_df
create_peak_df <- function(fluo4_data_test){
  n_cells = nrow(fluo4_data_test)
  peak_df = data.frame(well = NULL, cell_id = NULL, peak_id = NULL, left = NULL, max = NULL, right = NULL)
  for(i in 1:n_cells){
    well = unlist(fluo4_data_test$Well[i])
    cell_id = unlist( fluo4_data_test$Cell_ID[i] )
    y = unlist( fluo4_data_test[i, -(1:2)] )
    peaks = detect_peak(y = y)
    n_peaks = nrow(peaks)
    peak_df = rbind(peak_df, data.frame(well = rep(well, n_peaks), cell_id = rep(cell_id, n_peaks), peak_id = paste0("peak", 1:n_peaks), peaks))
  }
  return(peak_df)
}



# Single-peak cells
single_peak_signal <- function(peak_df){
  unique_cell_count <- peak_df %>% group_by(cell_id) %>% count(cell_id)
  single_peak_cells <- unique_cell_count[ which(unique_cell_count$n == 1) , 1]
  return(single_peak_cells)
}



# Single-peak signal elimination 
rid_single_peak_signal <- function(peak_df, single_peak_cells){
  peak_df <- dplyr::filter(peak_df, !(cell_id %in% unlist(single_peak_cells$cell_id)))
  return(peak_df)
}
rid_single_peak_cell <- function(fluo4_data_test, single_peak_cells){
  fluo4_data_test <- dplyr::filter(fluo4_data_test, !(Cell_ID %in% unlist(single_peak_cells$cell_id)))
  return( fluo4_data_test)
}



# Add phase to peak_df
add_phase <- function(peak_df){
  Peak_distance = NULL
  Peak_distance_median = NULL
  peak_distance = NULL
  peak_distance_median = NULL
  # select 3 columns from peak_amp_df to make dataframe 'phase' 
  phase = dplyr::select(peak_df, cell_id, peak_id, max)
  # run for-loop to obtain column 'time' per unique cell
  for (cellid in unique(phase$cell_id)) {
    peak_distance = NULL
    phase_n = filter(phase, cell_id == cellid) # filter out a unique cell from dataframe 'phase' 
    peak_n = sum(phase_n$cell_id == cellid) # number of peaks from a cell 
    if(peak_n == 1){
      peak_distance = 0
    }else{
      for (i in 1:peak_n) {
        if (i == 1) {
          peak_distance[i] <- phase_n$max[i+1] - phase_n$max[i]
        } else {
          peak_distance[i] <- phase_n$max[i] - phase_n$max[i-1]
        }
      }
    }
    peak_distance_median = median(unlist(peak_distance), na.rm=TRUE)
    peak_distance_median = rep(peak_distance_median, peak_n)
    Peak_distance_median = c(Peak_distance_median, peak_distance_median)
    Peak_distance = c(Peak_distance, peak_distance)
  }
  phase = cbind(phase, Peak_distance, Peak_distance_median)
  peak_df = cbind(peak_df, dplyr::select(phase, Peak_distance, Peak_distance_median))
  return(peak_df)
}



# Analyze peak abnormality in each signal 
analyze_peak_abnormality <- function(peak_df){
  N_peaks = nrow(peak_df)
  peak_df$peak_status_analytical_algorithm = "normal" # not yet
  peak_var_df = NULL
  for (cell_id in peak_df$cell_id) {
    n_peak = sum(peak_df$cell_id == cell_id)
    #i = which(peak_df$cell_id == cell_id)
    if(n_peak > 0){
      y = unlist( fluo4_data_test[fluo4_data_test$Cell_ID == cell_id, -(1:2)] )
      sub_peak_df = peak_df[peak_df$cell_id == cell_id, ]
      peak_vars = get_peak_var(y=y, peaks=sub_peak_df[, c("left", "max", "right")])
      peak_vars$cell_id = cell_id
      peak_vars$peak_id = paste0("peak", 1:n_peak)
      peak_var_df = rbind(peak_var_df, peak_vars) # save peak variables
      A_max = apply(peak_vars[, c("A_l", "A_r")], 1, max)
      A_min = apply(peak_vars[, c("A_l", "A_r")], 1, min)
      A_avg = mean(A_max)
      ## Irregular Phase assessment 
      # consider peaks that have not exhibited any other anomalies, except the peaks that exhibit double peak anomaly which are treated as a single peak that has the mean of the positions of double peaks as its reference value. Treating a double peak as a single peak prevents the overlap of irregularity and double peak anomaly. 
      # distance of the peaks differs by a user-defined percentage (e.g. 90%) from the median of the peak distances. 
      ud = 0.9 # user-defined percentage 
      for (k in 1:n_peak) {
        if ( (sub_peak_df$Peak_distance[k] < (1 - ud) * sub_peak_df$Peak_distance_median[k]) | (sub_peak_df$Peak_distance[k] > (1 + ud) * sub_peak_df$Peak_distance_median[k]) ) {
          sub_peak_df$peak_status_analytical_algorithm[k] = "abnormal"
        }
      }
      ## Assess peaks by amplitude
      if(A_max[1] < 0.5 * A_avg){
        sub_peak_df$peak_status_analytical_algorithm[1] = "abnormal"
      }
      for(k in 2:n_peak){
        if( (sub_peak_df$peak_status_analytical_algorithm[k-1] == "abnormal") & (A_max[k] < 0.5 * A_avg) ){
          sub_peak_df$peak_status_analytical_algorithm[k] = "abnormal"
        } else if( (sub_peak_df$peak_status_analytical_algorithm[k-1] == "normal") & (A_max[k] < 0.5 * A_max[k-1]) ){
          sub_peak_df$peak_status_analytical_algorithm[k] = "abnormal"
        }
      }
      ## Assess peaks by asymmetry
      for(k in 1:n_peak){
        if( (sub_peak_df$peak_status_analytical_algorithm[k] == "normal") & (A_min[k] < 0.80 * A_max[k]) ) {
          sub_peak_df$peak_status_analytical_algorithm[k] = "abnormal"
        }
      }
      # update peak status
      peak_df$peak_status_analytical_algorithm[peak_df$cell_id == cell_id] = sub_peak_df$peak_status_analytical_algorithm
    }
  }
  peak_info_df = merge(peak_df, peak_var_df, by = c("cell_id", "peak_id"), sort = FALSE)
  peak_info_df <- peak_info_df %>% mutate(cell_peak = paste(peak_info_df$cell_id, peak_info_df$peak_id, sep = ":"))
  peak_info_df = peak_info_df %>% distinct(cell_peak, .keep_all = TRUE)
  return(peak_info_df)
}



# Analyze cell abnormality 
analyze_cell_abnormality <- function(fluo4_data_test){
  cell_status_df = fluo4_data_test[, c(2,3)]
  names(cell_status_df)[1] <- paste("cell_id")
  names(cell_status_df)[2] <- paste("cell_status")
  cell_status_df[,"cell_status"] <- NA
  cell_status_df$cell_status = "normal"
  for(i in 1:nrow(cell_status_df)){
    cell_id = cell_status_df$cell_id[i]
    peak_status_analytical_algorithm = peak_df$peak_status_analytical_algorithm[peak_df$cell_id == cell_id]
    if(length(peak_status_analytical_algorithm) == 0){
      cell_status_df$cell_status[i] = "abnormal"
    }
    else if("abnormal" %in% peak_status_analytical_algorithm){
      cell_status_df$cell_status[i] = "abnormal"
    }
  }
  cell_status_df <- filter(cell_status_df, !(cell_id %in% c(single_peak_cells$cell_id)))
  return(cell_status_df)
}



# Update peak status
update_peak_status <- function(peak_info_df){
  peak_info_df$cell_peak_id = paste(peak_info_df$cell_id, peak_info_df$peak_id, sep = ":")
  peak_info_df = peak_info_df %>% distinct(cell_peak_id, .keep_all = TRUE)
  return(peak_info_df)
}



# Create test data with cell status from analytical method
create_test_data_with_status <- function(fluo4_data_test, cell_status_df){
  names(fluo4_data_test)[2] <- paste("cell_id")
  data_with_status = merge(cell_status_df, fluo4_data_test, by = c("cell_id"))
  data_with_status = data_with_status %>% distinct(cell_id, .keep_all = TRUE)
  return(data_with_status)
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



get_delta.var <- function(delta) {
  if(length(delta)>2) 
  {
    return(var(delta[-1]))
  }
  else{
    return(0)
  }
}



# peak svm
peak_data_svm <- function(svm_fit_peak, peak_df){
  # standardization of test peak-level variables
  peak_df[, eval(var_list)] <- scale(peak_df[, eval(var_list)])
  ## novel peak status prediction using svm_fit_peak
  peak_svm = as.vector(predict(svm_fit_peak, newdata = peak_df, decision.values = TRUE, probability = TRUE))
  peak_df$peak_svm <- peak_svm
  ## Peak prediction based on AND /  OR of ML/AnalyticalMethod
  peak_status_analytical_algorithm_or_svm <- rep("normal", nrow(peak_df))
  peak_status_analytical_algorithm_or_svm[peak_svm == "abnormal" |  peak_df$peak_status_analytical_algorithm == "abnormal"] = "abnormal"
  # Add predictions to Train Data
  peak_df <- peak_df %>% mutate(peak_svm = peak_svm, 
                                peak_status_analytical_algorithm_or_svm = factor(peak_status_analytical_algorithm_or_svm))
  return(peak_df)
}



# cell status prediction
cell_status_prediction <- function(svm_fit_cell, cell_var_list, peak_df){
  ## Based on MLpredPeak or PeakIdentify, produce cell_status_analytical_algorithm_or_svm
  test_pred_abnormal_cell <- peak_df %>% dplyr::filter(peak_svm == "abnormal") %>% dplyr::select(cell_id) %>% unique() %>% unlist() %>% as.character()
  test_pred_normal_cell <- peak_df %>% dplyr::filter(! cell_id %in% test_pred_abnormal_cell)  %>% dplyr::select(cell_id) %>% unique() %>% unlist() %>% as.character()
  test_cell_pred <- data.frame(cell_id = c(test_pred_abnormal_cell, test_pred_normal_cell), 
                               peak_status_analytical_algorithm_or_svm = c(rep("abnormal", length(test_pred_abnormal_cell)), rep("normal",length(test_pred_normal_cell))))
  Cell_var <- peak_df %>% group_by(cell_id) %>% summarise(n_peak = n(), 
                                                          n_abnormal = sum(peak_status_analytical_algorithm_or_svm == "abnormal"),
                                                          prop_abnormal = n_abnormal/n_peak,
                                                          var_A = var( (A_l + A_r)/2 ), 
                                                          var_delta = get_delta.var(delta), 
                                                          var_R = var(R))
  Cell_var = merge(test_cell_pred, Cell_var, by = c("cell_id"))
  Cell_var[, eval(cell_var_list)] <- scale(Cell_var[, eval(cell_var_list)])
  ## Cell status prediction by SVM
  pred_cell = predict(svm_fit_cell, newdata=Cell_var, decision.values = TRUE, probability = TRUE)
  fluo4_data_test$cell_status = as.factor(pred_cell)
  return(fluo4_data_test)
}


