#' Analyze peak abnormality in each signal
#'
#' @export
#' @keywords
#' @examples
#' analyze_peak_abnormality()

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
