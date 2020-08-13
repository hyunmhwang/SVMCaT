#' Add phase
#'
#' @export
#' @keywords
#' @examples
#' add_phase()

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
