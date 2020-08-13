#' peak svm
#'
#' @export
#' @keywords
#' @examples
#' peak_data_svm()

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
