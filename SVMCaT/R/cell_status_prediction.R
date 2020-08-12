#' cell status prediction
#'
#' @export
#' @keywords
#' @examples
#' cell_status_prediction()

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
