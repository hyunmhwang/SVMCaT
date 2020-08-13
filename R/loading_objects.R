#' loading objects
#'
#' @export
#' @keywords
#' @examples

# loading svm_fit_peak and svm_fit_cell for peak and cell status prediction
svm_fit_peak <- readRDS("./svm_fit_peak.rds")
svm_fit_cell <- readRDS("./svm_fit_cell.rds")



# 15 peak-level variables
var_list = c("A_l", "A_r", "A_d", "D_l", "D_r",
             "Dy_max", "Dy_min", "D2y_max", "D2y_min", "R",
             "delta", "delta_l2Dymax", "delta_m2Dymin", "Peak_distance", "Peak_distance_median")


# 4 cell-level variables
cell_var_list = c("n_peak", "n_abnormal", "prop_abnormal", "var_A", "var_delta", "var_R")
