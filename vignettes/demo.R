## -----------------------------------------------------------------------------
#library(SVMCaT)
source("/Users/hyunhwang/Desktop/Files/1/R/SVMCaT/R/Rfuncs_HH.R")

library(magrittr)
library(dplyr)
library(e1071)
library(reshape2)
library(ggplot2)

# loading svm_fit_peak and svm_fit_cell for peak and cell status prediction
svm_fit_peak <- readRDS("/Users/hyunhwang/Desktop/Files/1/R/SVMCaT/svm_fit_peak.rds")
svm_fit_cell <- readRDS("/Users/hyunhwang/Desktop/Files/1/R/SVMCaT/svm_fit_cell.rds")

# 15 peak-level variables 
var_list = c("A_l", "A_r", "A_d", "D_l", "D_r", 
             "Dy_max", "Dy_min", "D2y_max", "D2y_min", "R", 
             "delta", "delta_l2Dymax", "delta_m2Dymin", "Peak_distance", "Peak_distance_median")

# 4 cell-level variables
cell_var_list = c("n_peak", "n_abnormal", "prop_abnormal", "var_A", "var_delta", "var_R")

## -----------------------------------------------------------------------------
## Load data
fluo4_data_test <- read.csv("/Users/hyunhwang/Desktop/Files/1/R/SVMCaT/112119_ca0_1.csv")

## -----------------------------------------------------------------------------
# create peak_df
peak_df = create_peak_df(fluo4_data_test)

# single-peak signals 
single_peak_cells = single_peak_signal(peak_df)

# removal of single-peak signal
peak_df = rid_single_peak_signal(peak_df, single_peak_cells)
fluo4_data_test = rid_single_peak_cell(fluo4_data_test, single_peak_cells)

# add phase to peak_df
peak_df = add_phase(peak_df)

## -----------------------------------------------------------------------------
peak_df = analyze_peak_abnormality(peak_df)

## -----------------------------------------------------------------------------
data_with_status = analyze_cell_abnormality(fluo4_data_test)

# add data_with_status to fluo4_data_test
fluo4_data_test = create_test_data_with_status(fluo4_data_test, data_with_status)

## -----------------------------------------------------------------------------
# update peak status
peak_df = update_peak_status(peak_df)

# peak status prediction
peak_df = peak_data_svm(svm_fit_peak, peak_df)

# cell status prediction
fluo4_data_test = cell_status_prediction(svm_fit_cell, cell_var_list, peak_df)

## -----------------------------------------------------------------------------
n_cells = fluo4_data_test$cell_id # number of cells in calcium transient data
frame = 1:60 # number of frames per transient signal data

pdf("./Peak_plot_test.pdf")
for(i in n_cells){
  cell_id = i
  # print(cell_id)
  x = frame
  cells = subset(fluo4_data_test, cell_id == i)
  y = unlist(cells[, -(1:3)])
  n_peaks = sum(peak_df$cell_id == i)
  peaks = subset(peak_df, cell_id == i)
  
  p = ggplot(data = data.frame(x = x, y = y), 
    aes(x = x, y = y)) + 
    geom_line(size = 1) + 
    # geom_smooth(method = "loess", se = FALSE, span = 0.039) +
    geom_point(data = data.frame(x = c(peaks$left, peaks$right, peaks$max),
                                 y = y[c(peaks$left, peaks$right, peaks$max)], 
                                 peak = c(rep("S", 2*n_peaks), rep("M", n_peaks))),
               mapping = aes(x = x, y = y, colour = peak), size = 2) +
    annotate("text", x = peaks$max, y = y[peaks$max]*1.01, label = paste0("peak", 1:n_peaks)) +
    labs(title = cell_id, x = NULL, y = NULL)  +
    theme(text = element_text(size = 20), legend.position = "none")
  print(p)
}
dev.off()

## -----------------------------------------------------------------------------
pdf("./Signals_with_PredStatus_test.pdf")
lane_list = sort(unique(fluo4_data_test$Well))
for(lane in lane_list){
  temp = fluo4_data_test[fluo4_data_test$Well == lane, ]
  plot_data = melt(temp[, -(3)])
  plot_data$status = factor(plot_data$cell_status)
  p = ggplot(plot_data, aes(x = variable, y = value, 
                            group = lane,
                            color = status)) + 
    ylab("Intensity") +
    geom_line() + 
    facet_wrap(~cell_id, scales = "free_y") + 
    guides(color = guide_legend(title="SVMCaT Prediction")) +
    scale_x_discrete(name = "Frame", 
                   breaks = c("Frame1", "Frame20", "Frame40", "Frame60"), 
                   labels = c(1, 20, 40, 60))
  print(p)
}
dev.off()

