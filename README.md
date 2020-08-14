# SVMCaT
This R library is developed for analyzing Ca2+ signals with functions for 
- Identifying signal peaks
- Generating peak-level variables
- Training peak-level SVM model based on peak-level variables to differentiate normal and abnormal peaks
- Generating cell-level variables
- Training cell-level SVM models based on cell-level variables to differentiate normal and abnormal Ca2+ signals.

## Example Data
- Trained peak-level SVM model `svm_fit_peak.rds` and cell-level SVM model `svm_fit_cell.rds` are provided with this package
- Test Ca2+ signals are provided `test.csv` with this package
- See `/vignettes/demo.Rmd` for a script of demostration for using this tool

## License
[MIT](https://choosealicense.com/licenses/mit/)
