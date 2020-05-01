# calcium_transient_ml_project

Name: 
Identification of abnormal Ca2+ transient data from human induced pluripotent stem cell-derived cardiomyocytes by machine learning method

Description:
Human-induced pluripotent stem cell-derived cardiomyocytes (hiPSC-CMs) provide an excellent platform for potential clinical and research applications. Ca2+ transient analysis is a crucial method to evaluate cardiomyocyte function and has the potential to serve as a high-throughput assay. However, manual identification of abnormal Ca2+ transients is a labor-intensive and time-consuming process. By adapting the analytic algorithm proposed by Juhola et al. for analyzing Ca2+ transients, we first automatically identified Ca2+ transient peaks and determined the abnormality of Ca2+ transient peaks based on their amplitudes, symmetry properties, and peak distances, and profiled 15 peak specific quantitative variables to characterize identified peaks. Second, we implemented the support vector machine (SVM) method to learn more information from human expert assessment of peak normality and profiled peak variables in training data, and then used the trained SVM model to identify abnormal peaks in the test data. By implementing the SVM method with leave-one-out cross validation (LOOCV), we obtained SVM assessment about the normality of Ca2+ transient peaks and signals by taking these containing at least one abnormal peak as an abnormal signal. Third, we profiled the following four signal specific variables including peak assessment by the analytical method and cell assessment by expert. Lastly, we implemented SVM method using training data with signal specific variables, normality assessment based on SVM assessment of peak normality, and normality assessment by experts, to train a machine learning model that can be used to assess test Ca2+ transient signals. We tested the accuracy of identifying abnormal Ca2+ transient signals by applying our machine learning method with LOOCV to 201 cells and an independent data with 92 cells. We obtained 88% accuracy, 96% sensitivity, and 79% specificity for the LOOCV data, and 87% accuracy, 89% sensitivity, and 83% specificity for the independent data. Further, we aim to develop an R package for analyzing Ca2+ transient signals that can implement our proposed analysis procedure to detect peaks, profile peak specific variables, generate SVM assessment for peak normality, profile signal specific variables, and then train an SVM model for assessing signal normality.  Our R package can be employed in Ca2+ analysis to allow throughput analysis for drug discovery and basic disease mechanism study purposes.


Badges

Visuals

Installation: 
Simply running our script in RStudio should do. 

Usage

Support

Roadmap

Contributing

Authors and acknowledgment:
J.M. performed the Ca2+ transient assay. R.L. performed expert peak and cell assessments. H.H. and J.Y. performed the analysis. R.L. wrote the methods section of the manuscript. H.H. wrote the manuscript in consultation with R.L., J.Y., and C.X. 

License
