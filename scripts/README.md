Scripts to conduct data prep, outlier detection, home range calculation and analysis and rsfs. Files are run in their numeric order. 
In places where data formatting requires changing period lines are denoted with: <br>
#-----------------------CHANGE PERIOD------------------------------------------#

01.data_prep.r - reads in complete, raw dataset provided by the Knowledge Management Branch (complete_data_220526.csv), cleans mortality status, removes collar locations from penned animals, creates fields to use for the moving window analysis.

02.outlier.detection.r - uses output from previous file to format for ctmm functions to assess for outliers. Assigns outlier status to points based on movement speed, collar behaviour around mortality events, and points that do not have high DOP values, but plotting and visually inspecting reveals them as outliers. 

03.ctmm_moving_window.r - uses output from previous file to iterate each over collar and calculate home range over a 14 day moving window, advanced by 1 day to determine when each animal becomes range resident. Lines 380-394 show the estimates for late-winter residency based on this calculation. These periods are used for each herd for each period. Each herd is written to file separately.

05.Batch_Run.r - relies on Fit_Mods.r to fit movement models using ctmm::ctmm.select. This script creates subdirectories and generates movement models and home range estimates. Outputs include movemement model object (rda), akde of home range as tif and shapefile and a file of home range metrics for each period for each animal.

06.extract.raster.95.hr.estimate.r - extracts the mean value of each of the landscape rasters for each home range for each period. The results for each individual for each period are appended to the results from script 05 to create a new dataset named Caribou_reesults_with_error_covariates_230214.csv

07.HR_Regression - fits a generalized linear mixed model using the lme4 package to home range estimates.

08.diffusion.regression - fits a generalized linear mixed effects model using the lme4 package to diffusion estimates.
