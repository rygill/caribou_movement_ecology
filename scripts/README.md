Scripts to conduct data prep, outlier detection, home range calculation and analysis and rsfs. Files are run in their numeric order. 

01.data_prep.r - reads in complete, raw dataset provided by the Knowledge Management Branch, cleans mortality status, removes collar locations from penned animals, creates fields to use for the moving window analysis.

02.outlier.detection.r - uses output from previous file to format for ctmm functions to assess for outliers. Assigns outlier status to points based on movement speed, collar behaviour around mortality events, and points that do not have high DOP values, but plotting and visually inspecting reveals them as outliers. 

03.ctmm_moving_window.r - uses output from previous file to iterate each over collar and calculate home range over a 14 day moving window, advanced by 1 day to determine when each animal becomes range resident. Lines 380-394 show the estimates for late-winter residency based on this calculation. These periods are used for each herd for each period. Each herd is written to file separately.

05.Batch_Run.r - 
