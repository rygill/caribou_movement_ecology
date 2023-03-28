rm(list = ls())
gc()
#library(devtools)
#install_github("ctmm-initiative/ctmm", dependencies = TRUE)

library(tidyr)
library(lubridate)
library(geosphere)
library(ctmm)

getwd()
#setwd('../ABMI/caribou_anthropause')
source("./scripts/Fit_Mods.R")

#read calibration data:
calib = read.csv('./data/input_data/calibration_220221.csv')
calib$individual.local.identifier = as.factor(calib$individual.local.identifier)
levels(calib$individual.local.identifier)
calib = droplevels(calib[!(calib$individual.local.identifier %in% c('29094', '81361')),])

calib = as.telemetry(calib)

#Load in the tracking data
#read in each herd to account for differences in timing:
cn = read.csv('./data/input_data/home_range/CN_220605_clean_data_formatted_for_tele.csv')
cs = read.csv('./data/input_data/home_range/CS_220605_clean_data_formatted_for_tele.csv')
hr = read.csv('./data/input_data/home_range/HR_220605_clean_data_formatted_for_tele.csv')

data = rbind(cn, cs, hr)

#change for each herd
#data = cn
#change for each period:
data = data[data$after == 1,]

#Store the binomial
#match this with previous line:
BINOMIAL <- "after"
  
  #Convert to a telemetry object
DATA <- ctmm::as.telemetry(data)
#
###
  
#uere(cilla) = uere.fit(calib)

#Create the paths to the file directory
#these hold the output of Fit_Mods.r
dir.create(paste("./data/home_range/", BINOMIAL, sep = ""))
dir.create(paste("./data/home_range/", BINOMIAL, "/Fits_with_error", sep = ""))
Model_path = paste("./data/home_range/", BINOMIAL, "/Fits_with_error", sep = "")
  
dir.create(paste("./data/home_range/", BINOMIAL, "/UDs_with_error", sep = ""))
UD_path = paste("./data/home_range/", BINOMIAL, "/UDs_with_error", sep = "")
  
dir.create(paste("./data/home_range/", BINOMIAL, "/Figs_with_error", sep = ""))
Fig_Path = paste("./data/home_range/", BINOMIAL, "/Figs_with_error", sep = "")
  
  
  #Then walk through each individual sequentially
for(j in 1:length(DATA)){
    #j = 1
    cat("Working on individual ", j, " of ", length(DATA), "\n")
    
    #Extract the current individual
    cilla <- DATA[[j]]
    
   #uere(cilla) <- 10  #becomes calibration field from new data
    uere(cilla) = uere.fit(calib)
    

    RESULTS <- tryCatch(
      
      {
        
        RESULTS <- CTMM_FIT(cilla,
                            Model_path = Model_path,
                            UD_path = UD_path,
                            Fig_Path = Fig_Path,
                            speed = FALSE,
                            error = TRUE)
                
        
        
      }, error=function(err) {
        message(cat("Home range estimation failed, returning NAs \n"))
        
        RESULTS <- as.data.frame(t(unlist(c(BINOMIAL,
                                            cilla@info[1],
                                            rep(NA, 6)))))
        
        return(RESULTS)
      }
    )
    
 
    if(j == 1){
      
      write.table(RESULTS,
                  file = "./data/home_range/Caribou_Results_with_Error.csv",
                  row.names=FALSE,
                  col.names=TRUE,
                  sep=",",
                  append=TRUE)
      
    } else {
      
      write.table(RESULTS,
                  file = "./data/home_range/Caribou_Results_with_Error.csv",
                  row.names=FALSE,
                  col.names=FALSE,
                  sep=",",
                  append=TRUE)
      
    }
    
    
    
    
}#Closes the loop that runs over the as.telemetry object




