rm(list = ls())
gc()

#detach("package:ctmm", unload = TRUE)
#library(devtools)
#remotes::install_github("ctmm-initiative/ctmm", force = TRUE)
library(ctmm)
#library(terra)
#find out where tempdir is and make sure 'raster' is in there:
path = paste0(tempdir(),'/raster')
dir.create(path)

rm(list = ls())
gc()

#Earlier each winter residency period was determined for each herd for each period,
#for this analysis we can merge these individuals because we have determined that period.
#data1 = read.csv('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/home_range/CN_220605_clean_data_formatted_for_tele.csv')
#data2 = read.csv('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/home_range/CS_220605_clean_data_formatted_for_tele.csv')
#data3 = read.csv('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/home_range/HR_220605_clean_data_formatted_for_tele.csv')
#data = rbind(data1, data2, data3)
#write.csv(data, "C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/home_range/data_220811.csv", row.names = FALSE)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#read data in:
data = read.csv("./data/input_data/data_220811.csv")
data$resident = rowSums(data[,c("prior1", "prior2", "during", "after")])
data = data[data$resident == 1,]
data = data[complete.cases(data), ]

data$timestamp = as.POSIXct(data$timestamp, format = '%Y-%m-%d %H:%M', tz='UTC')
data$herd = as.factor(data$herd)
data$mort.status = as.factor(data$mort.status)
data$outlier = as.factor(data$outlier)

#make a period column to simplify things later:
data$period = ifelse(data$prior1 == 1, 'prior1',
                     ifelse(data$prior2 == 1, 'prior2',
                            ifelse(data$during == 1, 'during', 'after')))

data$period = as.factor(data$period)

#make sure they're all live locations:
levels(data$mort.status)
#-------------------------------------------------------------------------------
#set period:
period = 'after'
data = droplevels(data[data$period == period,]) 
#-------------------------------------------------------------------------------
#there are some collars which the UD failed, remove these so the loop runs
#prior 1: 81361
#prior 2: 30613, 81365
#during: 42144
#after: 44098

#In addition to those collars above, remove those collars we deemed as having 
#problematic outliers not detected during cleaning: 32611 and 44095. Even though
#problems occured for only one period we'll remove them completely:
data = data[!(data$individual.local.identifier %in% c(32611, 44095, 44098)),]

#running out of storage causes NA values in sampled rasters error. For during, I
#restarted at 81313
#data = data[data$individual.local.identifier > 81311,]

#running out of storage causes NA values in sampled rasters error. For prior1, I
#restarted at 81331
#data = data[data$individual.local.identifier > 81331,]
#data = data[data$individual.local.identifier > 81315,]


#read in landscape rasters:
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
elev = raster::raster('./data/rasters/august/elev_220809.tif')
slope = raster::raster('./data/rasters/august/slope_220809.tif')
#roughness = rast('./data/rasters/august/roughness_220809.tif')
proj.age = raster::raster('./data/rasters/august/proj_age_220809.tif')
#NAs throwing error, reclass to 0
m <- c(NA,NA,0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
proj.age <- raster::reclassify(proj.age, rclmat, include.lowest=TRUE)

heli.ten = raster::raster('./data/rasters/august/heli_ten_220809.tif')

#named list of rasters for rsf:
rlist = list(elev, slope, heli.ten, proj.age)
names(rlist) = c('elev', 'slope', 'heli.ten', 'proj.age')

#clean up:
gc()

#Convert to a telemetry object
DATA <- ctmm::as.telemetry(data)
#DATA = DATA[19:60]
#set directory where UDs from batch_run.r reside:
#-------------------------------------------------------------------------------
#UD_path = paste0("C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/home_range/", period, "/UDs_with_error")
UD_path = paste0("./data/home_range/", period, "/UDs_with_error")
#-------------------------------------------------------------------------------

#folder with UD
setwd(UD_path)
#load UDs from that folder:
ud.dir = list.files(pattern="*.rda", all.files=TRUE, 
                    full.names=FALSE)
uds <- list()

for(i in 1:length(ud.dir)){
  load(ud.dir[[i]])
  uds[[i]] <- cilla.akde

}
#CHECK that the records in DATA match the records in the uds object. If not
#add a folder to the UD_path called censor and move UDs not in DATA into that 
#folder to avoid crashes.

#set the wd back to the original so the other code works:
setwd('../../../../')
getwd()


#create empty list to append animals to as the loop progresses
reslist = list()

#*#set path for output of rsf models
#-------------------------------------------------------------------------------
rsf_path = paste0('./data/rsf/', 
                  period)
#-------------------------------------------------------------------------------


for(j in 1:length(DATA)){
  try({
  cat("Working on individual ", j, " of ", length(DATA))
  #j = 50
  #Extract the current individual
  cilla <- DATA[[j]] 
  
  #Extract corresponding UD
  cilla.ud = uds[[j]] 
  cilla.ud@info$identity
  cilla@info$identity
  #project so they are the same
  ctmm::projection(cilla) = ctmm::projection(cilla.ud)
  
      RESULTS <- ctmm::rsf.fit(cilla,
                         UD = cilla.ud,
                         R = rlist,
                         na.rm = TRUE,
                         integrated = TRUE,
                         #integrator="Riemann",
                         level.UD = 0.95, 
                         interpolate = TRUE, 
                         #isotropic = TRUE,
                         error = TRUE)
      
      #Get basic stats on the dataset
      res <- as.data.frame(period)
      res$ID <- cilla@info$identity

#      #Get HR results
      res$elev <- summary(RESULTS, units = FALSE)$CI[4,2]
      res$elev_min <- summary(RESULTS, units = FALSE)$CI[4,1]
      res$elev_max <- summary(RESULTS, units = FALSE)$CI[4,3]
#      #.
      res$slope <- summary(RESULTS, units = FALSE)$CI[3,2]
      res$slope_min <- summary(RESULTS, units = FALSE)$CI[3,1]
      res$slope_max <- summary(RESULTS, units = FALSE)$CI[3,3]
#      #.
      res$heli.ten <- summary(RESULTS, units = FALSE)$CI[2,2]
      res$heli.ten_min <- summary(RESULTS, units = FALSE)$CI[2,1]
      res$heli.ten_max <- summary(RESULTS, units = FALSE)$CI[2,3]
#      #.
      res$proj.age <- summary(RESULTS, units = FALSE)$CI[1,2]
      res$proj.age_min <- summary(RESULTS, units = FALSE)$CI[1,1]
      res$proj.age_max <- summary(RESULTS, units = FALSE)$CI[1,3]
      
     reslist[[j]] = res 
     #save(RESULTS, file = paste0(rsf_path, cilla@info[1],'_fit.rda')) ##################################NEW
     
     #Assign the file path
     rsf.out.path <- file.path(rsf_path, paste0("rsf_",period,'_', cilla@info[1], ".rda"))
     
     
     #And save
     save(RESULTS, file = rsf.out.path)
  
     },silent = FALSE)  
  }

#bind all results:
full.results = do.call(rbind, reslist)

#save results
#-------------------------------------------------------------------------------
result.dir = './data/rsf/'
#-------------------------------------------------------------------------------
write.csv(full.results, paste0(result.dir, period, '_rsf_results_230227.csv'), row.names = FALSE)
  

#read the csv files back in to merge them:
setwd('./caribou_movement_ecology')
p2 = read.csv("./data/rsf/prior2_rsf_results_230226.csv")
p1 = read.csv("./data/rsf/prior1_rsf_results_230227.csv")
du = read.csv("./data/rsf/during_rsf_results_230227.csv")
af = read.csv("./data/rsf/after_rsf_results_230228.csv")

rsf.output = rbind(p2, p1, du, af)
write.csv(rsf.output, "./data/rsf/ALL_PERIOD_RSF_results_230228.csv",
          row.names = FALSE)

#-------------------------------------------------------------------------------
#use ctmm::mean() on the list of Rdas:
rm(list = ls())
gc()

library(ctmm)

#set period
period = 'after'
setwd('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/git_caribou/caribou_movement_ecology')
#folder with UD
rsf_path = paste0('./data/rsf/', 
                  period)

setwd(rsf_path)
#load UDs from that folder:
rsf.dir = list.files(pattern="*.rda", all.files=TRUE, 
                     full.names=FALSE)

#double check directory
getwd()

rsf <- list()

for(i in 1:length(rsf.dir)){
  load(rsf.dir[[i]])
  rsf[[i]] <- RESULTS
  
}

#-------------------------------------------------------------------------------
#metafor - 2023-03-01
rsf.list = list()
#
#check field names
#summary(rsf[[1]])

for(i in 1:length(rsf)){
  #i = 1
  #summary(rsf[[i]])
  rsf.df <- as.data.frame(period)
  rsf.df$ID = rsf[[i]]@info$identity
  rsf.df$elev <- summary(rsf[[i]], units = FALSE)$CI['elev (1/elev)',2]
  rsf.df$elev_min <- summary(rsf[[i]], units = FALSE)$CI['elev (1/elev)',1]
  rsf.df$elev_max <- summary(rsf[[i]], units = FALSE)$CI['elev (1/elev)',3]
  rsf.df$elev.cov = rsf[[i]]$COV['elev','elev']
  #      #.
  rsf.df$slope <- summary(rsf[[i]], units = FALSE)$CI['slope (1/slope)',2]
  rsf.df$slope_min <- summary(rsf[[i]], units = FALSE)$CI['slope (1/slope)',1]
  rsf.df$slope_max <- summary(rsf[[i]], units = FALSE)$CI['slope (1/slope)',3]
  rsf.df$slope.cov = rsf[[i]]$COV['slope','slope']
  #      #.
  rsf.df$heli.ten <- summary(rsf[[i]], units = FALSE)$CI['heli.ten (1/heli.ten)',2]
  rsf.df$heli.ten_min <- summary(rsf[[i]], units = FALSE)$CI['heli.ten (1/heli.ten)',1]
  rsf.df$heli.ten_max <- summary(rsf[[i]], units = FALSE)$CI['heli.ten (1/heli.ten)',3]
  rsf.df$heli.ten.cov = rsf[[i]]$COV['heli.ten','heli.ten']
  #      #.
  rsf.df$proj.age <- summary(rsf[[i]], units = FALSE)$CI['proj.age (1/proj.age)',2]
  rsf.df$proj.age_min <- summary(rsf[[i]], units = FALSE)$CI['proj.age (1/proj.age)',1]
  rsf.df$proj.age_max <- summary(rsf[[i]], units = FALSE)$CI['proj.age (1/proj.age)',3]
  rsf.df$proj.age.cov = rsf[[i]]$COV['proj.age','proj.age']
  rsf.df$major.cov = rsf[[i]]$COV['major','major']
  
  rsf.list[[i]] = rsf.df
  }
full.results = do.call(rbind, rsf.list)
#

after = full.results
#clear these so there's no confusion for the next period
rm(list = c('period','rsf_path', 'rsf.dir', 'rsf', 'rsf.list'))
#
data = rbind(during, after, prior1, prior2)
#
write.csv(data, 
          '../all_periods_RSF_metafor_230301.csv', 
          row.names = FALSE)
#-------------------------------------------------------------------------------
#2023 EOF git_caribou
#-------------------------------------------------------------------------------
#check all rsfs have data
#this block finds rsfs with Inf in the high estimate for each parameter and
#prints them, but that does not seem to be causing the error:
#Error in eigen(S) : infinite or missing values in 'x'

#summary(rsf[[3]])
#z = c()
#for(k in 1:length(rsf)){
#
#  #if(is.infinite(summary(rsf[[k]], units = FALSE)$CI['error all (meters)',3])){
##  print(paste(k, rsf[[k]]@info$identity, summary(rsf[[k]], units = FALSE)$CI['error all (meters)',3], sep = '--'))
#  print(summary(k))
##  z <- c(z,k) #rsf[[k]]@info$identity
#}#}
#z
##prior2 missing some values (Inf) causing mean() to not work:
#rsf.inf = rsf[-(c(1,10,14,16,17))]
#q = c()
#for(l in 1:length(rsf.inf)){
#    print(paste(l, rsf.inf[[l]]@info$identity, summary(rsf.inf[[l]], units = FALSE)$CI['error all (meters)',3], sep = '--'))
#    q <- c(q,l)
#  }
##-------------------------------------------------------------------------------
##well all that above to remove Inf values did nothing, and for p2 26 and 35 were causing an error#
#
##the following subsets allow mean to run (2023-02-28):
#rsf.p2 = rsf[-c(26,35)]
#rsf.p1 = rsf[-c(55:58)]
#rsf.du = rsf[-c(4,47,49:51)]
#rsf.af = rsf[-c(3,6,7,29:31,38)]
#
##double check we're picking the correct values
#summary(rsf[[5]])
#
##this is tedious, but go through to see where this error is occurring..
#rsf.test = rsf[c(1,2,4,5,8:28,32:37,39)]
##change this to reflect each of the period lists above
#x = mean(rsf.af)
#
#per <- as.data.frame(period)
#
##      #Get HR results
##---------------------------------------------------------------------------------


