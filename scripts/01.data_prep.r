####
#this script does two things - initially it selects a broad range of dates to be used for the moving 
#window analysis and then takes those refined values to select that data to be used in future analyses.
#

rm(list = ls())

library(dplyr)
library(lubridate)

library(raster)
library(rgdal)
library(ctmm)

#read in points
getwd()
setwd('../OneDrive/ABMI/caribou_anthropause/git_caribou/caribou_movement_ecology')

#read in raw collar data
dat1 = read.csv('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/complete_data_220526.csv')
str(dat1)

#artifact row.name column, get rid of it.
#dat1 = dat1[,-1]
dat1$collar_id = as.factor(dat1$collar_id)
dat1$mort.status = as.factor(dat1$mort.status)
dat1$herd = as.factor(dat1$herd)
dat1$status = as.factor(dat1$status)
dat1$date_time = as.POSIXct(dat1$date_time, format = '%Y-%m-%d %H:%M', tz='UTC')


#an artifact from the original datasets was the different values for mortality. Those are updated here to be consistent:
levels(dat1$mort.status)

dat1 = dat1 %>%
  mutate(mort.status = recode(mort.status, 'Mortality no radius' = 'mortality', 
                              'Mortality No Radius' = 'mortality', 
                              'Nothing Detected' = 'normal'))

dat1$mort.status = as.factor(dat1$mort.status)
levels(dat1$mort.status)

#get rid of animals in the pen (I manually selected animals with locations in the pen and noted their ids), and dead animals:
#dat1 = dat1[dat1$status == 'wild' & dat1$mort.status == 'normal',]
#2022-02-11 keep mortality status to calibrate DOP values for error estimates:
dat1 = dat1[dat1$status == 'wild',]

##2022-02-11 update - using moving window to determine start and end dates for each
#animal for each year so expand the temporal window to catch all possible variability:
#earliest start:
#julian(as.Date('2020-12-01')) #18597
#latest end:
#julian(as.Date('2021-04-15')) #18732



#use days since January 1 to define the period for all years (so I can compare the same period among years):
#create a temp date with a common year:
#dat1$temp_date = update(dat1$date_time,year = 2000)

#October 31 = 303
#May 31 = 151

#December 1 = 335 days
#January 15 = 15
#January 20 = 20
#March 25 = 83
#April 1 = 90
#April 10 = 100 days
#April 15 = 115 days
#dat1$day_num = as.numeric(floor(difftime(dat1$temp_date, as.Date('2000-01-01'), units="days")))

#dat1 = dat1[dat1$collar_id== 30276,]
dat1$prior2 = ifelse(dat1$date_time >= '2018-10-31 00:00:00' & dat1$date_time <= '2019-05-31 00:00:00', 1, 0)
dat1$prior1 = ifelse(dat1$date_time >= '2019-10-31 00:00:00' & dat1$date_time <= '2020-05-31 00:00:00', 1, 0)
dat1$during = ifelse(dat1$date_time >= '2020-10-31 00:00:00' & dat1$date_time <= '2021-05-31 00:00:00', 1, 0)
dat1$after = ifelse(dat1$date_time >= '2021-10-31 00:00:00' & dat1$date_time <= '2022-05-31 00:00:00', 1, 0)

#define 1 column with the before, during and after periods:
#dat1$period = as.factor(ifelse(dat1$prior1 == 1, 'prior1',
#                               ifelse(dat1$prior2 == 'prior2',
#                                    ifelse(dat1$during == 1, 'during', 'post')))


#get those records that fall within the broader moving window specified above:
dat1$within_mov_window = rowSums(cbind(dat1$prior1, dat1$prior2, dat1$during, dat1$after), na.rm = T )


#select only those records among all years that fall within the October 31 - May 31 period:
dat.period = dat1[dat1$within_mov_window == 1,]

#write this csv to be used for the moving window:
dat.period.write = dat.period[,c(-9)]
#write data to be used in moving window analysis to determine timing:
write.csv(dat.period.write,'./data/input_data/220526_moving_window_dat.csv', row.names = FALSE)


#-------------------------------------------------------------------------------
#EOF