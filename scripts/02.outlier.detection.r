#library(devtools)
#install_github("ctmm-initiative/ctmm", dependencies = TRUE)

rm(list = ls())

library(dplyr)
library(ctmm)

getwd()
setwd('../OneDrive/ABMI/caribou_anthropause/git_caribou/caribou_movement_ecology')

#read in data for moving window analysis:
dat1 = read.csv('./data/input_data/220526_moving_window_dat.csv')

str(dat1)

#-------------------------------------------------------------------------------
#just some exploration:
#how many collars are there available?
#collar.year = as.data.frame(table(dat1$sy,dat1$collar_id))
#collar.year = collar.year[collar.year$Freq > 0,] #66 collars active during the lockdown
#
#collar.use = collar.year %>%
#  group_by(Var2) %>%
#  summarize(n_distinct())
#
#dat3 = dat1[,c('collar_id', 'herd')]
#dat3 = dat3 %>%
#  group_by(herd) %>%
#  summarize(n_distinct(collar_id, herd))
#-------------------------------------------------------------------------------
#carpentry:
dat1$collar_id = as.factor(dat1$collar_id)
dat1$date_time = as.POSIXct(dat1$date_time, format = '%Y-%m-%d %H:%M', tz='UTC')
dat1$mort.status = as.factor(dat1$mort.status)
dat1$herd = as.factor(dat1$herd)
#dat1$period = as.factor(dat1$period)
#dat1$sy = as.factor(dat1$sy)

#check the date range is correct:
max(dat1$date_time)
min(dat1$date_time)


#format names to match required:
dat2 = plyr::rename(dat1, c('collar_id' = 'individual.local.identifier','date_time' = 'timestamp', 'latitude' = 'location-lat', 'longitude' = 'location-long'))
str(dat2)

#duplicated positions?
summary(duplicated(dat2))

#save this file that's formatted for the telemetry, but unfiltered for outliers:
#write.csv(dat2,'C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/220526_data_formatted_for_tele.csv', row.names = FALSE)


#########################****************************#######################
#ALL SY
telem = as.telemetry(dat2, timeformat = '%Y-%m-%d %H:%M:%S', timezone = 'UTC')

#ALL SY
#CHECK FOR OUTLIERS:
telem.out = outlie(data = telem, plot = TRUE, by = 'd')
telem.x = plyr::ldply(telem.out, rbind)
telem.x$speed = telem.x$speed * 86.4
telemetry.check = cbind(dat2, telem.x)

#ALL SY
#make some rules for outliers:
telemetry.check$outlier = ifelse(telemetry.check$speed >= 10, 1,
                                 ifelse((telemetry.check$speed < 10 & telemetry.check$speed >= 5), 2, 0))

#also identify DOP values >= 10
telemetry.check$outlier = ifelse(telemetry.check$DOP > 10, 3, telemetry.check$outlier)
#outliers are then: 0: normal, 1: outlier, 2: possible outlier, 3: outlier based on DOP value


#mort status check
#check for collars with strange mortality signals interspersed, or for collars that are only alive for a short time.
ggdat = telemetry.check[telemetry.check$mort.status == 'mortality',c(2,3,7)]
ggdat = droplevels(distinct(ggdat, individual.local.identifier, .keep_all = TRUE))
levels(ggdat$individual.local.identifier)
mort.check = telemetry.check[telemetry.check$individual.local.identifier %in% 
                               c("22561", "25315", "29094", "37670", "37672", "42295", "81289", "81297", "81318", "81347"),]
mort.check$date = as.Date(mort.check$timestamp)
library(ggplot2)
#geom point is an inobvious choice, but it shows intermittent mortality signals.
ggplot() +
  geom_point(data = mort.check, aes(x = date, y = individual.local.identifier, color = mort.status))
ggsave(file = "mortality_check.jpg", path = './figures', width = 12, height = 8)
#22561 REMOVE - only 5 days alive within period
#37672 REMOVE - unpredictable mortality reports
#37670 remove for anthropause due to intermittent mortality signals
#81289 has some anomalous mort signals but in May so can keep

#remove the erroneous ones:
telemetry.check = telemetry.check[!(telemetry.check$individual.local.identifier %in% c(22561, 37672)),]
#remove 37670 for anthropause period but keep it for others:
telemetry.check = telemetry.check[!(telemetry.check$individual.local.identifier == 37670 & telemetry.check$during == 1),]
telemetry.check = telemetry.check[!(telemetry.check$individual.local.identifier == 37670 & telemetry.check$after == 1),]

#calculating home ranges reveals a couple of outliers in the HR during data not captured by their DOP values. 
#Collar 44095 has 2 points 60km and 205km away from main points over 24 hrs with an outlier vaue of 3.
telemetry.check = telemetry.check[!(telemetry.check$id %in% c(32597, 32598)),]

#write raw data:
#write.csv(telemetry.check, './data/input_data/220720_telemetry_check_raw_dat.csv')

#ALL SY
filtered.data = telemetry.check[telemetry.check$outlier %in% c(0,2,3) & telemetry.check$mort.status == 'normal',]
#write.csv(filtered.data, './collar_data/dataset/analysis/220219_clean_data_formatted_for_tele.csv', row.names = FALSE)

#write filtered data for moving window analysis
write.csv(filtered.data, './data/input_data/220526_moving_window_formatted_for_tele_dat.csv', row.names = FALSE)



#eof

