#Loop through windows of X days of locations, compute HR size for each window
#Plot ML point estimate of HR size and 95% CIs
#Mark Bidwell, Chris Fleming (mostly Chris!)

#Modified by Michael Noonan

#Last updated: Feb 14th 2022

rm(list = ls())
gc()

library(ctmm)
library(lutz)
library(tidyr)
library(dplyr)
#library(plyr)
library(stringr)
library(ggplot2)
library(lubridate)


getwd()

#read in data from 02.outlier.detection.r to establish when the arrive on winter range:

dat = read.csv('./data/input_data/220526_moving_window_formatted_for_tele_dat.csv')
dat = dat[dat$mort.status == 'normal',]
dat$timestamp = as.POSIXct(dat$timestamp, format = '%Y-%m-%d %H:%M', tz='UTC')

#calculating home ranges reveals a couple of extreme outliers in the HR during data not captured by their DOP values. 
#Collar 44095 has 2 points 60km and 205km away from main points over 24 hrs with an outlier vaue of 3.
#these were removed in 02, but double check
dat = dat[!(dat$id %in% c(32597, 32598)),] 

#-------------------------------------------------------------------------------
#what are the monitoring periods for each animal for each year:
#mon.times = dat %>%
#  group_by(individual.local.identifier) %>%
#  dplyr::summarise(first.date = min(timestamp), last.date = max(timestamp))
#mon.times$win = as.numeric(floor(difftime(mon.times$last.date, mon.times$first.date, units = c("days"))))
#-------------------------------------------------------------------------------
#dat.collars.d = dat %>%
#  group_by(herd, individual.local.identifier) %>%
#  distinct(individual.local.identifier)
#-------------------------------------------------------------------------------

#-----------------------CHANGE PERIOD------------------------------------------#
#change this to each period and run: prior2, prior1, during, after
period = 'prior2'
#-----------------------CHANGE PERIOD------------------------------------------#
tel.dat = droplevels(dat[dat$during == 1 & !(dat$individual.local.identifier %in% c(44095, 81365)),])
tel.dat = droplevels(dat[dat$during == 1 & dat$individual.local.identifier > 29103,])

data = as.telemetry(tel.dat, timeformat = '%Y-%m-%d %H:%M:%S', timezone = 'UTC') 
#prior1 fail: 81330 #29098 #81331 #81365
#prior2 fail: 22563, 22565, 29092, 29094, 30613, 81321
#during fail: #42144, 44094, 81365
#after fail: 42139, 81295
collars = tel.dat %>% distinct(individual.local.identifier)

#create directories
dir.create(paste0("./data/input_data/moving_window/", period, "/Fits"), recursive = TRUE)
dir.create(paste0("./data/input_data/moving_window/", period, "/UDs"), recursive = TRUE)

window.HR <- function(data, dt, win, Fig_Path, Result_Path, UD_path){
  DATA = data
  tryCatch({
  times <- seq(DATA$t[1], DATA$t[nrow(DATA)], by = dt)
  
  # chop off end of times
  times <- times[-((length(times)-win/dt+1):length(times))]
  })
  #initialize arrays for results
  AREAS_lo <- rep(NA,length(times))
  AREAS_ml <- rep(NA,length(times))
  AREAS_hi <- rep(NA,length(times))

  for(i in 1:length(times))
  {tryCatch
    ({
    print(paste((i),"of",length(times),"iterations"))
  
    # subset times t to t+win
    SUBSET <- DATA[times[i]<=DATA$t & DATA$t<=times[i]+win,]
  
    # analyse subset
    GUESS <- ctmm.guess(SUBSET,interactive=FALSE)
    FIT <- try(ctmm.select(SUBSET,GUESS))
    if(class(FIT)=="ctmm")
    {
      AKDE <- akde(SUBSET,FIT)
      
      #to not overwrite each iteration, need to add the date of each window as file name:
      fname = DATA[1]$timestamp
      f.name = lubridate::date(fname[i])
      
      #save AKDE 'date_animal_period_UD'
      akde.path <- file.path(UD_path, paste0(f.name, '_', period, "_UD_", DATA@info[1], ".rda"))
      save(AKDE, file = akde.path)
      
      #save fitted model 'date_animal_period.rda'
      mod_path <- file.path(Result_Path, paste0(f.name, '_', period,"_Fits_", DATA@info[1], ".rda"))
      save(FIT, file = mod_path)
      
      # store results in arrays
      SUMMARY <- summary(AKDE,units=FALSE)
      AREAS_lo[i] <- SUMMARY$CI[[1]]
      AREAS_ml[i] <- SUMMARY$CI[[2]]
      AREAS_hi[i] <- SUMMARY$CI[[3]]
    }
  })
  }
  # Plot ML point estimate of area and 95% CIs
  tryCatch({
  TIMES <- as.POSIXct(times, origin="1970-01-01",
                      tz = lutz::tz_lookup_coords(DATA$latitude[1],
                                            DATA$longitude[1]))
  RESULTS<-data.frame(TIMES,
                       AREAS_lo,
                       AREAS_ml,
                       AREAS_hi)
  names(RESULTS)[1]<-"Time"
  })
  ############################################
  #Plot all the results and save them as a png
  ############################################

  fig.path <- file.path(Fig_Path,
                        paste(DATA@info[1], ".png", sep = ""))
  
  
  #Save the graphic device's output as a png
  png(file=fig.path,
      type="cairo",
      units = "in",
      width = 6.81, height = 3,
      pointsize = 10,
      res = 600) #

  #Set the par to plot all on same screen
  par(mfrow=c(1,2),
      mgp = c(1.5, 0.5, 0),
      oma=c(0,0,0,0),
      mar=c(3,3,2,2),
      cex.lab=1.2,
      cex.main = 1,
      family = "serif") 

  #Plot the relocation data, coloured by time
  #Create a function that scales colours between red and blue
  rbPal <- colorRampPalette(c('#FF0000','#046C9A'))
  #Then create a variable that scales from red to blue between the two times
  DATA$Col <- rbPal(nrow(DATA))[as.numeric(cut(DATA$t,breaks = nrow(DATA)))]
  plot(DATA,
       col.grid = NA,
       pch = 20,
       cex = 0.2,
       col.DF = "#669543",
       col = DATA$Col,
       labels=FALSE)
  title(main = "a)", adj = 0)


  #Plot of the range estimates over time
  plot(RESULTS$AREAS_ml~RESULTS$Time,
       pch=19,
       cex=1.25,
       ylim=c(0,max(RESULTS$AREAS_hi,na.rm=TRUE)),
       ylab="Home Range Area",
       xlab="Date")
  arrows(RESULTS$Time,
         RESULTS$AREAS_lo,
         RESULTS$Time,
         RESULTS$AREAS_hi,
         length=0.05,
         angle=90,
         code=3)
  title(main = "b)", adj = 0)

  dev.off()

return(RESULTS)
}
  


results = lapply(data,
                  window.HR,
                  dt = 1 %#% 'day',
                  win <- 14 %#% 'day',
                  Fig_Path = paste0(getwd(),"/data/input_data/moving_window/", period),
                  Result_Path = paste0(getwd(),"/data/input_data/moving_window/", period,'/Fits'),
                  UD_path = paste0(getwd(),"/data/input_data/moving_window/", period,'/UDs'))


#-------------------------------------------------------------------------------
#load rda: 
getwd()
setwd(paste0(getwd(),'/data/input_data/moving_window/',period,'/UDs'))
ud.period = list.files(pattern="*.rda", all.files=TRUE, 
                     full.names=FALSE)

ud.period.uds <- list()

for(i in 1:length(ud.period)){
  load(ud.period[[i]])
  ud.period.uds[[i]] <- AKDE
}


hr.list = list()

for(i in 1:length(ud.period.uds)){
  #i = 10
  hr <- as.data.frame(period)
  hr$ID = ud.period.uds[[i]]@info$identity
  hr$start = str_sub(print(ud.period[[i]]),1,10)
  hr$hr_low = summary(ud.period.uds[[i]], units = FALSE)$CI[1,"low"]/1000000
  hr$hr_est = summary(ud.period.uds[[i]], units = FALSE)$CI[1,"est"]/1000000
  hr$hr_high = summary(ud.period.uds[[i]], units = FALSE)$CI[1,"high"]/1000000
  hr.list[[i]] = hr
  
}

ud.period.window = do.call(rbind, hr.list)
setwd('../../')
write.csv(ud.period.window, paste0(period,'_220604_moving_window.csv'), row.names = FALSE)

ggplot(ud.period.window) +
  geom_point(aes(x = start, y = hr_est))

#-------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

#read in those files saved above:
prior1 = read.csv('./data/input_data/moving_window/prior1_220415_moving_window.csv')
prior2 = read.csv('./data/input_data/moving_window/prior2_220416_moving_window.csv')
during = read.csv('./data/input_data/moving_window/during_220604_moving_window.csv')
after = read.csv('./data/input_data/moving_window/after_220527_moving_window.csv')

data = rbind(prior1, prior2, during, after)
data$start = as.Date(data$start)
data$period = as.factor(data$period)


herd = read.csv('./data/input_data/herd.csv')
herd$herd = as.factor(herd$herd)
data = merge(data, herd, by = 'ID')
data = data[data$herd %in% c('cn', 'cs', 'hr'),]
rm(herd)
rm(prior1)
rm(prior2)
rm(during)
rm(after)

write.csv(data, 'data/input_data/moving_window/MERGED_periods_mw_220624.csv', row.names = FALSE)
#-------------------------------------------------------------------------------

rm(list = ls())
gc()

library(tidyr)
library(lubridate)
library(ggplot2)

#read in line 252 again:
data = read.csv('./data/input_data/moving_window/MERGED_periods_mw_220624.csv')
data$ID = as.factor(data$ID)
data$period = as.factor(data$period)
data$start = as.Date(data$start)
data$herd = as.factor(data$herd)

dat = data %>%
  pivot_longer(!c(period, ID, herd, start), names_to = "HR_type", values_to = "95_estimate")
#dat$start = as.POSIXct(dat$start,  format = '%Y-%m-%d %H:%M', tz='UTC')
dat$week = week(dat$start)
dat$HR_type = as.factor(dat$HR_type)

dat$winter.week = ifelse(dat$week == 44, 1, 
                         ifelse(dat$week == 45, 2,
                                ifelse(dat$week == 46, 3,
                                       ifelse(dat$week == 47, 4,
                                              ifelse(dat$week == 48, 5, 
                                                     ifelse(dat$week == 49, 6,
                                                            ifelse(dat$week == 50, 7,
                                                                   ifelse(dat$week == 51, 8,
                                                                          ifelse(dat$week == 52, 9, 
                                                                                 ifelse(dat$week == 53, 10, 
                                                                                        dat$week + 10))))))))))
#1)
ggplot() +
  geom_smooth(dat = dat[dat$HR_type == "hr_est",], aes(x = winter.week, y = `95_estimate`, color = period), linetype = 'solid') +
  scale_x_discrete(name ="Weeks from October 31", limits = c(1:30)) +
  scale_y_continuous(name = '95% Home Range Estimate') +
  geom_vline(xintercept = c(9,25), linetype = 'dashed') #week 9 = Dec 24, week 25 = Apr 8
ggsave("C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/presentation/home_range.jpg")

#2)
ggplot() +
  geom_smooth(dat = dat[dat$HR_type == "hr_est",], aes(x = winter.week, y = `95_estimate`, color = period), linetype = 'solid') +
  scale_x_discrete(name ="Weeks from October 31", limits = c(1:30)) +
  scale_y_continuous(name = '95% Home Range Estimate') +
  geom_vline(xintercept = c(9,25), linetype = 'dashed') +  #week 9 = Dec 24, week 25 = Apr 8
  facet_wrap(period~herd, ncol = 2)

#-------------------------------------------------------------------------------

#average home range for period, herd and start time:
dat1 = data %>% 
  dplyr::group_by(period, herd, start) %>%
  dplyr::summarise(hr_ave_low = mean(hr_low), hr_ave_estimate = mean(hr_est), hr_ave_high = mean(hr_high))


dat2 = dat1 %>%
  pivot_longer(!c(period, herd, start), names_to = "HR_type", values_to = "HR_estimate")
dat2$week = week(dat2$start)
dat2$HR_type = as.factor(dat2$HR_type)
#calculate weeks from October 31 (winter.week)
dat2$winter.week = ifelse(dat2$week == 44, 1, 
                         ifelse(dat2$week == 45, 2,
                                ifelse(dat2$week == 46, 3,
                                       ifelse(dat2$week == 47, 4,
                                              ifelse(dat2$week == 48, 5, 
                                                     ifelse(dat2$week == 49, 6,
                                                            ifelse(dat2$week == 50, 7,
                                                                   ifelse(dat2$week == 51, 8,
                                                                          ifelse(dat2$week == 52, 9, 
                                                                                 ifelse(dat2$week == 53, 10, 
                                                                                        dat2$week + 10))))))))))

#update period so we have years instead of descriptors:
dat2$period = as.character(dat2$period)
dat2[which(dat2$period == "prior2"),"period"] <- "2019"
dat2[which(dat2$period == "prior1"),"period"] <- "2020"
dat2[which(dat2$period == "during"),"period"] <- "2021"
dat2[which(dat2$period == "after"),"period"] <- "2022"
dat2$period = as.factor(dat2$period)

#3)
ggplot() +
  geom_smooth(data = dat2[dat2$HR_type == 'hr_ave_estimate',], aes(x = winter.week, y = HR_estimate, color = period), linetype = 'solid') +
  geom_smooth(data = dat2[dat2$HR_type == "hr_ave_low",], aes(x = winter.week, y = HR_estimate, color = period), linetype = 'dashed', se = FALSE) +
  geom_smooth(data = dat2[dat2$HR_type == "hr_ave_high",], aes(x = winter.week, y = HR_estimate, color = period), linetype = 'dashed', se = FALSE) +
  scale_x_discrete(name ="Weeks from October 31", limits = c(1:30)) +
  scale_y_continuous(name = 'Log-scaled 95% Home Range Estimate', trans = 'log10') +
  geom_vline(xintercept = c(4,25), linetype = 'dashed') #week 9 = Dec 24, week 25 = Apr 8

#4)
ggplot() +
  geom_smooth(data = dat2[dat2$HR_type == "hr_ave_estimate" & dat2$herd == 'hr',], aes(x = winter.week, y = HR_estimate, color = period), linetype = 'solid') +
  scale_x_discrete(name ="Weeks from October 31", limits = c(1:30)) +
  scale_y_continuous(name = 'Log Scaled 95% Home Range Estimate', trans = 'log10') +
  geom_vline(xintercept = c(11,22), linetype = 'dashed', size = 1) +  #week 9 = Dec 24, week 25 = Apr 8 / 11-22 is original Jan 20 - Mar 25
  geom_vline(xintercept = c(6,15), linetype = 'dashed', color = '#F8766D', size = 0.75) + #hr during
  geom_vline(xintercept = c(6,15), linetype = 'dashed', color = '#00BA38', size = 0.75) + #hr prior1
  geom_vline(xintercept = c(6,15), linetype = 'dashed', color = '#619CFF', size = 0.75)   #hr prior2
  #facet_wrap(~herd, ncol = 1)

#-------------------------------------------------------------------------------
#check each herd's home range settlement period for each period:

#5)
ggplot() +
  geom_smooth(data = dat2[dat2$HR_type == "hr_ave_estimate",], aes(x = winter.week, y = HR_estimate, color = period),span = 0.1, linetype = 'solid') +
  geom_point(data = dat2[dat2$HR_type == "hr_ave_estimate",], aes(x = winter.week, y = HR_estimate, color = period), alpha = 0.5) +
  scale_x_discrete(name ="Weeks from October 31", limits = c(1:30)) +
  scale_y_continuous(name = 'Log Scaled 95% Home Range Estimate', trans = 'log10') +
  facet_wrap(period~herd, ncol = 2, scales = 'free')
ggsave("C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/presentation/home_range_by_period_herd_log_scaled_2022-06-03.jpg", width = 10, height = 12)
#-------------------------------------------------------------------------------
#5)
ggplot() +
  geom_smooth(data = dat2[dat2$HR_type == "hr_ave_estimate" & dat2$herd == 'cn',], aes(x = winter.week, y = HR_estimate, color = period), span = 0.5, linetype = 'solid') +
  geom_point(data = dat2[dat2$HR_type == "hr_ave_estimate" & dat2$herd == 'cn',], aes(x = winter.week, y = HR_estimate, color = period), alpha = 0.5) +
  scale_x_discrete(name ="Weeks from October 31", limits = c(1:30)) +
  scale_y_continuous(name = 'Log Scaled 95% Home Range Estimate', trans = 'log10') +
  ggtitle("Columbia North") +
  facet_wrap(~period, ncol = 2, scales = 'free')
ggsave("C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/presentation/home_range_CN_by_period_herd_log_scaled_2022-06-03.jpg", width = 10, height = 12)

#Estimated breakpoints (week from October 31) from figure 5
#HR Prior1: 8-15 (2019-12-17 to 2020-02-04)
#HR Prior2: 7-13 (2018-12-10 to 2019-01-21)
#HR During: 7-17 (2020-12-09 to 2021-02-18)
#HR After: 14-18 (2022-01-22 to 2022-02-25)

#CN P1: 13-23 (2020-01-17 to 2020-03-26)
#CN P2: 18-23 (2019-02-19 to 2019-04-01)
#CN DU: 15-22 (2021-01-29 to 2021-03-25)
#CN AF: 12-24 (2022-01-08 to 2022-04-08)

#CS P1: 10-20 (2019-12-31 to 2020-03-10)
#CS P2: 19-24 (2019-02-26 to 2019-04-08)
#CS DU: 10-24 (2020-12-30 to 2021-04-08)
#CS AF: 18-24 (2022-02-19 to 2022-04-08)

target = dat %>% 
  group_by(period, winter.week) %>%
  dplyr::summarise(min.week.date = min(start), max.week.date = max(start))

#6) recreate previous figure with v_lines for each period defined from above lines
#in excel I created a file that reflects those break dates:
dat_migr = read.csv('./data/input_data/moving_window/moving_window_weeks_220604.csv')
dat_migr$herd = as.factor(dat_migr$herd)
#change period values to be years, instead of descriptors:
dat_migr$period = as.character(dat_migr$period)
dat_migr[which(dat_migr$period == "prior2"),"period"] <- "2019"
dat_migr[which(dat_migr$period == "prior1"),"period"] <- "2020"
dat_migr[which(dat_migr$period == "during"),"period"] <- "2021"
dat_migr[which(dat_migr$period == "after"),"period"] <- "2022"
dat_migr$period = as.factor(dat_migr$period)

dat_migr$week = as.numeric(dat_migr$week)

herd = 'hr'
titl = 'Hart Ranges'
#set the data:
ggdat = dat2[dat2$HR_type == "hr_ave_estimate" & dat2$herd == herd,]
ggmigr= dat_migr[dat_migr$herd == herd,]

blue1 = "#3A8EB7"
blue2 = "#03AAFC"
blueA = "#BCD5EE"
red = "#FF0000"

mp = ggplot() +
  geom_smooth(data = ggdat, 
              aes(x = winter.week, y = HR_estimate, color = period), span = 0.5, linetype = 'solid') +
    geom_point(data = dat2[dat2$HR_type == "hr_ave_estimate" & dat2$herd == herd,], 
             aes(x = winter.week, y = HR_estimate, color = period), alpha = 0.7) +
  scale_x_continuous(name ="Weeks from October 31") +
  geom_point(data = ggdat[ggdat$HR_type == "hr_ave_estimate" & ggdat$herd == herd,], 
             aes(x = winter.week, y = HR_estimate), color = 'gray', pch = 21, alpha = 0.8) +
  scale_y_continuous(name = 'Log Scaled 95% Home Range Estimate (km^2)', trans = 'log10') +
  ggtitle(titl) +
  theme(legend.position = 'none') +
  scale_color_manual(values = c(blue2, blue1, red, blueA),
                      labels = c("2021", "2020", "2019", "2022")) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~period, ncol = 2, scales = 'free_x')

mp +
  theme(legend.position="none",
        axis.title.y = element_text(size=12, family = "sans"),
        axis.title.x = element_text(size=12, family = "sans"),
        axis.text.y = element_text(size=10, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans")) +
  geom_vline(data = ggmigr,
             aes(xintercept = week), linetype = 'dashed')


getwd()
ggsave(paste0("./figures/home_range_",herd,"_by_period_herd_log_scaled_2023-01-27.jpg"), 
       units = 'px', width = 2200, height = 1400)

#-------------------------------------------------------------------------------
#subset to windows determined above:
run.dat = read.csv('./data/input_data/220526_moving_window_formatted_for_tele_dat.csv')

#HART RANGES
run.hr = run.dat[run.dat$herd == 'hr',]
#HR prior1: 2019-12-17 to 2020-02-04
run.hr$prior1 = ifelse(run.hr$timestamp >= '2019-12-17' & run.hr$timestamp <= '2020-02-04', 1, 0)

#HR prior2: 2018-12-10 to 2019-01-21
run.hr$prior2 = ifelse(run.hr$timestamp >= '2018-12-10' & run.hr$timestamp <= '2019-01-21', 1, 0)

#HR during: 2020-12-09 to 2021-02-18
run.hr$during = ifelse(run.hr$timestamp >= '2020-12-09' & run.hr$timestamp <= '2021-02-18', 1, 0)

#HR after: 2022-01-22 to 2022-02-25
run.hr$after = ifelse(run.hr$timestamp >= '2022-01-22' & run.hr$timestamp <= '2022-02-25', 1, 0)

#COLUMBIA NORTH
run.cn = run.dat[run.dat$herd == 'cn',]
#CN prior1: 2020-01-17 to 2020-03-26
run.cn$prior1 = ifelse(run.cn$timestamp >= '2020-01-17' & run.cn$timestamp <= '2020-03-26', 1, 0) 

#CN prior2: 2019-02-19 to 2019-04-01
run.cn$prior2 = ifelse(run.cn$timestamp >= '2019-02-19' & run.cn$timestamp <= '2019-04-01', 1, 0)

#CN during: 2021-01-29 to 2021-03-25
run.cn$during = ifelse(run.cn$timestamp >= '2021-01-29' & run.cn$timestamp <= '2021-03-25', 1, 0)

#CN after: 2022-01-08 to 2022-04-08
run.cn$after = ifelse(run.cn$timestamp >= '2022-01-08' & run.cn$timestamp <= '2022-04-08', 1, 0)
  
#CENTRAL SELKIRKS  
run.cs = run.dat[run.dat$herd == 'cs',]
#CS prior1: 2019-12-31 to 2020-03-10
run.cs$prior1 = ifelse(run.cs$timestamp >= '2019-12-31' & run.cs$timestamp <= '2020-03-10', 1, 0)

#CS prior2: 2019-02-26 to 2019-04-08
run.cs$prior2 = ifelse(run.cs$timestamp >= '2019-02-26' & run.cs$timestamp <= '2019-04-08', 1, 0)

#CS during: 2020-12-30 to 2021-04-08
run.cs$during = ifelse(run.cs$timestamp >= '2020-12-30' & run.cs$timestamp <= '2021-04-08', 1, 0)

#CS after: 2022-02-19 to 2022-04-08
run.cs$after = ifelse(run.cs$timestamp >= '2022-02-19' & run.cs$timestamp <= '2022-04-08', 1, 0)

#subset for just the period within the winter range windows
#each of the dates is subset in 06.Batch_Run.r so these last 6 lines are not required, but added to calculate sample size.:
run.cn$resident = rowSums(run.cn[,c("prior1", "prior2", "during", "after")])
run.cn = run.cn[run.cn$resident == 1,]
run.cs$resident = rowSums(run.cs[,c("prior1", "prior2", "during", "after")])
run.cs = run.cs[run.cs$resident == 1,]
run.hr$resident = rowSums(run.hr[,c("prior1", "prior2", "during", "after")])
run.hr = run.hr[run.hr$resident == 1,]


#write each file separately to retain those dates specific to each herd:
write.csv(run.hr, './data/input_data/home_range/HR_220605_clean_data_formatted_for_tele.csv', row.names = FALSE)
write.csv(run.cn, './data/input_data/home_range/CN_220605_clean_data_formatted_for_tele.csv', row.names = FALSE)
write.csv(run.cs, './data/input_data/home_range/CS_220605_clean_data_formatted_for_tele.csv', row.names = FALSE)


#EOF



