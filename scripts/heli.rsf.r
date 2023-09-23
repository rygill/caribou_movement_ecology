#
#code to run glm on human from Strava
#

#-------------------------------------------------------------------------------
#2023-02-05
#-------------------------------------------------------------------------------
#VERSION 2 - using strava data:
library(ggplot2)
library(dplyr)

library(sf)
library(terra)

library(mgcv)
library(caret)

rm(list = ls())
gc()

#this dataset was subset to remove points within 100m of a stream to remove those
#locations where skiers congregate for pickup. In addition, a 200 m buffer was 
#drawn around the points, and the location of the maximum elevation was extracted.
#This maximum location was then buffered by 250m to remove locations on mountain tops,
#or on the opposite side of ridges from where skiers descended
#These were done to select only points that fell within the run itself, to avoid those
#locations that were in the pickup and dropoff zones.
#different subset data
#20m
#strava.verts = st_read('./data/strava_run_points_20m_creek_250m_peak_230307.shp')
#100 m
#strava.verts = st_read('./data/strava_run_points_100m_creek_250m_peak_230307.shp')

#decided to use all data without removing those locations
#all data:
strava.verts = st_read('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/disturbance/heli_runs/digitized_points/strava_centroids_win_tenures_Albers.shp')
#check for empty geometries
any(is.na(st_dimension(strava.verts)))
st_is_valid(strava.verts)
#clean up empty geometries (only full dataset has t his problem)
strava.verts = strava.verts %>% filter(!st_is_empty(.))


head(strava.verts)
strava.verts$type = as.factor(strava.verts$type)
strava.verts$tenure = as.factor(strava.verts$tenure)
levels(strava.verts$tenure)

#take only those points with home ranges assigned
strava.verts = strava.verts[!(is.na(strava.verts$tenure)),]
levels(strava.verts$type)
#limit to those points that we know are not ski touring lodges that do not fall
#within a heli-ski tenure
strava.verts = strava.verts[is.na(strava.verts$type),]

min(strava.verts$DN)

#read in heli ski tenures to clip everything to:
tenures = st_read('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/disturbance/heliski_tenures.shp')

#-------------------------------------------------------------------------------
#buffer tenures by 40m before clipping - this allows us to clip to the tenures
#for the prediction which removes spurious values along the edge of the raster.
clp = st_buffer(tenures, 40)

#avcan = st_read('C:/Users/Ryan/OneDrive/geodatabase/boundaries/avcan_regions.shp')
#rasterize heli home range
#avcan.raster = rasterize(avcan, elev, field = 'ras_value', background = 0, filename = './data/rasters/large_extent/regions.tif')

#-------------------------------------------------------------------------------
#load rasters:
#large extent unstandardised rasters:
elev = rast('./data/rasters/large_extent/elev_2023.tif')
slope = rast('./data/rasters/large_extent/slope_2023.tif')
roughness = rast('./data/rasters/large_extent/roughness_2023.tif')
aspect = rast('./data/rasters/large_extent/aspect_2023.tif')
tpi = rast('./data/rasters/large_extent/tpi_2023.tif')
dist.ldge = rast('./data/rasters/large_extent/heli_lodge_prox_2023.tif')
heli.hr.rast = rast('./data/rasters/large_extent/heli_home_range_2023.tif')
region = rast('./data/rasters/large_extent/regions.tif')
slope_length = rast('./data/rasters/large_extent/slope_length.tif')
wind = rast('./data/rasters/large_extent/wind.tif')

names(slp.len)
hist(region$ras_value)

#tpi = terrain(elev, 'TPI')
#alternative to terrain()
#tpi <- focal(elev, w=matrix(7/7, nc=7, nr=7), fun = \(elev) elev[5] - mean(elev[-5]))
#names(tpi) = 'tpi'
#writeRaster(tpi, './data/rasters/large_extent/tpi_2023.tif')

#flow = terrain(elev, 'flowdir')

#res(elev)
#aggregate to 100m
#elev = aggregate(elev, fact = 4, fun = 'mean')
#slope = aggregate(slope, fact = 4, fun = 'mean')
#roughness = aggregate(roughness, fact = 4, fun = 'mean')
#aspect = aggregate(aspect, fact = 4, fun = 'mean')
#dist.ldge = aggregate(dist.ldge, fact = 4, fun = 'mean')
#heli.hr.rast = aggregate(heli.hr.rast, fact = 4, fun = 'mean')
#heli.hr.rast = aggregate(heli.hr.rast, fact = 4, fun = 'mean')

#hist(flow)


#operations on input rasters for gam
#elevation long right tailed so square root
#elev = sqrt(elev)
#slope sqr root
#slope = sqrt(slope)
#deal with long tails in tpi
#tpi = sqrt((tpi + 150) / 100)

#stack and extract
rstack = c(elev, slope, roughness, aspect, tpi, dist.ldge, heli.hr.rast, region, slope_length, wind) #, flow
rstack = crop(rstack, clp) #clp is tenures buffered by 40 m
names(rstack)

#extract
strava.attrs = extract(rstack, strava.verts, bind = TRUE, na.rm = TRUE)

names(strava.attrs)

#available
stk.df = terra::as.data.frame(rstack)
names(stk.df)

#clean up datasets to join them:
strava.attrs = strava.attrs[,c("DN",'ras_value',
                               "elev_2023","slope_2023","roughness_2023", "aspect", 'tpi', 
                               "heli_lodge_prox_2023", "heli_home_range_2023", "slope_length", "wind")]

#assign DN value of 0 for non-detections:
stk.df$DN = 0
stk.df = stk.df[,c("DN",'ras_value',
                   "elev_2023","slope_2023","roughness_2023", "aspect", 'tpi',
                   "heli_lodge_prox_2023", "heli_home_range_2023", "slope_length", "wind")]

#check names again
names(strava.attrs)
names(stk.df)

#sample randomly within each area
stk.df1 <- stk.df[stk.df$ras_value > 0,] %>% 
  group_by(ras_value) %>% 
  slice_sample(n=300000)

#some data clean up and formatting
stk.df = as.data.frame(stk.df1)
rm(stk.df1)
strava.attrs = as.data.frame(strava.attrs)


#merge known and available:
data = rbind(strava.attrs, stk.df)
data = na.omit(data)
data$ras_value = as.factor(data$ras_value)

#examine data
head(data)
tail(data)
max(data$heli_lodge_prox_2023)

gc()

#subsample to get train and test:
dt = sort(sample(nrow(data), nrow(data)*.8))
train = data[dt,]
test = data[-dt,]


hist(train$slope_2023, xlab = "sqrt( slope )", main = "Histogram of Slope")
hist(train$elev_2023, xlab = 'sqrt( elevation )', main = "Histogram of elevation")
hist(train$tpi, xlab = "sqrt((tpi + 150) / 100)", main = 'Histogram of topographic\nposition index')
hist(train$heli_lodge_prox_2023)

#-------------------------------------------------------------------------------
#plots of use as a function of each predictor:
#*plots of use as a function of each predictor:
plot.g = function(a,b,c)
{
  ggplot() +
    geom_point(data = a, aes(x = b, y = c)) +
    geom_smooth(data = a, aes(x = b, y = c), 
                method = 'loess', se = FALSE) +
    xlab(b) +
    ylab("Intensity of Use")
}

plot.g(strava.attrs, strava.attrs$elev_2023, strava.attrs$DN)
plot.g(strava.attrs, strava.attrs$slope_2023, strava.attrs$DN)
plot.g(strava.attrs, strava.attrs$aspect_2023, strava.attrs$DN)
plot.g(strava.attrs, strava.attrs$slope_length, strava.attrs$DN)


max(train$elev)
names(train)

names(rstack)

#collinearity
library(corrplot)
dat.cor = data[,c("elev_2023","slope_2023","roughness_2023", "aspect",
                  "tpi", "heli_lodge_prox_2023", "slope_length")]
names(dat.cor) = c("Elevation","Slope","Roughness", "aspect",
                   "tpi", "Lodge Distance", "slope length")
corrplot(cor(dat.cor), method = 'number')

##------------------------------------------------------------------------------
#GAM
##------------------------------------------------------------------------------
library(mgcv)
m_1 <-
  bam(DN ~ 
        s(elev_2023, k = 5, bs = "tp") + #fs for random effects
        s(slope_2023, k = 5, bs = "tp") + 
        s(aspect, k = 5, bs = "tp") +
        s(tpi, k = 5, bs = "tp") + 
        s(heli_lodge_prox_2023, k = 10, bs = "tp") + #k = 10 results in better performance
        s(slope_length, k = 10, bs = "tp") +
        s(wind, k = 5, bs = "tp") +
        heli_home_range_2023, 
      family = poisson(link = "log"),
      data = train,
      method = 'fREML', # fast REML
      discrete = TRUE,  # discretize the posterior for faster computation
      # use multiple threads, print progress (excessive threads slow it down)
      #parallel::detectCores(logical = FALSE)
      control = gam.control(nthreads = 7, trace = TRUE), 
      na.action = 'na.fail')

summary(m_1)

#plot(m_1$fitted.values, m_1$residuals, xlab = "dfd", ylab = "dfdd")

#m_1$deviance-m_1$null.deviance
#plot model parameters
plot.gam(m_1, rug = FALSE, pages = 1, scale = 0,
     residuals = FALSE, ylab = "smoothed term")

#validation
# % deviance explained by our best model, often anything >30% is good.
#import 20% holdout - cross validation points
#####################################
test$pred.gam <- predict(m_1,newdata=test ,type="response")
#DN is at least 1, so we know 0 are all available points. Recode to make a boxplot
test$pt_type = ifelse(test$DN == 0,"Available","Used")



#####################################
ggplot() +
  geom_histogram(data = test, aes(pred.gam), bins=200)

#confusion matrix:
cm = table(test$pred.gam > 1, test$pt_type)
cm[2,2]/sum(cm[,2]) #32%
cm[1,1]/sum(cm[,1]) #97%

#precision 0.25
pr = sum(cm[,2]) / (sum(cm[,2]) + cm[2,1])
#accuracy 0.78
ac = (sum(cm[,2]) + sum(cm[,1])) / (sum(cm[,2]) + sum(cm[,1]) + cm[2,1] + cm[1,2])
####################################

pred.gam = predict(rstack, m_1, type = 'response')
#plot  rug = false pages = 1 scale = 0
setMinMax(pred.gam)
minmax(pred.gam)

pred.gam.clp = crop(pred.gam, tenures)
#pred.gam.clp1 = exp(pred.gam.clp)

writeRaster(pred.gam.clp, './data/rasters/predictions/run_pred_230420_gam_FULL_DATA.tif',
            overwrite = TRUE)

pred.val = 
ggplot() +
  geom_boxplot(data = test, aes(x = pt_type, y = pred.gam)) +
  ggtitle("Model Evaluation") +
  coord_cartesian(ylim = c(0,15)) +
    theme_bw() +
  ylab(expression(paste("Predicted value of ski terrain"))) +
  xlab("Available and Used locations") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=10, family = "sans"),
          axis.title.x = element_text(size=10, family = "sans"),
          axis.text.y = element_text(size=8, family = "sans"),
          axis.text.x  = element_text(size=8, family = "sans"),
          plot.title = element_text(size = 10, family = "sans"))

ggsave(pred.val,
       width = 2000, height = 1400, units = "px", dpi = 300,
       file="./figures/model_evaluation.png")

#normalize the values so they're on the same scale:
test$norm.DN = test$DN/max(test$DN)
test$norm.pred = test$pred.gam/max(test$pred.gam)

#scatter plot of predicted vs observed:
ggplot() +
  geom_point(data = test, aes(x = norm.DN, y = norm.pred), alpha = 0.25) +
  geom_smooth(data = test, aes(x = norm.DN, y = norm.pred),method = 'lm') +
  ggtitle("Model Evaluation") +
  theme_bw() +
  ylab(expression(paste("Normalized predicted value of ski terrain"))) +
  xlab("Normalized color intensity") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans"),
        axis.title.x = element_text(size=10, family = "sans"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(size = 10, family = "sans"))

ggsave(width = 2000, height = 1400, units = "px", dpi = 300,
       file="./figures/pred_vs_obs_230412.png")

####################################################

#read points
#pts = st_read('C:/Users/Ryan/OneDrive/temp/pointsID.shp')
#pts = pts %>%
#  group_by(id) %>%
#  summarise(max_elev = max(VALUE))
#
#write_sf(pts, 'C:/Users/Ryan/OneDrive/temp/points_max.shp')

#-------------------------------------------------------------------------------
#try again with binomial
#subsample to get train and test:
dt = sort(sample(nrow(data), nrow(data)*.8))
train2 = data[dt,]
test2 = data[-dt,]

#set DN to reflect used/available:
train2$detect = ifelse(train2$DN == 0,0,1)
test2$detect = ifelse(test2$DN == 0,0,1)

m_2 <-
  bam(detect ~ 
        s(elev_2023, k = 5, bs = "tp") + #fs for random effects
        s(slope_2023, k = 5, bs = "tp") + 
        #s(flowdir, k = 5, bs = "tp") + 
        s(aspect, k = 5, bs = "tp") +
        s(tpi, k = 5, bs = "tp") + 
        s(heli_lodge_prox_2023, k = 10, bs = "tp") + #k = 10 results in better performance
        heli_home_range_2023, 
      family = binomial(link = "logit"),
      data = train2,
      method = 'fREML', # fast REML
      discrete = TRUE,  # discretize the posterior for faster computation
      # use multiple threads, print progress (excessive threads slow it down)
      #parallel::detectCores(logical = FALSE)
      control = gam.control(nthreads = 7, trace = TRUE), 
      na.action = 'na.fail')

summary(m_2)

#plot model parameters
plot.gam(m_2, rug = FALSE, pages = 1, scale = 0,
         residuals = FALSE, ylab = "smoothed term")

#import 20% holdout - cross validation points
#####################################
test2$pred.gam2 <- predict(m_2,newdata=test2 ,type="response")
#DN is at least 1, so we know 0 are all available points. Recode to make a boxplot
test2$pt_type = ifelse(test2$detect == 0,"Available","Used")

ggplot() +
  geom_boxplot(data = test2, aes(x = pt_type, y = pred.gam2)) +
  ggtitle("Model Evaluation") +
  theme_bw() +
  ylab(expression(paste("Predicted value of ski terrain"))) +
  xlab("Available and Used locations") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans"),
        axis.title.x = element_text(size=10, family = "sans"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(size = 10, family = "sans"))


#normalize the values so they're on the same scale:
#test2$norm.DN = test2$DN/max(test2$DN)
#test2$norm.pred2 = test2$pred.gam/max(test2$pred.gam)

#evaluate accuracy
tab = table(test2$detect) #we have 28 1's so lets make a sample of 0 and 1 both 28 long
#how many rows to include?:
tab[2]

p <- sample(test2$pred.gam2[test2$detect ==1],tab[2])
a <- sample(test2$pred.gam2[test2$detect ==0],tab[2])
mv <- wilcox.test(p,a)
auc <- as.numeric(mv$statistic) / (length(p) * length(a))
auc  #0.5 = a random guess... I get 0.904 


#####################################
#confusion matrix:
cm = table(test2$pred.gam2 > 0.4, test2$pt_type)
cm[2,2]/sum(cm[,2]) #32%
cm[1,1]/sum(cm[,1]) #97%

#precision 0.77
pr = sum(cm[,2]) / (sum(cm[,2]) + cm[2,1])
#accuracy 0.92
ac = (sum(cm[,2]) + sum(cm[,1])) / (sum(cm[,2]) + sum(cm[,1]) + cm[2,1] + cm[1,2])
#proportion correct

pred.gam = predict(rstack, m_2, type = 'response')
#plot  rug = false pages = 1 scale = 0

pred.gam.clp = crop(pred.gam, tenures)
#pred.gam.clp1 = exp(pred.gam.clp)

writeRaster(pred.gam.clp, './data/rasters/predictions/run_pred_230414_gam_BINOMIAL.tif',
            overwrite = TRUE)
