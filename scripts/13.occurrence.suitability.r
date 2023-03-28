rm(list = ls())
gc()


#library(devtools)
#install_github("ctmm-initiative/ctmm", dependencies = TRUE, force = TRUE)
library(ctmm)

#find out where tempdir is and make sure 'raster' is in there:
path = paste0(tempdir(),'/raster')
dir.create(path)


#-------------------------------------------------------------------------------
#read in landscape rasters:
#-------------------------------------------------------------------------------
#rpath = "C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/spatial/rasters_22/wgs/"
#-------------------------------------------------------------------------------
elev = raster('./data/rasters/large_extent/elev_2023.tif')
slope = raster('./data/rasters/large_extent/slope_2023.tif')
proj.age = raster('./data/rasters/large_extent/proj_age_numeric_RCL_2023.tif')
heli.ten = raster('./data/rasters/large_extent/heli_tenure_2023.tif')
rsf = raster('./data/rasters/predictions/RATA_rsf/RATA_suitability_dplyr.tif')
#NAs throwing error, reclass to 0
#m <- c(NA,NA,0)
#rclmat <- matrix(m, ncol=3, byrow=TRUE)
#proj.age <- raster::reclassify(proj.age, rclmat, include.lowest=TRUE)
#raster::writeRaster(proj.age, './data/rasters/large_extent/proj_age_numeric_RCL_2023.tif')

#named list of rasters for rsf:
rlist = list(elev, slope, heli.ten, proj.age, rsf)
names(rlist) = c('elev', 'slope', 'heli.ten', 'proj.age', 'rsf')

#now that the rasters are in a list I can remove them:
rm(list = c('elev', 'slope', 'heli.ten', 'proj.age', 'rsf'))

#clean up:
gc()

rm(data)
rm(DATA)

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

#set period:
period = 'after'
data = droplevels(data[data$period == period,]) 

#there are some collars which the UD failed, remove these so the loop runs
#problematic outliers not detected during cleaning: 32611 and 44095.
#additional collars for which rsfs failed
#prior2
#data = data[!(data$individual.local.identifier %in% c(32611,44095,29102,30613,81288,81320,81321,81328,81344,81355,81365)),]
#prior1
#data = data[!(data$individual.local.identifier %in% c(32611,44095,81332,81361,81362)),]
#during
#data = data[!(data$individual.local.identifier %in% c(32611,44095,42144,81332,81355,81365,83782)),]
#after
data = data[!(data$individual.local.identifier %in% c(32611,44095,44098)),]


#Convert to a telemetry object
DATA <- ctmm::as.telemetry(data)


#-------------------------------------------------------------------------------
#model fits (run OCCURENCE without rlist):
#FITS_path = paste0("C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/home_range/", period, "/Fits_with_error")

#rsf:
FITS_path = paste0("./data/rsf/", period)
#-------------------------------------------------------------------------------

#folder with UD
setwd(FITS_path)
getwd()
#load UDs from that folder:
fit.dir = list.files(pattern="*.rda", all.files=TRUE, 
                    full.names=FALSE)
fits <- list()
#if using model fits use fits[[i]] <- FIT
#if using rsf use fits[[i]] <- RESULTS
for(i in 1:length(fit.dir)){
  load(fit.dir[[i]])
  fits[[i]] <- RESULTS
  
}
#set the wd back to the original so the other code works:
setwd('../../../')
getwd()

#create empty list to append animals to as the loop progresses
occur.list = list()

#*#set path for output of rsf models
#-------------------------------------------------------------------------------
occur_path = paste0('./data/occurrence/', 
                  period)
fig_fold = paste0(occur_path,'/figs')
#UD.fold = 'C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/analysis/occurrence/UD'
#-------------------------------------------------------------------------------
#create directories
dir.create(paste("./data/occurrence/", period, sep = ""))
dir.create(paste("./data/occurrence/", period, "/figs", sep = ""))


time.start = Sys.time()

for(j in 1:length(DATA)){
  try({
    cat("Working on individual ", j, " of ", length(DATA), "\n")
    #j = 10
    #Extract the current individual
    cilla <- DATA[[j]] 
    
    #Extract corresponding UD
    cilla.fit = fits[[j]] 
    cilla.fit@info$identity
    cilla@info$identity
    #project so they are the same
    #ctmm::projection(cilla) = ctmm::projection(cilla.fit)

    OCCUR <- ctmm::occurrence(cilla, cilla.fit, R = rlist)
   
    #save the plot:
    fig.path <- file.path(fig_fold,
                          paste('occur_plot_',period,'_',cilla@info[1], ".jpeg", sep = ""))
    
    #Save the graphic device's output as a png
    jpeg(file = fig.path,
        type = "cairo",
        units = "in",
        width = 6.81, height = 3,
        pointsize = 10,
        res = 600)
    plot(cilla, UD = OCCUR, level.UD = .95)
    dev.off() 
    
    #Assign the file path
    occur.out.path <- file.path(occur_path, paste0("occur_",period,'_', cilla@info[1], ".rda"))
    
    #And save:
    #rda:
    save(OCCUR, file = occur.out.path)
    
    #UD as raster:
    writeRaster(OCCUR, occur.out.path, format = 'GTiff',
                overwrite = TRUE, options=c("COMPRESS=NONE", "TFW=YES"), DF = 'PDF')
    
    #save 95% range estimate:
    shp.path = file.path(occur_path, paste0("occur_",period,'_', cilla@info[1], ".shp"))
    ctmm::writeShapefile(OCCUR, shp.path, level.UD=0.95, level=0.95)
    
    #multiply OCCUR by underlying rasters
    #convert UD to raster:
   ##########################################
    OCCURr = raster(OCCUR, DF = 'PDF')
    
    #read in the rasters in the loop:
    elev = raster('./data/rasters/large_extent/elev_2023.tif')
    slope = raster('./data/rasters/large_extent/slope_2023.tif')
    proj.age = raster('./data/rasters/large_extent/proj_age_numeric_RCL_2023.tif')
    heli.ten = raster('./data/rasters/large_extent/heli_tenure_2023.tif')
    rsf = raster('./data/rasters/predictions/RATA_rsf/RATA_suitability_dplyr.tif')
    
    #check projections:
    projection(OCCURr)
    projection(elev)
    #check the extents
    raster::extent(OCCURr) == raster::extent(elev)
    
    #nothing lining up so reproject.
    elev = raster::projectRaster(elev, OCCURr)
    slope = raster::projectRaster(slope, OCCURr)
    heli.ten =raster::projectRaster(heli.ten, OCCURr)
    proj.age = raster::projectRaster(proj.age, OCCURr)
    rsf = raster::projectRaster(rsf, OCCURr)
    
    #check
    projection(elev)
    projection(OCCURr)
    #check the extents
    raster::extent(OCCURr) == raster::extent(elev)
    
    #not working so set the extents:
    #elev.1 = raster::setExtent(elev,OCCURr)
    #slope.1 = raster::setExtent(slope, OCCURr)
    #heli.ten.1 = raster::setExtent(heli.ten, OCCURr)
    #proj.age.1 = raster::setExtent(proj.age, OCCURr)
    
    #also need to set the resolutions to match:
    OCCUR_res <- raster::resample(OCCURr, elev, method = 'bilinear')
    
    writeRaster(OCCUR_res, occur.out.path, format = 'GTiff',
                overwrite = TRUE, options=c("COMPRESS=NONE", "TFW=YES"), DF = 'PDF')
    
    #normalize raster so values fall 0-1:
    OCCUR_res = OCCUR_res/raster::cellStats(OCCUR_res, 'sum')
    
    #calculate:
    occur.df <- as.data.frame(period)
    occur.df$ID = cilla@info$identity
    occur.df$model = do.call(rbind,summary(cilla.fit)[1])
    occur.df$elev = raster::cellStats(OCCUR_res*elev, 'sum')
    occur.df$slope = raster::cellStats(OCCUR_res*slope, 'sum')
    occur.df$heli.ten = raster::cellStats(OCCUR_res*heli.ten, 'sum')
    occur.df$proj.age = raster::cellStats(OCCUR_res*proj.age, 'sum')
    occur.df$rsf = raster::cellStats(OCCUR_res*rsf, 'sum')
    
    #add them to a list to append others to:
    occur.list[[j]] = occur.df 
    
    #CHECK
    #Mean elevation at the GPS locations
    #buff <- SpatialPoints.telemetry(cilla)
    #raster::projection(buff) <- raster::projection(elev)
    #val = mean(raster::extract(x = elev, y = buff)) #1820 here, vs 1757 using occurrence
    #plot(elev)
    #plot(cilla, add = T)
    
    #rm('elev.1', 'slope.1', 'heli.ten.1', 'proj.age.1',
    #   'elev', 'slope', 'heli.ten', 'proj.age', 'rsf')
  
  },silent = FALSE)  
}

#bind all results:
full.results = do.call(rbind, occur.list)
write.csv(full.results, paste0('./data/occurrence/',period,'_OD_extract_results_230323.csv'), row.names = FALSE)

rm(full.results)

#read those saved results in:
p2 = read.csv("./data/occurrence/prior2_OD_extract_results_230323.csv")
p1 = read.csv("./data/occurrence/prior1_OD_extract_results_230323.csv")
du = read.csv("./data/occurrence/during_OD_extract_results_230323.csv")
af = read.csv("./data/occurrence/after_OD_extract_results_230323.csv")

dat = rbind(af, du, p1, p2)
write.csv(dat, './data/occurrence/OD_data_rsf_230323.csv', row.names = FALSE)

time.end = Sys.time()
difftime(time.end, time.start, units='mins')

#-------------------------------------------------------------------------------

#do some plotting:
library(dplyr)
library(ggplot2)
data = read.csv('./data/occurrence/OD_data_rsf_230323.csv')

herd = read.csv('C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause/collar_data/dataset/herd.csv')
herd$ID = as.factor(herd$ID)
data = merge(data, herd, by = 'ID')

#remove OUf because currently that model predicts incorrectly:
data = data[!(data$model == 'OUf error'),]

#wide to long:
rsf.param = tidyr::pivot_longer(data, cols = c(3:6), names_to = 'parameters',  values_to = 'param.est')

blue1 = "#4774A0"
blue2 = "#7BA2C9"
blueA = "#BCD5EE"
red = "#FF0000"

#ggplot() +
#  geom_boxplot(data = rsf.param, aes(x = period, y = param.est, color = period)) +
#  scale_colour_manual(values = c(blueA, red, blue1, blue2), labels = c('2022', '2021', '2020', '2019')) +
#  facet_wrap(~parameters, scales = 'free') 

library(lme4)
#for p values:
library(lmtest)
library(dplyr)
library(plyr)
library(MuMIn)
library(sjPlot)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#first rename after, so during is what the other levels are compared to:
data[which(data$period == "after"),'period'] <- "return"

#standardise slope and elevation:
data$elev <- scale(as.numeric(data$elev))
data$slope <- scale(as.numeric(data$slope))

#reorder the period levels so they plot in chronological order:
lvl.ordr.x = c('prior2','prior1','during','return')

#Forest age
age_mod <- glmer(proj.age ~ period + (1|ID:herd),
                 data = data,
                 family = Gamma('log'),
                 na.action = "na.fail")
age_null <- glmer(proj.age ~ 1 + (1|ID:herd),
                  data = data,
                  family = Gamma('log'),
                  na.action = "na.fail")
age_res <- anova(age_mod, age_null)
summary(age_mod)
dredge(age_mod)

age_P <- paste("p = ",round(age_res$`Pr(>Chisq)`[2],2), sep = "")

a <- 
  ggplot(data) +
  geom_boxplot(aes(x = factor(period, level = lvl.ordr.x), y = proj.age, fill = period),  size = 0.2, outlier.size = 0.2) +
  ggtitle("a)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=10, family = "sans"),
        plot.title = element_text(size=10, hjust = 0, family = "sans"),
        axis.text.y  = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.position = "none") +
        #legend.background = element_blank(),
        #legend.title = element_text(size=5, family = "sans", vjust = -1),
        #legend.title = element_blank(),
        #legend.text = element_text(size=8, family = "sans"),
        #legend.key.size = unit(0.2, "cm"),
        #legend.key = element_blank()) +
  ylab("Projected Age") +
  scale_fill_manual(labels = c("2021", "2020", "2019", "2022"), values = c(red, blue1, blue2, blueA)) +
  scale_x_discrete(labels = c("2019", "2020", "2021", "2022")) +
  geom_text(x=0.95, y=1.4, label=age_P, family = "sans", size = 3)


#Slope
slope_mod <- lmer(slope ~ period + (1|ID:herd),
                  data = data,
                  na.action = "na.fail",
                  REML = F)

slope_null <- lmer(slope ~ 1 + (1|ID:herd),
                   data = data,
                   na.action = "na.fail",
                   REML = F)

slope_res <- anova(slope_mod, slope_null)
summary(slope_mod)
dredge(slope_mod)

slope_P <- paste("p = ",round(slope_res$`Pr(>Chisq)`[2],3), sep = "")


b <- 
  ggplot(data) +
  geom_boxplot(aes(x = factor(period, level = lvl.ordr.x), y = slope, fill = period),  size = 0.2, outlier.size = 0.2) +
  ggtitle("b)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=10, family = "sans"),
        plot.title = element_text(size=10, hjust = 0, family = "sans"),
        axis.text.y  = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.position = 'none') + #"top",
        #legend.background = element_blank(),
        #legend.title = element_text(size=5, family = "sans", vjust = -1),
        #legend.title = element_blank(),
        #legend.text = element_text(size=8, family = "sans"),
        #legend.key.size = unit(0.2, "cm"),
        #legend.key = element_blank()) +
  ylab("Standardised slope") +
  scale_fill_manual(labels = c("2021", "2020", "2019", "2022"), values = c(red, blue1, blue2, blueA)) +
  scale_x_discrete(labels = c("2019", "2020", "2021", "2022")) +
  geom_text(x=0.95, y=-3.5, label=slope_P, family = "sans", size = 3)



#Elevation
elevation_mod <- lmer(elev ~ period + (1|ID:herd),
                      data = data,
                      na.action = "na.fail",
                      REML = F)

elevation_null <- lmer(elev ~ 1 + (1|ID:herd),
                       data = data,
                       na.action = "na.fail",
                       REML = F)

elevation_res <- anova(elevation_mod, elevation_null)
summary(elevation_mod)
dredge(elevation_mod)

elevation_P <- paste("p = ",round(elevation_res$`Pr(>Chisq)`[2],2), sep = "")


c <- 
  ggplot(data) +
  geom_boxplot(aes(x = factor(period, level = lvl.ordr.x), y = elev, fill = period),  size = 0.2, outlier.size = 0.2) +
  ggtitle("c)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=10, family = "sans"),
        plot.title = element_text(size=10, hjust = 0, family = "sans"),
        axis.text.y  = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.position="none") +
        #legend.background = element_blank(),
        #legend.title = element_text(size=5, family = "sans", vjust = -1),
        #legend.title = element_blank(),
        #legend.text = element_text(size=8, family = "sans"),
        #legend.key.size = unit(0.2, "cm"),
        #legend.key = element_blank()) +
  ylab("Standardised Elevation") +
  scale_fill_manual(labels = c("2021", "2020", "2019", "2022"), values = c(red, blue1, blue2, blueA)) +
  scale_x_discrete(labels = c("2019", "2020", "2021", "2022")) +
  geom_text(x=1, y=-3.5, label=elevation_P, family = "sans", size = 3)




#Tenures
heli_mod <- glmer(heli.ten ~ period + (1|ID:herd),
                  family = binomial(link = "logit"),
                  data = data,
                  na.action = "na.fail")

heli_null <- glmer(heli.ten ~ 1 + (1|ID:herd),
                   family = binomial(link = "logit"),
                   data = data,
                   na.action = "na.fail")

heli_res <- anova(heli_null, heli_mod)
summary(heli_mod)
dredge(heli_mod)
heli_res$`Pr(>Chisq)`[2]

heli_P <- paste("p = ",round(heli_res$`Pr(>Chisq)`[2],2), sep = "")


d <- 
  ggplot(data) +
  geom_boxplot(aes(x = factor(period, level = lvl.ordr.x), y = heli.ten, fill = period),  size = 0.2, outlier.size = 0.2) +
  ggtitle("d)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=10, family = "sans"),
        plot.title = element_text(size=10, hjust = 0, family = "sans"),
        axis.text.y  = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        legend.position="none") +
        #legend.background = element_blank(),
        #legend.title = element_text(size=5, family = "sans", vjust = -1),
        #legend.title = element_blank(),
        #legend.text = element_text(size=8, family = "sans"),
        #legend.key.size = unit(0.2, "cm"),
        #legend.key = element_blank()) +
  ylab("Heli-ski tenures") +
  scale_fill_manual(labels = c("2021", "2020", "2019", "2022"), values = c(red, blue1, blue2, blueA)) +
  scale_x_discrete(labels = c("2019", "2020", "2021", "2022")) +
  geom_text(x=0.95, y=0.1, label=heli_P, family = "sans", size = 3)


FIGS <- ggarrange(a,b,
                  c,d,
                  legend = 'none',
                  #common.legend = TRUE,
                  ncol = 2,
                  nrow = 2)

ggsave(FIGS,
       file="./UBC/thesis/figures/2022/Habitat_boxplots_221125.png",
       width = 6.23,
       height=4,
       units = "in",
       dpi = 600)



