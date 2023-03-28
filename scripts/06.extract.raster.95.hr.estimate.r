#script to extract raster covariates and assign to each of the 95% HR polygons
#for each period home range estimates for each animal were merged as single
#features into one shapefile for simplicity of extracting raster covariates

#using 'mean' as the function.

library(sf)
library(terra)
library(dplyr)
library(stringi)

#-------------------------------------------------------------------------------
#create a folder to hold temp files as they take up much HD space.
path = paste0(tempdir(),'//raster')
dir.create(path)

rm(list = ls())
gc()

#landscape
#rpath = "C:/Users/Ryan/OneDrive/ABMI/caribou_anthropause"

#original rasters
#elev = rast(paste0(rpath,'/spatial/dem/derived/cropped_20k/elev_20k.tif'))
#slope = rast(paste0(rpath,'/spatial/dem/derived/cropped_20k/slope_20k.tif'))
#proj.age = rast(paste0(rpath,'/spatial/dem/derived/cropped_20k/proj_age_AGE_20k.tif'))
#heli.ten = rast(paste0(rpath,'/spatial/dem/derived/cropped_20k/heli_tenures_20k.tif'))

  
#large extent to capture 2022  
elev = rast('./data/rasters/august/elev_220809.tif')
slope = rast('./data/rasters/august/slope_220809.tif')
#age
proj.age = rast('./data/rasters/august/proj_age_220809.tif')
#max(proj.age)
#disturbances:
heli.ten = rast('./data/rasters/august/heli_ten_220809.tif')

scl.elev = scale(elev)
scl.slope = scale(slope)
scl.age = scale(proj.age)


rstack1 = c(elev, scl.elev, scl.slope, scl.age, heli.ten)

#-------------------------------------------------------------------------------
period = 'prior1'
#to get this next file I merged all HR shapes in QGIS for a given period.
#shapefiles were selected for the 95% HR estimate:
dat.shp = st_read(dsn = './data/home_range/merged_95_HR',
                  layer = paste0(period,'_merged_95_HR_Albers'))

#dat.shp = dat.shp[dat.shp$est == 'est',]

#all raster values:
dat.extract = extract(rstack1, dat.shp, fun = 'mean', sp = TRUE, bind = TRUE, na.rm = TRUE)

dat.df = as.data.frame(dat.extract)
dat.df$BINOMIAL = period
#dat.df$name = as.character(dat.df$layer)
#length(dat.df$name)
dat.df$ID = stri_sub(dat.df$name, from = 0, to = 5)
dat.cols = as.data.frame(colnames(dat.df))

dat.a = dat.df[,c(8,9,4:7)]

write.csv(dat.a, paste0('./data/home_range/',period,'/',period,'_covariates_95_230214.csv'), row.names = FALSE)
rm(dat.shp)
#-------------------------------------------------------------------------------
rm(list = ls())
gc()

#join all those files together:
p1 = read.csv("./data/home_range/prior1/prior1_covariates_95_230214.csv")
p2 = read.csv("./data/home_range/prior2/prior2_covariates_95_230214.csv")
po = read.csv("./data/home_range/after/after_covariates_95_230214.csv")
du = read.csv("./data/home_range/during/during_covariates_95_230214.csv")
results = read.csv('./data/home_range/Caribou_Results_with_Error_230213.csv')

dat.b = rbind(p1, p2, du, po)

#join spatial covariates with home range output:
dat.b$joiner = paste(dat.b$BINOMIAL, '-',dat.b$ID)
results$joiner = paste(results$BINOMIAL, '-', results$ID)
dat.c = merge(results, dat.b, by = 'joiner')
dat.c$period = dat.c$BINOMIAL

#clean up:
dat.c = dat.c %>%
  dplyr::rename('period' = 'BINOMIAL.x',
         'ID'= 'ID.x')

#drop columns
dat.c$joiner = dat.c$ID.y = dat.c$BINOMIAL.y = NULL

#get herd for each animal:
herd = read.csv('./data/input_data/herd.csv')
dat.d <- merge(dat.c, herd, by = "ID")

dat.d$Year = as.numeric(round(dat.d$Year,0))

#clean up one more time:
col.dat = as.data.frame(colnames(dat.d))
dat.d = dat.d[,c(1,2,28,3:27)]

write.csv(dat.d, './data/home_range/Caribou_results_with_error_covariates_230214.csv', row.names = FALSE)

#EOF

