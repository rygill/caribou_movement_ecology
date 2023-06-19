library(sf)
library(stringr)

getwd()
#-----------------------CHANGE PERIOD------------------------------------------#
period = 'after'
#-----------------------CHANGE PERIOD------------------------------------------#

shp_path = paste0("./data/home_range/", period, "/UDs_with_error")
#-------------------------------------------------------------------------------

#folder with UD
setwd(shp_path)
#load UDs from that folder:
shp.dir = list.files(pattern="*.shp", all.files=TRUE, 
                    full.names=FALSE)

#read the shapefiles in
shps = lapply(shp.dir, st_read)

#change the directory back to the project directory
getwd()
setwd('../../../../')
getwd()

period.shp = dplyr::bind_rows(shps)
period.shp = period.shp[str_detect(period.shp$name, "est"),]

period.shp = st_transform(period.shp, crs = st_crs('+init=EPSG:3005'))
st_crs(period.shp)
st_write(period.shp, dsn = './data/home_range/merged_95_HR', paste0(period, "_merged_95_HR_Albers_230619.shp"), 
                                            driver = 'ESRI Shapefile')

#EOF