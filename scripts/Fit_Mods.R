
# This script is a function for fitting movement models to relocation data
# It uses the best fit model to then estimate a an AKDE HR
# The fitted model, and UD are then exported for future use, and the area estimate +/- CIs is returned
# Finally, a 4 panel figure is returned for diagnostic purposes


#It requires as inputs: i) tracking data on 1 individual, as a telemetry object
#                       ii) A path to the location where the movement models will be saved
#                       iii) A path to the location where the UDs will be saved
#                       iv) A path to the location where the figures will be saved



#Written by Michael Noonan

#Last updated: Aug 6th 2020



library(ctmm)
library(geosphere)

#This function is used to fit movement models to gps data and calculate an akde home range estimate
#It requires ctmm

CTMM_FIT <- function(cilla,
                    Model_path,
                    UD_path,
                    Fig_Path,
                    speed = TRUE,
                    error = TRUE, #works with this set as error = TRUE
                    captures = NA){
  
  message("Fitting movement model for: ", unlist(cilla@info[1]))

  #First fit the movement model using the usual workflow
  #Generate the variogram
  vg.cilla <- variogram(cilla)
  
  #fit the variogram
  GUESS <- variogram.fit(vg.cilla, interactive= FALSE)
  
  if(error){GUESS$error <- TRUE}
  
  
  #Fit the movement models
  cilla.mods <- tryCatch(
    {
      cilla.mods <- ctmm.select(cilla,
                                CTMM = GUESS,
                                method = "pHREML",
                                control=list(method="pNewton",
                                             cores = -1),
                                trace = TRUE
                                )
    }, error=function(err) {
      message(paste("Model fitting using pHREML and optim failed, defaulting to ML"))
      
      #If pREML and/or the optimiser fail, this will try it again
      #using the standard ML without any optimisation
      
      cilla.mods <- ctmm.select(cilla,
                                CTMM = GUESS,
                                method = "ML")
      
      #The return here is essential, otherwise results aren't returned and it fails
      return(cilla.mods)
    }
  )
  
  #Save the best fit model
  #First get the best fit model if more than 1 have been fit
  if(class(cilla.mods) == "list") {FIT <- cilla.mods[[1]]} else {
    
    FIT <- cilla.mods
  }
  
  #Assign the file path
  ctmm.path <- file.path(Model_path, paste("Fits_", cilla@info[1], ".rda", sep = ""))
  
  
  #And save
  save(FIT, file = ctmm.path)
  
  
  
  ############################################
  #Calculate the AKDE home range area
  ############################################
  
  
  message("Estimating AKDE range area")

  
  #calculate the akde based on the best fit model
  cilla.akde <- akde(cilla,
                     FIT,
                     res = 50)
  
  
  #Assign the file path
  akde.path <- file.path(UD_path, paste("UD_", cilla@info[1], ".rda", sep = ""))
  shp.path = file.path(UD_path, paste0("UD_", cilla@info[1], ".shp"))
  
  #And save
  #save UD as raster:
  writeRaster(cilla.akde, akde.path, format = 'GTiff',
              overwrite = TRUE, options=c("COMPRESS=NONE", "TFW=YES"))
  #save rda:
  save(cilla.akde, file = akde.path)
  
  #save 95% range estimate:
  writeShapefile(cilla.akde, shp.path, level.UD=0.95, level=0.95)
  
  
  ############################################
  #Calculate mean speed
  ############################################
  if(speed == TRUE){
    message("Estimating median speed")
  
  #Estimate median speed (may be a mix of OU and OUF, so using robust = true)
  tryCatch(
    {
      SPEED <- speed(cilla, FIT, robust = TRUE, units = FALSE)
    }, error=function(err) {
      SPEED <- rep(NA, 3)
    })
  
  }

  ############################################
  #Plot all the results and save them as a pdf
  ############################################
  
  message("\n", "Saving the figures")
  
  #Assign the file path and name for saving the results
  fig.path <- file.path(Fig_Path,
                      paste("ctmm_", cilla@info[1], ".png", sep = ""))
  
  #Save the graphic device's output as a pdf
  png(file=fig.path,
      type="cairo",
      units = "in",
      width = 6.81, height = 6,
      pointsize = 10,
      res = 400) #
  
  #Set the par to plot all on same screen
  par(mfrow=c(2,2),
      mgp = c(1.5, 0.5, 0),
      oma=c(0,0,0,0),
      mar=c(3,3,2,2),
      cex.lab=1.5,
      cex.main = 2) 
  
  
  
  #Plot the zoomed in variogram 
#  tryCatch(
#    {
#      plot(vg.cilla, CTMM=FIT,family="serif", fraction = 0.005) 
#      title(main = "a)", family = "serif", adj = 0)
#      
#    }, error=function(err) {
#      
#      plot(vg.cilla, CTMM=FIT,family="serif", fraction = 0.05) 
#      title(main = "a)", family = "serif",  adj = 0)
#      
#    }
#    
#  )
  
  
  
  #Plot the variogram of the full time series
  plot(vg.cilla, CTMM=FIT,family="serif") 
  title(main = "a)", family = "serif",  adj = 0)
  
  
  #Plot a check for outliers
  OUTLIERS <- outlie(cilla, family="serif")
  title(main = "b)", family = "serif", adj = 0)
  
  
  plot(OUTLIERS, family="serif")
  title(main = "c)", family = "serif",  adj = 0)
  
  
  #Plot the AKDE range estimate, with the relocation data, coloured by time
  #Create a function that scales colours between red and blue
  rbPal <- colorRampPalette(c('#FF0000','#046C9A'))
  #Then create a variable that scales from red to blue between the two times
  cilla$Col <- rbPal(nrow(cilla))[as.numeric(cut(cilla$t,breaks = nrow(cilla)))]
  
  
  #Plot of the range estimate
  plot(cilla,
       UD=cilla.akde,
       col.grid = NA,
       family = "serif",
       pch = 20,
       cex = 0.2,
       col.DF = "#669543",
       col = cilla$Col,
       labels=FALSE)
  
  title(main = "d)", family = "serif",  adj = 0)
  
  

  dev.off()
  
  
  
  ##########################################################################################
  ##########################################################################################
  # Export all of the results
  ##########################################################################################
  ##########################################################################################
  
  #Get basic stats on the dataset
  res <- as.data.frame(BINOMIAL)
  res$ID <- cilla@info$identity
  res$Year <- median(as.numeric(format(as.Date(cilla$timestamp),"%Y")))
  res$Lat <- median(cilla$latitude)
  res$Long <- median(cilla$longitude)
  res$Frequency <- median(diff(cilla$t))/60
  res$Duration <- (tail(cilla$t,1) - head(cilla$t,1))/60/60/24
  res$n <- nrow(cilla)
  
    
  #Get tau_p
  if(any(grepl("position", row.names(summary(FIT, units = FALSE)$CI)))) {
    
  res$tau_p <- summary(FIT, units = FALSE)$CI[grepl("position", row.names(summary(FIT, units = FALSE)$CI)),2] 
  res$tau_p_min <- summary(FIT, units = FALSE)$CI[grepl("position", row.names(summary(FIT, units = FALSE)$CI)),1]
  res$tau_p_max <- summary(FIT, units = FALSE)$CI[grepl("position", row.names(summary(FIT, units = FALSE)$CI)),3]
  } else{res$tau_p <- NA 
    res$tau_p_min <- NA 
    res$tau_p_max <- NA
  }
  
  #Get tau_v
  if(any(grepl("velocity", row.names(summary(FIT, units = FALSE)$CI)))) {
    
    res$tau_v <- summary(FIT, units = FALSE)$CI[grepl("velocity", row.names(summary(FIT, units = FALSE)$CI)),2] 
    res$tau_v_min <- summary(FIT, units = FALSE)$CI[grepl("velocity", row.names(summary(FIT, units = FALSE)$CI)),1]
    res$tau_v_max <- summary(FIT, units = FALSE)$CI[grepl("velocity", row.names(summary(FIT, units = FALSE)$CI)),3]
  } else{res$tau_v <- NA 
  res$tau_v_min <- NA 
  res$tau_v_max <- NA
  }
  
  #Get diffusion estimates
  if("diffusion (square meters/second)" %in% row.names(summary(FIT, units = FALSE)$CI)) {
    res$diffusion <- summary(FIT, units = FALSE)$CI["diffusion (square meters/second)",2] 
    res$diffusion_min <- summary(FIT, units = FALSE)$CI["diffusion (square meters/second)",1] 
    res$diffusion_max <- summary(FIT, units = FALSE)$CI["diffusion (square meters/second)",3] 
  } else{res$diffusion <- NA 
      res$diffusion_min <- NA 
      res$diffusion_max <- NA
      }
  
  #Get HR results
  res$HR <- summary(cilla.akde, units = FALSE)$CI[2]
  res$HR_min <- summary(cilla.akde, units = FALSE)$CI[1]
  res$HR_max <- summary(cilla.akde, units = FALSE)$CI[3]
  
  #Get Speeds results
  if(speed == TRUE){
    res$Speed <- SPEED[2]
    res$Speed_min <- SPEED[1]
    res$Speed_max <- SPEED[3]
  }else{res$Speed <- NA
  res$Speed_min <- NA
  res$Speed_max <- NA
    
  }


  res

  }


