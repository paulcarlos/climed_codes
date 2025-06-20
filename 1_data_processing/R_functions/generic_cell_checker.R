######################################
# Function: Plotting selected cells and the accompanying polygon(s)/point(s) of administrative border 
# Author: Paul LC Chua 
# R package(s): terra version 1.8-5 
# Version: 2025-01-31


######################################
# REQUIRED INPUTS
# cellnum = dataframe of selected cells from generic_cell_extractor function
# polygon = Path or SpatVector of polygon of locations or multiple locations
# raster = File path or SpatRaster of the gridded data to be extracted (for single variable rasters only)


######################################
# OUTPUT(S) OF THE FUNCTION
# Plot of grid(s) containing proportion of coverage within polygon of administrative border overlapped to the actual borders
# NOTE: The grids with color(s) are grids with values or non NAs


######################################
# FUNCTION
generic_cell_checker <- function(cellnum, polygon, raster) {
  
  # check R packages
  reqlib <- c("terra")
  if (sum(reqlib%in%rownames(installed.packages()))!=length(reqlib)) {
    stop(paste0("install ",paste0(reqlib,collapse=" ")," R packages!"))
  }
  
  # check polygon
  if (class(polygon)[1]=="SpatVector") {
    poly1 <- terra::project(polygon,y="OGC:CRS84")
  } else if (file.exists(polygon)) {
    poly1 <- terra::vect(polygon)
    poly1 <- terra::project(poly1,y="OGC:CRS84")
  } else {
    stop("polygon not a SpatVector or a correct path/file name")
  }
  name_location <- attr(cellnum,which="name_location")
  if (sum(names(poly1)==name_location)==0) {
    stop("name_location is incorrect!")
  }
  
  # check raster
  if (class(raster)[1]=="SpatRaster") {
    if (terra::crs(raster)=="OGC:CRS84") {
      ras1 <- raster
    } else {
      ras1 <- terra::project(raster,y="OGC:CRS84")
    }
  } else if (file.exists(raster)) {
    ras1 <- terra::rast(raster)
    if (terra::crs(ras1)=="OGC:CRS84") {
      ras1 <- terra::project(ras1,y="OGC:CRS84")
    }
  } else {
    stop("raster not a SpatRaster or a correct path/file name")
  }
  
  # get min cover
  min_cover <- attr(cellnum,which="min_cover")
  message("\nselected minimum grid/cell coverage is ",min_cover,"\n\n")
  
  # get locations
  loc <- unique(cellnum$locations) # locations
  message("\nplotting proportion of cell coverages for ",length(loc)," locations... \n \n") # send notice/message
  
  # create alternative empty raster 
  ras2 <- terra::rast(terra::ext(ras1),crs="OGC:CRS84", resolution=terra::res(ras1))
  ras2[] <- NA
  ras2[which(!is.na(ras1[]))] <- -0.1
  
  # loop plots
  if (length(loc)>1) {
    par(mfrow=c(2,2),mar=c(3,5,4,2))
  } else {
    par(mar=c(3,5,4,2))
  }
  for (x in loc) {
    poly2 <- poly1[which(poly1[[name_location]]==x)] # get polygon
    extnt <- terra::ext(poly2) # extent of plot
    extnt[c(1,3)] <- floor(extnt[c(1,3)])
    extnt[c(2,4)] <- ceiling(extnt[c(2,4)])
    #ras2[cellnum$cell] <- -0.1
    ras3 <- ras2
    ras3[cellnum$cell[cellnum$locations==x]] <- cellnum$coverage[cellnum$locations==x]
    ras4 <- terra::crop(ras3,extnt)
    terra::plot(ras4,col=c("#B2BEB5","#219EBC","#023047","#FFB703","#FB8500"),
                breaks=c(-0.1,0,0.25,0.5,0.75,1.1),
                main=paste0(x," \n(proportion of cell coverage)"),cex.main=0.9)
    terra::lines(poly2,col="red",lwd=3)
  }
  message("\ngrey colored cells have values and white colored cells are empty cells\n")
  
}


#rm(list=ls());gc()
