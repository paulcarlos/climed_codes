######################################
# Function: Plotting selected cells and the accompanying polygon(s)/point(s) of administrative border for ERA5/ERA5-Land
# Author: Paul LC Chua 
# R package(s): terra version 1.8-5 
# Version: 2025-01-31


######################################
# REQUIRED INPUTS
# cellnum = dataframe of selected cells within ERA5 or ERA5-Land using cell selector
# polygon = Path or SpatVector of polygon of locations or multiple locations


######################################
# OUTPUT(S) OF THE FUNCTION
# Plot of grid(s) containing proportion of coverage within polygon of administrative border overlapped to the actual borders
# NOTE: The grids with color(s) are grids with value either in ERA5 or ERA5-Land


######################################
# FUNCTION
era5_cell_checker <- function(cellnum, polygon) {
  reqlib <- c("terra")
  if (sum(reqlib%in%rownames(installed.packages()))!=length(reqlib)) {
    stop(paste0("install ",paste0(reqlib,collapse=" ")," R packages!"))
  }
  
  # check if cell dataframe exists
  if (!class(cellnum)=="data.frame") {
    stop("cellnum input not data.frame!")
  }
  
  if (sum(colnames(cellnum)%in%c("locations","cell","coverage"))!=3) {
    stop("cellnum variables/column names incorrect!")
  }
  
  # check if data type is filled
  if (length(attr(cellnum,which="data_type"))>0) {
    data_type <- attr(cellnum,which="data_type")
  } else {
    stop("data_type does not exist in cellnum!")
  }
  if (length(data_type)==1) {
    if (tolower(data_type)!="era5-land" & tolower(data_type)!="era5") {
      stop("data_type is incorrect!")
    }
  } else {
    stop("data_type is incorrect!")
  }
  
  # upload sample raster
  file_source <- attr(cellnum,which="file_source")
  if (substr(file_source,nchar(file_source),nchar(file_source))=="/") {
    #file_source1 <- paste0(file_source,tolower(data_type),"/")
    file_source1 <- file_source
  } else {
    #file_source1 <- paste0(file_source,"/",tolower(data_type),"/")
    file_source1 <- paste0(file_source,"/")
  }
  if (!dir.exists(file_source1)) {
    stop("file_source does not exist!")
  }
  ras1 <- terra::rast(paste0(tolower(file_source1),"sample-raster_",tolower(data_type),".tif"))
  
  # create alternative empty raster 
  ras2 <- ras1
  #ras2[] <- NA
  ras2[!is.na(ras2)] <- -0.1
  
  # check polygon
  if (class(polygon)[1]=="SpatVector") {
    poly1 <- terra::project(polygon,y=terra::crs(ras1))
  } else if (file.exists(polygon)) {
    poly1 <- terra::vect(polygon)
    poly1 <- terra::project(poly1,y=terra::crs(ras1))
  } else {
    stop("polygon not a SpatVector or a correct path/file name")
  }
  name_location <- attr(cellnum,which="name_location")
  if (sum(names(poly1)==name_location)==0) {
      stop("name_location is incorrect or does not match with polygon!")
  }
  
  # get locations
  loc <- unique(cellnum$locations) # locations
  message("\nplotting proportion of cell coverages for ",length(loc)," locations... \n \n") # send notice/message
  
  # loop plots
  if (length(loc)>1) {
    par(mfrow=c(2,2),mar=c(3,5,4,2))
  } else {
    par(mar=c(3,5,4,2))
  }
  for (x in loc) {
    buf <- unique(cellnum$buffer[cellnum$locations==x])
    poly2 <- poly1[which(poly1[[name_location]]==x)] # get polygon
    
    if (terra::geomtype(poly2)%in%c("polygon","polygons")) {
      extnt <- terra::ext(terra::buffer(poly2,width=buf)) # extent of plot
    } else if (terra::geomtype(poly2)%in%c("point","points")) {
      extnt <- terra::ext(poly2)
    }
    
    extnt[c(1,3)] <- floor(extnt[c(1,3)])
    extnt[c(2,4)] <- ceiling(extnt[c(2,4)])
    #ras2[cellnum$cell] <- -0.1
    ras3 <- ras2
    ras3[cellnum$cell[cellnum$locations==x]] <- cellnum$coverage[cellnum$locations==x]
    ras4 <- terra::crop(ras3,extnt)
    terra::plot(ras4,col=c("#B2BEB5","#219EBC","#023047","#FFB703","#FB8500"),
                breaks=c(-0.1,0,0.25,0.5,0.75,1.1),
                main=paste0(x," \n(proportion of cell coverage)"),cex.main=0.9)
    
    if (terra::geomtype(poly2)%in%c("polygon","polygons")) {
      terra::lines(poly2,col="red",lwd=3)
    } else if (terra::geomtype(poly2)%in%c("point","points")) {
      terra::points(poly2,col="red")
    }
    
  }
  message("\ngrey colored cells have ",data_type," values and white colored cells are empty cells\n")
}

#rm(list=ls());gc()