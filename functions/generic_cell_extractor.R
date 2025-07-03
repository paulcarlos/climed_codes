######################################
# Function: Select cells/grids and extract values within/intersecting administrative border (polygon/point) of a location or multiple locations
# Author: Paul LC Chua 
# R package(s): terra version 1.8-54 
# Version: 2025-06-03


######################################
# REQUIRED INPUTS
# polygon = File path or SpatVector of polygon(s) or point(s) of locations or multiple locations 
# name_location = Column name containing the name of the location or add vector of location names (No Defaults)
# raster = File path or SpatRaster of the gridded data to be extracted (for single variable rasters only)
# min_cover = Minimum proportion of coverage of a cell/grid to consider (Default is no minimum)
# no_overlap = Assign a single location to a grid; TRUE or FALSE (Default is FALSE) selecting highest proportion of coverage
# buffer_dist = Distance based on chosen CRS to widen the search for most nearest grids with values, if default "OGC:CRS84", distance is in meters 
# NOTE: CRS used is "OGC:CRS84"


######################################
# OUTPUT(S) OF THE FUNCTION
# DATAFRAME OF SEVERAL VARIABLES:
# "ID" = row number from input polygon(s)/point(s)
# "locations" = name of the area/location/subnational unit
# "cell" = cell numbers intersecting with polygon borders
# "coverage" = proportion of cell/grid within/intersecting with polygon borders (values are 0-1.00, 1.00 means entire cell/grid within border)
# attributes: "name_location"


######################################
# FUNCTION 
generic_cell_extractor <- function(polygon, name_location, raster, min_cover=0, no_overlap=FALSE, buffer_dist=0) {
  
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
  if (sum(names(poly1)==name_location)==0) {
    stop("name_location is incorrect!")
  }
  
  # check raster
  if (class(raster)[1]=="SpatRaster") {
    if (terra::crs(raster)==terra::crs("OGC:CRS84")) {
      ras1 <- raster
    } else {
      ras1 <- terra::project(raster,y="OGC:CRS84")
    }
  } else if (file.exists(raster)) {
    ras1 <- terra::rast(raster)
    if (terra::crs(ras1)==terra::crs("OGC:CRS84")) {
      ras1 <- terra::project(ras1,y="OGC:CRS84")
    }
  } else {
    stop("raster not a SpatRaster or a correct path/file name")
  }
  
  # check min_cover
  if (min_cover<0 | min_cover>1) {
    stop("min_cover beyond 0.00 to 1.00!")
  }
  
  # extract values numbers
  lname <- unlist(poly1[[name_location]])
  celln <- terra::extract(ras1,poly1,method="simple",cells=TRUE,ID=TRUE,weights=TRUE,exact=TRUE)
  #celln <- data.frame(terra::cells(ras1,poly1,method="simple",weights=TRUE,exact=TRUE,small=FALSE))
  if (terra::geomtype(poly1)%in%c("polygon","polygons")) {
    if (!is.null(min_cover)) {
      celln1 <- celln[celln$weight>=min_cover,]
    } else {
      celln1 <- celln
    }
  } else if (terra::geomtype(poly1)=="points") {
    celln$weights <- 1
    celln1 <- celln
  }
  celln1$locations <- lname[match(celln1$ID,1:length(lname))]
  celln1$overlap <- FALSE
  celln1$overlap[duplicated(celln1$cell)] <- TRUE
  #sum(celln$overlap)
  valnm <- setdiff(colnames(celln),c("ID","cell","weight","locations","overlap"))
  celln1$values <- celln1[,valnm]
  
  if (no_overlap==TRUE) {
    celln0 <- aggregate(weight~ID+cell+locations,data=celln1,"max",na.rm=T)
    celln2 <- celln1[celln1$cell%in%celln0$cell,]
  } else {
    celln2 <- celln1
  }
  celln3 <- celln2[,c("ID","locations","cell","weight","values")]
  
  # NA values
  celln3$na_value <- is.na(celln3$values)
  # value within raster
  #celln3$values <- unlist(ras1[celln3$cell])
  #celln3$na_value <- is.na(celln3$values)
  
  # for all NAs, replace with nearest grids within buffer zone
  na_grids <- aggregate(na_value~locations,data=celln3,FUN=all)
  na_loc <- na_grids$locations[na_grids$na_value==TRUE]
  for (n in na_loc) {
    message(paste0(length(na_loc)," location(s)/area(s) used cell value(s) from nearest cell(s)/grid(s) within ",buffer_dist," m"))
    poly2 <- poly1[unlist(poly1[[name_location]])==n]
    pbuff <- terra::buffer(poly2,width=buffer_dist)
    #cellbuff <- data.frame(terra::cells(ras1,pbuff,weights=T,exact=T,small=T))
    #cellbuff$values <- unlist(ras1[cellbuff$cell])
    cellbuff <- terra::extract(ras1,pbuff,method="simple",cells=TRUE,ID=TRUE,weights=TRUE,exact=TRUE)
    cellbuff$values <- cellbuff[,valnm]
    cellbuff$na_value <- is.na(cellbuff$values)
    if (all(cellbuff$na_value)) {
      #stop("buffer distance results to NAs, increase the buffer distance")
      message("buffer distance results to NAs, increase the buffer distance")
    } else {
      cellbuff2 <- cellbuff[cellbuff$na_value==FALSE,]
      cellbuff2$weight <- 0
      cellbuff2$locations <- n
      cellbuff2 <- cellbuff2[,colnames(celln3)]
      celln3 <- rbind(celln3,cellbuff2)
    }
    
  }
  
  # remove areas without value
  if (length(na_loc)>0 & buffer_dist==0) {
    message(paste0(length(na_loc)," location(s)/area(s) was removed because no cell(s) have value(s)"))
  }
  celln4 <- celln3[celln3$na_value==F,]
  
  # remove column
  cnm_rm <- !grepl("na_value",colnames(celln4))
  celln4 <- celln4[,cnm_rm]
  renam <- grep("weight",colnames(celln4))
  colnames(celln4)[renam] <- "coverage"
  
  # final dataframe
  celln5 <- celln4[,c("ID","locations","cell","coverage","values")]
  attr(celln5,which="name_location") <- name_location
  attr(celln5,which="min_cover") <- min_cover
  attr(celln5,which="no_overlap") <- no_overlap
  return(celln5)
}

#rm(list=ls());gc()