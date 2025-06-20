######################################
# Function: Select ERA5/ERA5-Land cells/grids within/intersecting administrative border (polygon/point) of a location or multiple locations
# Author: Paul LC Chua 
# R package(s): terra version 1.8-54, lutz version 0.3.2 
# Version: 2025-06-03


######################################
# REQUIRED INPUTS
# polygon = Path or SpatVector of polygon(s) or point(s) of locations or multiple locations 
# name_location = Column name containing the name of the location or add vector of location names (No Defaults)
# data_type = "ERA5" or "ERA5-Land" (No Defaults)
# file_source = Path of data repository for base rasters
# min_cover = Minimum proportion of coverage of a cell/grid to consider (Default is no minimum)
# no_overlap = Assign a single location to a grid; TRUE or FALSE (Default is FALSE) selecting highest proportion of coverage
# buffer_dist = For ERA5-Land, distance in meters to widen the search for most nearest grids with values because some small islands are NA
# time_zone = TRUE or FALSE, FALSE uses UTC 0:00, brings back Olson time zones
# tz_pop_centroid = 2000, 2005, 2010, 2015, 2020, or NULL, selects time zone of the most populated grid in the chosen year, NULL uses polygon/boundary centroid


######################################
# OUTPUT(S) OF THE FUNCTION
# DATAFRAME OF SEVERAL VARIABLES:
# "ID" = row number from input polygon(s)/point(s)
# "locations" = name of the area/location/subnational unit
# "cell" = cell numbers intersecting with polygon borders
# "coverage" = proportion of cell/grid within/intersecting with polygon borders (values are 0-1.00, 1.00 means entire cell/grid within border)
# "pop2000" = population count in year 2000
# "pop2005" = population count in year 2005
# "pop2010" = population count in year 2010
# "pop2015" = population count in year 2015
# "pop2020" = population count in year 2020
# "tzone" = time zone
# attributes: "data_type","file_source","name_location", "time_zone" (to carry over to plot and extraction functions)  
# NOTE: population weights derived from Gridded Population of the World, Version 4 (GPWv4): Population Count Adjusted to Match 2015 Revision of UN WPP Country Totals, Revision 11)



######################################
# FUNCTION 
era5_cell_selector <- function(polygon, name_location, data_type, file_source, min_cover=0, no_overlap=FALSE, buffer_dist=0, time_zone=FALSE, tz_pop_centroid=NULL) {
  
  # check R packages
  reqlib <- c("terra","lutz")
  if (sum(reqlib%in%rownames(installed.packages()))!=length(reqlib)) {
    stop(paste0("install ",paste0(reqlib,collapse=" ")," R packages!"))
  }
  
  # check if data type is filled
  if (length(data_type)==1) {
    if (tolower(data_type)!="era5-land" & tolower(data_type)!="era5") {
      stop("data_type is incorrect!")
    }
  } else {
    stop("data_type is incorrect!")
  }
  
  # upload sample raster
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
  
  # check polygon
  if (class(polygon)[1]=="SpatVector") {
    poly1 <- terra::project(polygon,y=terra::crs(ras1))
  } else if (file.exists(polygon)) {
    poly1 <- terra::vect(polygon)
    poly1 <- terra::project(poly1,y=terra::crs(ras1))
  } else {
    stop("polygon not a SpatVector or a correct path/file name")
  }
  if (sum(names(poly1)==name_location)==0) {
    stop("name_location is incorrect!")
  }
  
  # check min_cover
  if (min_cover<0 | min_cover>1) {
    stop("min_cover beyond 0.00 to 1.00!")
  }
  
  # get cell numbers via extract function
  lname <- unlist(poly1[[name_location]])
  celln <- terra::extract(ras1,poly1,method="simple",cells=TRUE,ID=TRUE,weights=TRUE,exact=TRUE)
  if (terra::geomtype(poly1)%in%c("polygon","polygons")) {
    if (!is.null(min_cover)) {
      celln1 <- celln[celln$weight>=min_cover,]
    } else {
      celln1 <- celln
    }
  } else if (terra::geomtype(poly1)=="points") {
    celln$weight <- 1
    celln1 <- celln
  }
  celln1$locations <- lname[match(celln1$ID,1:length(lname))]
  celln1$overlap <- FALSE
  celln1$overlap[duplicated(celln1$cell)] <- TRUE
  
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
  
  # for all NAs, replace with nearest grids within buffer zone
  na_grids <- aggregate(na_value~locations,data=celln3,FUN=all)
  na_loc <- na_grids$locations[na_grids$na_value==TRUE]
  if (tolower(data_type)=="era5-land" & buffer_dist>0) {
    for (n in na_loc) {
      poly2 <- poly1[unlist(poly1[[name_location]])==n]
      pbuff <- terra::buffer(poly2,width=buffer_dist)
      cellbuff <- terra::extract(ras1,pbuff,method="simple",cells=TRUE,ID=TRUE,weights=TRUE,exact=TRUE)
      cellbuff$values <- cellbuff[,valnm]
      cellbuff$na_value <- is.na(cellbuff$values)
      if (all(cellbuff$na_value)) {
        #stop("buffer distance results to NAs, increase the buffer distance")
        message("buffer distance results to NAs, increase the buffer distance")
      } else {
        cellbuff$ID <- unique(celln3$ID[celln3$locations==n])
        cellbuff2 <- cellbuff[cellbuff$na_value==FALSE,]
        cellbuff2$weight <- 0
        cellbuff2$locations <- n
        cellbuff2 <- cellbuff2[,colnames(celln3)]
        celln3 <- rbind(celln3,cellbuff2)
      }
      
    }
    celln3$buffer <- 0
    celln3$buffer[celln3$locations%in%na_loc] <- buffer_dist
    message(paste0(length(na_loc)," location(s)/area(s) used cell value(s) from nearest cell(s)/grid(s) within ",buffer_dist," m"))
  } else if (tolower(data_type)=="era5" & buffer_dist>0) {
    message("ERA5 does not need buffer_dist")
  } else if (buffer_dist==0) {
    celln3$buffer <- 0
  }
  
  # remove areas without value
  if (length(na_loc)>0 & tolower(data_type)=="era5-land" & buffer_dist==0) {
    message(paste0(length(na_loc)," location(s)/area(s) was removed because no cell(s) have value(s)"))
  }
  celln4 <- celln3[celln3$na_value==FALSE,]
  
  # get population counts as weights
  for (i in c(2000,2005,2010,2015,2020)) {
    if (!file.exists(paste0(file_source1,"gpw411_pop-count_",i,"_",tolower(data_type),".tif"))) {
      stop("population raster within file_source do not exist!")
    }
    popras <- terra::rast(paste0(file_source1,"gpw411_pop-count_",i,"_",tolower(data_type),".tif"))
    celln4[,paste0("pop",i)] <- popras[celln4$cell][,1] 
  }
  
  # remove column
  cnm_rm <- !grepl("na_value",colnames(celln4))
  celln4 <- celln4[,cnm_rm]
  renam <- grep("weight",colnames(celln4))
  colnames(celln4)[renam] <- "coverage"
  
  # time zone finder
  if (time_zone==TRUE) {
    if (tz_pop_centroid%in%c(2000,2005,2010,2015,2020)) {
      celln4$pcent <- celln4[,paste0("pop",tz_pop_centroid)]
      cent1 <- aggregate(pcent~ID+locations,data=celln4,FUN="max")
      cent1$cell <- celln4$cell[match(paste0(cent1$locations,"_",cent1$pcent),paste0(celln4$locations,"_",celln4$pcent))]
      cent1$lon <- terra::xFromCell(ras1,cent1$cell)
      cent1$lat <- terra::yFromCell(ras1,cent1$cell)
      cent1$tzone <- lutz::tz_lookup_coords(cent1$lat,cent1$lon,method="accurate",warn=FALSE)
      celln4$tzone <- cent1$tzone[match(celln4$locations,cent1$locations)]
      celln4$tzone <- gsub(";.*","",celln4$tzone) # eliminate multiple time zones from disputed areas
    } else {
      cent1 <- terra::centroids(poly1)
      cent2 <- as.data.frame(terra::crds(cent1))
      cent2$locations <- unlist(cent1[[name_location]])
      cent2$tzone <- lutz::tz_lookup_coords(cent2$y,cent2$x,method="accurate",warn=FALSE)
      celln4$tzone <- cent2$tzone[match(celln4$locations,cent2$locations)]
      celln4$tzone <- gsub(";.*","",celln4$tzone) # eliminate multiple time zones from disputed areas
    }
    # combine same time zones
    tzdf <- data.frame("tz"=unique(celln4$tzone),"offset"=NA)
    if (nrow(tzdf)>1) {
      tzdf$offset <- sapply(tzdf$tz,function(z)lutz::tz_offset(as.Date("2000-01-01",tz="UTC"),tz=z)[["utc_offset_h"]])
      ofst <- unique(tzdf$offset)
      if (length(ofst)>1) {
        tzdf$newtz <- NA
        for (ofs in ofst) {
          tzdf2 <- tzdf[tzdf$offset==ofs,]
          #tzdf$newtz[tzdf$offset==ofs] <- tzdf2$tz[1]
          celln4$tzone[celln4$tzone%in%tzdf2$tz] <- tzdf2$tz[1]
        }
      } else {
        celln4$tzone <- tzdf$tz[1]
      }
    }
    
  } else if (time_zone==FALSE) {
    celln4$tzone <- "UTC"
  } else if (time_zone%in%OlsonNames()) {
    celln4$tzone <- time_zone
  } else {
    stop("time_zone is incorrect!")
  }
  
  
  # final dataframe
  celln5 <- celln4[,c("ID","locations","cell","coverage","buffer","pop2000","pop2005","pop2010","pop2015","pop2020","tzone")]
  celln5 <- celln5[order(celln5$ID),]
  rownames(celln5) <- 1:nrow(celln5)
  attr(celln5,which="data_type") <- data_type
  attr(celln5,which="file_source") <- file_source
  attr(celln5,which="name_location") <- name_location
  attr(celln5,which="time_zone") <- time_zone
  return(celln5)
}

#rm(list=ls());gc()
