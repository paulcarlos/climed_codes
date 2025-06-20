######################################
# Function: Cells selector and tropical cyclone identifier (Dec 1999-Dec 2021) from ISIMIP3a dataset
# Author: Paul LC Chua 
# R packages: terra version 1.8-54
# Version: 2025-06-03


######################################
# REQUIRED INPUTS
# polygon = Path or SpatVector of polygon(s) or point(s) of locations or multiple locations
# name_location = Column name containing the name of the location or add vector of location names (No Defaults)
# file_source = Path of data repository for baseline rasters of tropical cyclones
# variable = "wind" or "rain", no default
# model = "H08" or "ER11", no default (H08 is Holland 2008; ER11 is Emanuel and Rotunno 2011)
# year_start = starting year for selection of tropical cyclone events, if NULL selects the lowest possible
# year_end =  end year for selection of tropical cyclone events, if NULL selects the highest possible
# min_cover = Minimum proportion of coverage of a cell/grid to consider (Default is no minimum)
# no_overlap = Assign a single location to a grid; TRUE or FALSE (Default is FALSE) selecting highest proportion of coverage
# time_zone = TRUE or FALSE, FALSE uses UTC 0:00, brings back Olson time zones
# tz_pop_centroid = 2000, 2005, 2010, 2015, 2020, or NULL, selects time zone of the most populated grid in the chosen year, NULL uses polygon/boundary centroid


######################################
# OUTPUT(S) OF THE FUNCTION
# List of dataframes of intersecting cells, population weights, and tropical cyclone IDs
# "cells" = dataframe containing:
#    "ID" = row number from input polygon(s)/point(s)
#    "locations" = name of the area/location/subnational unit
#    "cell" = cell numbers intersecting with polygon borders
#    "coverage" = proportion of cell/grid within/intersecting with polygon borders (values are 0-1.00, 1.00 means entire cell/grid within border)
#    "pop2000" = population count in year 2000
#    "pop2005" = population count in year 2005
#    "pop2010" = population count in year 2010
#    "pop2015" = population count in year 2015
#    "pop2020" = population count in year 2020
#    "tzone" = time zone for each location (not by grid cell)
# "trop_id" = dataframe containing:
#    "tcid" = Tropical cyclone ID numbers
#    "yrmo" = Year and month when the TC occurred
#    "basin" = basin/committee responsible for the first recorded event
# basins/committees are: NA=North Atlantic; EP=Eastern North Pacific; WP=Western North Pacific; NI=North Indian; SI=South Indian; SP=Southern Pacific; SA=South Atlantic
# attributes: "variable", "model", "min_cover", "time_zone", "file_source", "name_location" from REQUIRED INPUTS


######################################
# FUNCTION 
isimip_tc_id <- function(polygon, name_location, file_source, variable, model, year_start=NULL, year_end=NULL, min_cover=0, no_overlap=FALSE, time_zone=FALSE, tz_pop_centroid=NULL) {
  
  # check libraries
  reqlib <- c("terra","lutz")
  if (sum(reqlib%in%rownames(installed.packages()))!=length(reqlib)) {
    stop(paste0("install ",paste0(reqlib,collapse=" ")," R packages!"))
  }
  
  # check polygon
  if (class(polygon)[1]=="SpatVector") {
    poly1 <- terra::project(polygon,y=terra::crs("OGC:CRS84"))
  } else if (file.exists(polygon)) {
    poly1 <- terra::vect(polygon)
    poly1 <- terra::project(poly1,y=terra::crs("OGC:CRS84"))
  } else {
    stop("polygon not a SpatVector or a correct path/file name")
  }
  
  # check if file source exists
  if (substr(file_source,nchar(file_source),nchar(file_source))=="/") {
    file_source1 <- file_source
  } else {
    file_source1 <- paste0(file_source,"/")
  }
  file_source_ind <- paste0(file_source1,"individual/")
  file_source_pool <- paste0(file_source1,"pooled/")
  file_source_pop <- paste0(file_source1,"population/")
  if (sum(dir.exists(c(file_source_ind,file_source_pool,file_source_pop)))<3) {
    stop("file_source does not exist or incorrect!")
  }
  
  # validate if variable and model are correct
  if (model=="ER11") {
    modname <- "Emanuel and Rotunno 2011"
  } else {
    modname <- "Holland 2008"
  }
  varnm <- c("wind_er11","rain_er11","wind_h08","rain_h08")
  name_var <- paste0(tolower(variable),"_",tolower(model))
  if (!name_var%in%varnm) {
    stop("variable and model are incorrect!")
  }
  message("\nselected variable: ",tolower(variable),"\n","selected model: ",modname,"\n")
  
  # get files related to TCs
  reftab <- readRDS(paste0(file_source,"tc_dataframe_isimip3a_2000-2021.rds"))
  #files_pooled <- grep(".tif",list.files(file_source_pool,pattern=name_var),value=TRUE)
  files_pooled_ref <- readRDS(paste0(file_source_pool,name_var,"_list-tcid.rds"))
  
  # check years
  if (is.null(year_start)) {
    year_start <- min(as.integer(reftab$year))
  } else {
    if (year_start < min(as.integer(reftab$year))) {
      stop(paste0("year_start is lesser than the minimum year ",min(as.integer(reftab$year)),"!"))
    }
  }
  if (is.null(year_end)) {
    year_end <- max(as.integer(reftab$year))
  } else {
    if (year_end > max(as.integer(reftab$year))) {
      #stop(paste0("year_end is beyond maximum year ",min(as.integer(reftab$year)),"!"))
      year_end <- max(as.integer(reftab$year))
      message("year_end is replaced with the maximum year of ",year_end)
    } else if (year_end < min(as.integer(reftab$year))) {
      year_end <- min(as.integer(reftab$year))
      message("year_end is replaced with the minimum year of ",year_end)
    }
  }
  yrange <- year_start:year_end # create range of years
  
  # get intersecting cells
  lname <- unlist(poly1[[name_location]])
  ras1 <- terra::rast(paste0(file_source_pool,name_var,"_overall.tif"))
  celln <- terra::extract(ras1,poly1,method="simple",cells=TRUE,ID=TRUE,weights=TRUE,exact=TRUE)
  #celln <- data.frame(terra::cells(ras1,poly1,weights=T,exact=T,small=T))
  #celln$locations <- unname(lname[celln$ID])
  celln$locations <- lname[match(celln$ID,1:length(lname))]
  celln$overlap <- FALSE
  celln$overlap[duplicated(celln$cell)] <- TRUE
  #celln$tcval <- unlist(ras1[celln$cell])
  
  valnm <- setdiff(colnames(celln),c("ID","cell","weight","locations","overlap"))
  celln$tcval <- celln[,valnm]
  
  if (sum(celln$tcval,na.rm=TRUE)<1) { # STOP IF NO TC VALUE
    stop("polygon(s)/point(s) do not intersect with any tropical cylone!")
  }
  celln <- celln[celln$tcval>0 & !is.na(celln$tcval),]
  
  if (terra::geomtype(poly1)%in%c("polygon","polygons")) {
    if (min_cover>0) {
      celln1 <- celln[celln$weight>=min_cover,]
    } else {
      celln1 <- celln
    }
  } else if (terra::geomtype(poly1)=="points") {
    celln$weight <- 1
    celln1 <- celln
  }
  
  #celln2 <- aggregate(weights~ID+cell,data=celln1,"max",na.rm=T)
  if (no_overlap==TRUE) {
    celln2 <- aggregate(weight~ID+cell+locations,data=celln1,"max",na.rm=T)
  } else {
    celln2 <- celln1
  }
  celln2$coverage <- celln2$weight
  celln3 <- celln2[,c("ID","locations","cell","coverage")]
  
  # get population counts as weights
  for (p in c(2000,2005,2010,2015,2020)) {
    if (!file.exists(paste0(file_source_pop,"gpw411_pop-count_",p,"_isimip3a-tc.tif"))) {
      stop("population raster within file_source do not exist!")
    }
    popras <- terra::rast(paste0(file_source_pop,"gpw411_pop-count_",p,"_isimip3a-tc.tif"))
    celln3[,paste0("pop",p)] <- popras[celln3$cell][,1] 
    #celln3[is.na(celln3[,paste0("pop",p)]),paste0("pop",p)] <- 0
    celln3[is.na(celln3[,paste0("pop",p)]),paste0("pop",p)] <- 1 # avoid locations with all zeroes
  }
  
  # time zone finder
  if (time_zone==TRUE) {
    if (!is.null(tz_pop_centroid)) {
      if (tz_pop_centroid%in%c(2000,2005,2010,2015,2020)) {
        celln3$pcent <- celln3[,paste0("pop",tz_pop_centroid)] # pop centroid
        cent1 <- aggregate(pcent~ID+locations,data=celln3,FUN="max")
        cent1$cell <- celln3$cell[match(paste0(cent1$locations,"_",cent1$pcent),paste0(celln3$locations,"_",celln3$pcent))]
        cent1$lon <- terra::xFromCell(ras1,cent1$cell)
        cent1$lat <- terra::yFromCell(ras1,cent1$cell)
        cent1$tzone <- lutz::tz_lookup_coords(cent1$lat,cent1$lon,method="accurate",warn=FALSE)
        celln3$tzone <-cent1$tzone[match(celln3$locations,cent1$locations)]
        celln3$tzone <- gsub(";.*","",celln3$tzone) # eliminate multiple time zones from disputed areas
      } else {
        stop("tz_pop_centroid is not 2000, 2005, 2010, 2015, or 2020!")
      }
    } else {
      cent1 <- terra::centroids(poly1)
      cent2 <- as.data.frame(terra::crds(cent1))
      cent2$locations <- unlist(cent1[[name_location]])
      cent2$tzone <- lutz::tz_lookup_coords(cent2$y,cent2$x,method="accurate",warn=FALSE)
      celln3$tzone <- cent2$tzone[match(celln3$locations,cent2$locations)]
      celln3$tzone <- gsub(";.*","",celln3$tzone) # eliminate multiple time zones from disputed areas
    }
    # combine same time zones
    tzdf <- data.frame("tz"=unique(celln3$tzone),"offset"=NA)
    if (nrow(tzdf)>1) {
      tzdf$offset <- sapply(tzdf$tz,function(z)lutz::tz_offset(as.Date("2000-01-01",tz="UTC"),tz=z)[["utc_offset_h"]])
      ofst <- unique(tzdf$offset)
      if (length(ofst)>1) {
        tzdf$newtz <- NA
        for (ofs in ofst) {
          tzdf2 <- tzdf[tzdf$offset==ofs,]
          #tzdf$newtz[tzdf$offset==ofs] <- tzdf2$tz[1]
          celln3$tzone[celln3$tzone%in%tzdf2$tz] <- tzdf2$tz[1]
        }
      } else {
        celln3$tzone <- tzdf$tz[1]
      }
    }
    
  } else if (time_zone==FALSE) {
    celln3$tzone <- "UTC"
  } else if (time_zone%in%OlsonNames()) {
    celln3$tzone <- time_zone
  }
  
  #message("\nextracted intersecting cells, population weights, and time zones... \n \n")
  
  # loop intersecting rasters with polygon
  message("\nsearching TC basins: NI=North Indian, SI=South Indian, SP=Southern Pacific,\n",
          "WP=Western North Pacific, EP=Eastern North Pacific, NA=North Atlantic, SA=South Atlantic \n \n")
  tcdf <- data.frame(matrix(NA,nrow=0,ncol=length(reftab),dimnames=list(NULL,colnames(reftab))))
  for (f in seq(files_pooled_ref)) {
    basin <- names(files_pooled_ref)[f]
    idlist <- files_pooled_ref[[f]]
    for (g in seq(idlist)) {
      rgrp <- terra::rast(paste0(file_source_pool,name_var,"_",basin,"_grp",g,".tif"))
      if (sum(rgrp[celln2$cell],na.rm=TRUE)>0) {
        tcid <- idlist[[g]]
        if (any(tcid$year%in%yrange)) {
          tcid <- tcid[tcid$year%in%yrange,]
          for (h in 1:nrow(tcid)) {
            tcras <- terra::rast(paste0(file_source_ind,name_var,"_",basin,"_",tcid$tcid[h],".tif"))
            if (sum(tcras[celln2$cell],na.rm=TRUE)>0) {
              tcdf <- rbind(tcdf,tcid[h,])
            }
          }
        }
      }
    }
    message(paste0("searched ",basin,"...")," ")
  }
  if (nrow(tcdf)>0) {
    df1 <- tcdf[order(paste0(tcdf$year,"-",tcdf$month)),]
    rownames(df1) <- 1:nrow(df1)
  } else {
    df1 <- tcdf
  }
  
  
  # return list
  dlist <- list()
  dlist[["cells"]] <- celln3[,c("ID","locations","cell","coverage","pop2000","pop2005","pop2010","pop2015","pop2020","tzone")]
  dlist[["trop_id"]] <- df1
  attr(dlist,which="variable") <- variable
  attr(dlist,which="model") <- model
  attr(dlist,which="time_zone") <- time_zone
  attr(dlist,which="file_source") <- file_source_pool
  attr(dlist,which="name_location") <- name_location
  message("\n\n",nrow(df1)," tropical cyclones/storms were identified\n\n")
  return(dlist)
  
}

#rm(list=ls());gc()