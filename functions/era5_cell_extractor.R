######################################
# Function: Extracting ERA5/ERA5-Land cell values from 2000-2023 (minimum processing is by month) 
# Author: Paul LC Chua 
# Year: 2025
# R packages: terra version 1.8-5; lubridate version 1.9.4; lutz version 0.3.2 
# Version: 2025-01-31


######################################
# REQUIRED INPUTS
# cellnum = dataframe of selected cells within ERA5 or ERA5-Land using cell selector
# file_source = Path of data repository for hourly ERA5 TIFs
# start_time = YYYY/MM, only from 2000-01 to 2023-12 (ISO 2000 wk1 starts Jan 3; Epi 2000 wk1 starts Jan 2)  
# end_time = YYYY/MM, only from 2000-01 to 2023-12 (ISO 2023 wk52 ends Dec 31; Epi 2023 wk52 ends Dec 30)
# variable = "Temperature" or "Precipitation"
# pop_weights = TRUE or FALSE, TRUE uses population weights using GPWv4 
# time_res = "Daily" or "Weekly"
# format_week = If weekly, "ISO" or "EPI", ISO is ISO 8601 while EPI is Epiweek of US-CDC


######################################
# OUTPUT(S) OF THE FUNCTION
# Dataframe of time series of extracted variable
# "location" = name of location(s)
# "date" or "yrwk" = date (YYYY-MM-DD) or year-week (YYYY-WW)
# "tp" or "t2m" = "tp" is total precipitation and "t2m" is temperature


######################################
# FUNCTION
era5_cell_extractor <- function(cellnum, file_source, start_time, end_time, variable, pop_weights=FALSE, time_res, format_weekly=NULL) {
  # check if libraries are installed
  reqlib <- c("terra","lubridate","lutz")
  if (sum(reqlib%in%rownames(installed.packages()))!=length(reqlib)) {
    stop(paste0("install ",paste0(reqlib,collapse=" ")," R packages!"))
  }
  
  # check if cell dataframe exists
  if (!class(cellnum)=="data.frame") {
    stop("cellnum input not data.frame!")
  }
  if (pop_weights==TRUE) {
    if (sum(colnames(cellnum)%in%c("locations","cell","pop2000","pop2005","pop2010","pop2015","pop2020"))!=7) {
      stop("cellnum does not contain population variables!")
    }
  } else if (time_zone==TRUE)  {
    if (sum(colnames(cellnum)%in%c("locations","cell","tzone"))!=3) {
      stop("cellnum does not contain time zone(s)!")
    }
  } else {
    if (sum(colnames(cellnum)%in%c("locations","cell"))!=2) {
      stop("cellnum variables/column names incorrect!")
    }
  }
  
  # check data_type and file source
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
  if (substr(file_source,nchar(file_source),nchar(file_source))=="/") {
    #file_source1 <- paste0(file_source,"tif_",tolower(data_type),"/")
    file_source1 <- paste0(file_source,"rast_",tolower(data_type),"/")
  } else {
    #file_source1 <- paste0(file_source,"/tif_",tolower(data_type),"/")
    file_source1 <- paste0(file_source,"/rast_",tolower(data_type),"/")
  }
  if (tolower(variable)=="precipitation") {
    vnm <- "tp"
    message("total precipitation (tp) is expressed in millimeters")
  } else if (tolower(variable)=="temperature") {
    vnm <- "t2m"
    message("2-meter temperature (t2m) is expressed in Celsius")
  }
  file_source2 <- paste0(file_source1,vnm,"/")
  if (!dir.exists(file_source2)) {
    stop("file_source does not exist!")
  }
  
  # check if time is available
  minyear <- 2000 # UPDATE THE YEAR IF NEEDED
  maxyear <- 2023 # UPDATE THE YEAR IF NEEDED
  full_time <- sort(do.call(paste0,expand.grid(minyear:maxyear,"-",sprintf("%02d",1:12))))
  if (which(full_time==start_time)>which(full_time==end_time)) {
    stop("start_time or end_time is incorrect!")
  }
  if (all(full_time!=start_time) & all(full_time!=end_time)) {
    stop("start_time and/or end_time is not between 2000-01 to 2023-12")
  }
  if (tolower(time_res)=="weekly") {
    if (is.null(format_weekly)) {
      stop("if time resolution is Weekly, indicate format_weekly either ISO or Epi!")
    }
    if (!tolower(format_weekly)%in%c("iso","epi")) {
      stop("incorrect format_weekly! please put either ISO or Epi")
    }
  }
  full_time2 <- c(paste0(minyear-1,"-12"),full_time,paste0(maxyear+1,"-01"))
  
  # check if time zone available (lowest is UTC-12, highest is UTC+14)
  time_zone <- attr(cellnum,which="time_zone")
  if (time_zone==TRUE | time_zone%in%OlsonNames()) {
    tz1 <- unique(cellnum$tzone)
    if (!all(tz1%in%OlsonNames())) {
      stop("time zone is incorrect or not in Olson format!")
    }
    tz2 <- data.frame("tz"=tz1,"dif"=unname(sapply(tz1,function(z) lutz::tz_offset("2000-12-20",tz=z)$utc_offset_h)))
    if (length(tz2$dif[tz2$dif<0])>0) { 
      hrneg <- min(tz2$dif[tz2$dif<0],na.rm=TRUE)
    } else {
      hrneg <- 0
    }
    if (length(tz2$dif[tz2$dif>0])>0) {
      hrpos <- max(tz2$dif[tz2$dif>0],na.rm=TRUE)
    } else {
      hrpos <- 0
    }
  } else if (time_zone==FALSE) {
    tz1 <- "UTC"
    hrneg <- hrpos <- 0
  } 
  
  # year-month sequence
  time_seq <- full_time[which(full_time==start_time):which(full_time==end_time)]
  
  # loop extraction
  colnam <- c("locations","date",vnm) # all extracted by daily
  dfsav <- as.data.frame(matrix(NA,nrow=0,ncol=length(colnam),dimnames=list(NULL,colnam))) # blank dataframe to save all values
  message("\nextracting month by month... (format: YYYY-MM)\n")
  
  for (t in seq(time_seq)) {
    
    ymon1 <- time_seq[t] # year-month
    message(ymon1," ")
    
    #if (t==1) {
    #  extval1 <- cellnum[,c("locations","cell")]
    #} else {
    #  extval1 <- extval1[,c(1:2,(length(extval1)-47):length(extval1))] # retain latest 48 hours 
    #}
    extval1 <- cellnum[,c("locations","cell")] # starting dataframe
    
    if (hrneg>0 & hrpos>0) {
      ym1 <- full_time2[which(full_time2==start_time)-1]
      ym3 <- full_time2[which(full_time2==end_time)+1]
      ymras1 <- terra::rast(list.files(paste0(file_source2,substr(ym1,1,4),"/"),pattern=paste0(ym1,"|",gsub("-","_",ym1)),full.names=TRUE))
      ymras2 <- terra::rast(list.files(paste0(file_source2,substr(ymon1,1,4),"/"),pattern=paste0(ymon1,"|",gsub("-","_",ymon1)),full.names=TRUE))
      ymras3 <- terra::rast(list.files(paste0(file_source2,substr(ym3,1,4),"/"),pattern=paste0(ym3,"|",gsub("-","_",ym3)),full.names=TRUE))
      ymras4 <- c(ymras1,ymras2,ymras3)
      hseqneg <- terra::time(ymras1)
      hseq1 <- hseqneg[(length(hseqneg)-(hrneg-1)):length(hseqneg)]
      hseqpos <- terra::time(ymras3)
      hseq3 <- hseqpos[1:hrpos]
      hseqfull <- c(terra::time(ymras1),terra::time(ymras2),terra::time(ymras3))
      hseq <- c(hseq1,terra::time(ymras2),hseq3)
      ymon_rast <- ymras4[[which(hseqfull%in%hseq)]]
    } else if (hrneg>0 & hrpos==0) {
      ym1 <- full_time2[which(full_time2==start_time)-1]
      ymras1 <- terra::rast(list.files(paste0(file_source2,substr(ym1,1,4),"/"),pattern=paste0(ym1,"|",gsub("-","_",ym1)),full.names=TRUE))
      ymras2 <- terra::rast(list.files(paste0(file_source2,substr(ymon1,1,4),"/"),pattern=paste0(ymon1,"|",gsub("-","_",ymon1)),full.names=TRUE))
      ymras4 <- c(ymras1,ymras2)
      hseqneg <- terra::time(ymras1)
      hseq1 <- hseqneg[(length(hseqneg)-(hrneg-1)):length(hseqneg)]
      hseqfull <- c(terra::time(ymras1),terra::time(ymras2))
      hseq <- c(hseq1,terra::time(ymras2))
      ymon_rast <- ymras4[[which(hseqfull%in%hseq)]]
    } else if (hrneg==0 & hrpos>0) {
      ym3 <- full_time2[which(full_time2==end_time)+1]
      ymras2 <- terra::rast(list.files(paste0(file_source2,substr(ymon1,1,4),"/"),pattern=paste0(ymon1,"|",gsub("-","_",ymon1)),full.names=TRUE))
      ymras3 <- terra::rast(list.files(paste0(file_source2,substr(ym3,1,4),"/"),pattern=paste0(ym3,"|",gsub("-","_",ym3)),full.names=TRUE))
      ymras4 <- c(ymras2,ymras3)
      hseqpos <- terra::time(ymras3)
      hseq3 <- hseqpos[1:hrpos]
      hseqfull <- c(terra::time(ymras2),terra::time(ymras3))
      hseq <- c(terra::time(ymras2),hseq3)
      ymon_rast <- ymras4[[which(hseqfull%in%hseq)]]
    }
    
    #ymon_rast <- terra::rast(list.files(paste0(file_source2,substr(ymon1,1,4),"/"),pattern=paste0(ymon1,"|",gsub("-","_",ymon1)),full.names=TRUE))
    #hseq <- terra::time(ymon_rast)
    #dmon <- unname(lubridate::days_in_month(as.Date(paste0(ymon1,"-01","%Y-%m-%d"))))
    #hseq <- seq(as.POSIXct(paste0(ymon1,"-01 00:00"),tz="UTC")-lubridate::days(2),as.POSIXct(paste0(ymon1,"-",dmon," 23:00"),tz="UTC")+lubridate::days(2),"hour")
    
    for (h in 1:length(hseq)) {
      hrval <- hseq[h]
      hnam <- paste0("H",gsub("-","",as.character(as.Date(hrval))),"_",sprintf("%02d",lubridate::hour(hrval)),"H")
      #if (!any(hnam==colnames(extval1))) {
      #  extrastr <- terra::rast(paste0(file_source2,lubridate::year(hrval),"/",sprintf("%02d",lubridate::month(hrval)),"/",
      #                                 tolower(data_type),"_",vnm,"_",as.character(as.Date(hrval)),"_",sprintf("%02d",lubridate::hour(hrval)),"H.tif"))
      #  extval1[,hnam] <- extrastr[extval1$cell]
      #}
      extrastr <- ymon_rast[[h]]
      extval1[,hnam] <- extrastr[extval1$cell]
      #cat(h," ")
    }
    
    # aggregate values
    if (tolower(data_type)=="era5-land" & vnm=="tp") { # accumulation from Day 1 hour 0:00 to Day 2 hour 0:00 or 25 hours 
      
      mat1 <- as.matrix(unname(extval1[,3:(length(hseq)+2)]))  # create matrix of hourly values
      mat2 <- matrix(NA,nrow=nrow(mat1),ncol=ncol(mat1)) # matrix of 1 hour delay
      mat2[,2:dim(mat1)[2]] <- mat1[,2:dim(mat1)[2]] - mat1[,1:(dim(mat1)[2]-1)]
      mat2[,which(lubridate::hour(hseq)==1)] <- mat1[,which(lubridate::hour(hseq)==1)] - mat2[,which(lubridate::hour(hseq)==0)]
      
      for (b in tz1) { # loop by time zone groups
        rntz <- which(cellnum$tzone==b) # row numbers with same time zones
        loc1 <- split(1:length(rntz),cellnum$locations[rntz]) # split row numbers by location
        
        if (all(unlist(lapply(loc1,length))==1)) {
          loc3 <- loc1 # locations are all single grids
        } else if (any(unlist(lapply(loc1,length))==1)) { # split further if there are subnational locations with a single grid
          #sploc <- which(unlist(lapply(loc1,length))==1)
          loc2 <- loc1[unlist(lapply(loc1,length))!=1] # multiple grids
          loc3 <- loc1[unlist(lapply(loc1,length))==1] # single grids
        }
        
        hdif <- tz2$dif[tz2$tz==b] # UTC difference hour
        daytz <- as.Date(lubridate::with_tz(hseq,tzone=b))
        mat3 <- t(mat2[rntz,]) # transpose so x are hours and y are grids
        mat4 <- rowsum(mat3,group=daytz,na.rm=TRUE) # sum of precipitation per day by grid
        mat3 <- t(mat4) # transpose back so x are grids and y are days
        if (pop_weights==TRUE) {
          yrpop <- seq(2000,2020,5)
          yrdif <- as.numeric(substr(ymon1,1,4))-yrpop
          if (any(yrdif%in%-4:0)) {
            yrval <- paste0("pop",yrpop[which(yrdif%in%-4:0)])
          } else {
            yrval <- paste0("pop",yrpop[which(yrdif==min(yrdif))])  # years beyond 2020 uses 2020
          }
          popval <- unlist(split(cellnum[rntz,yrval],cellnum$locations[rntz]))
          popval[popval<1 & !is.na(popval)] <- 1
          popval[is.na(popval)] <- 0
          
          if (all(unlist(lapply(loc1,length))==1)) {
            mat4 <- do.call(rbind,lapply(loc3,function(z)sapply(mat3[z,],function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
          } else if (any(unlist(lapply(loc1,length))==1)) {
            mat41 <- do.call(rbind,lapply(loc2,function(z)apply(mat3[z,],2,function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
            mat42 <- do.call(rbind,lapply(loc3,function(z)sapply(mat3[z,],function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
            mat4 <- rbind(mat41,mat42)
          } else {
            mat4 <- do.call(rbind,lapply(loc1,function(z)apply(mat3[z,],2,function(n)weighted.mean(n,w=popval[z],na.rm=TRUE)))) # issue with length 1 loc1 344/343 for China
          }
          
          
        } else { # no population weights
          mat4 <- do.call(rbind,lapply(loc1,function(z)colMeans(mat3[z,],na.rm=TRUE)))
        }
        
        mat5 <- t(mat4)
        for (a in 1:dim(mat5)[2]) {
          df2 <- data.frame(matrix(NA,nrow=nrow(mat5),ncol=length(colnam),dimnames=list(NULL,colnam)))
          df2$locations <- dimnames(mat5)[[2]][a]
          df2$date <- as.Date(dimnames(mat5)[[1]])
          df2[,vnm] <- mat5[,a]
          dfsav <- rbind(dfsav,df2)
        } 
      }
      
    } else { # values not accumulated
      
      mat2 <- as.matrix(unname(extval1[,3:(length(hseq)+2)])) # create matrix of hourly values
      
      for (b in tz1) { # loop time zone groups
        
        rntz <- which(cellnum$tzone==b) # row numbers with same time zones
        loc1 <- split(1:length(rntz),cellnum$locations[rntz]) # split row numbers by location
        
        if (all(unlist(lapply(loc1,length))==1)) { 
          loc3 <- loc1 # locations are all single grids
        } else if (any(unlist(lapply(loc1,length))==1)) { # split further if there are any subnational locations with a single grid
          #sploc <- which(unlist(lapply(loc1,length))==1)
          loc2 <- loc1[unlist(lapply(loc1,length))!=1] # multiple grids
          loc3 <- loc1[unlist(lapply(loc1,length))==1] # single grids
        }
        
        hdif <- tz2$dif[tz2$tz==b] # hourly difference
        daytz <- as.Date(lubridate::with_tz(hseq,tzone=b)) # days with righ time zone
        
        mat3 <- t(mat2[rntz,]) # switch hours to rows, grids to columns
        if (vnm=="t2m") { # temperature ERA5 and ERA5-Land
          mat4 <- apply(mat3,2,function(z)do.call(cbind,lapply(split(z,daytz),mean,na.rm=TRUE))) # take mean by grid
          dimnames(mat4)[[1]] <- as.character(unique(daytz))
          
        } else if (tolower(data_type)=="era5" & vnm=="tp") { # rainfall for ERA5
          mat4 <- rowsum(mat3,group=daytz,na.rm=TRUE) # sum for precipitation per grid
        }
        
        mat3 <- t(mat4) # switch days to columns, grids to rows
        
        if (pop_weights==TRUE) { # population weighted
          yrpop <- seq(2000,2020,5)
          yrdif <- as.numeric(substr(ymon1,1,4))-yrpop
          if (any(yrdif%in%-4:0)) {
            yrval <- paste0("pop",yrpop[which(yrdif%in%-4:0)])
          } else {
            yrval <- paste0("pop",yrpop[which(yrdif==min(yrdif))])  # years beyond 2020 uses 2020
          }
          popval <- unlist(split(cellnum[rntz,yrval],cellnum$locations[rntz]))
          popval[popval<1 & !is.na(popval)] <- 1
          popval[is.na(popval)] <- 0
          
          if (all(unlist(lapply(loc1,length))==1)) {
            mat4 <- do.call(rbind,lapply(loc3,function(z)sapply(mat3[z,],function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
          } else if (any(unlist(lapply(loc1,length))==1)) { # adjustment for subnational locations with single grid
            mat41 <- do.call(rbind,lapply(loc2,function(z)apply(mat3[z,],2,function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
            mat42 <- do.call(rbind,lapply(loc3,function(z)sapply(mat3[z,],function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
            mat4 <- rbind(mat41,mat42)
          } else {
            mat4 <- do.call(rbind,lapply(loc1,function(z)apply(mat3[z,],2,function(n)weighted.mean(n,w=popval[z],na.rm=TRUE)))) # potential issues/errors
          }
          
        } else { # regular mean, no population weights
          mat4 <- do.call(rbind,lapply(loc1,function(z)colMeans(mat3[z,],na.rm=TRUE))) 
        }
        
        mat5 <- t(mat4)
        for (a in 1:dim(mat5)[2]) {
          df2 <- data.frame(matrix(NA,nrow=nrow(mat5),ncol=length(colnam),dimnames=list(NULL,colnam)))
          df2$locations <- dimnames(mat5)[[2]][a]
          df2$date <- as.Date(dimnames(mat5)[[1]])
          df2[,vnm] <- mat5[,a]
          dfsav <- rbind(dfsav,df2)
        } 
      }
    } 
  }
  
  # convert to weekly if needed
  if (tolower(time_res)=="weekly") {
    if (tolower(format_weekly)=="iso") {
      dfsav$yrwk <- paste0(lubridate::isoyear(dfsav$date),"-",sprintf("%02d",lubridate::isoweek(dfsav$date)))
    } else if (tolower(format_weekly)=="epi") {
      dfsav$yrwk <- paste0(lubridate::epiyear(dfsav$date),"-",sprintf("%02d",lubridate::epiweek(dfsav$date)))
    }
    if (vnm=="tp") {
      wk.agg <- aggregate(tp~locations+yrwk,data=dfsav,FUN=sum) # weekly sums of precipitation
      wk.agg$tp <- wk.agg$tp*1000 # convert to meters
    } else if (vnm=="t2m") {
      wk.agg <- aggregate(t2m~locations+yrwk,data=dfsav,FUN=mean,na.rm=TRUE) # weekly average
      wk.agg$t2m <- wk.agg$t2m-273.15 # convert to Celsius
    }
    # reorganize
    wk.agg <- wk.agg[order(wk.agg$locations),]
    # clip to start and end year
    yrange <- as.integer(substr(time_seq[1],1,4)):as.integer(substr(time_seq[length(time_seq)],1,4))
    wk.agg <- wk.agg[substr(wk.agg$yrwk,1,4)%in%as.character(yrange),]
    row.names(wk.agg) <- 1:nrow(wk.agg)
    
    # return value
    return(wk.agg)
    
  } else if (tolower(time_res)=="daily") {
    # reorganize by location
    dfsav2 <- dfsav[order(dfsav$locations),]
    # clip to start and end year-month
    dfsav2 <- dfsav2[paste0(lubridate::year(dfsav2$date),"-",sprintf("%02d",lubridate::month(dfsav2$date)))%in%time_seq,] 
    row.names(dfsav2) <- 1:nrow(dfsav2)
    if (vnm=="tp") {
      dfsav2$tp <- dfsav2$tp*1000 # convert to meters
    } else if (vnm=="t2m") {
      dfsav2$t2m <- dfsav2$t2m-273.15 # convert to Celsius
    }
    # return value
    return(dfsav2)
  }
  
}

#rm(list=ls());gc()