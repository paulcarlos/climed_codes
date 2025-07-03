######################################
# Function: Extractor of tropical cyclone wind or rain from Dec 1999 to Dec 2021 from ISIMIP3a dataset
# Author: Paul LC Chua 
# R packages: terra version 1.8-5; ncdf4 version 1.23; lubridate version 1.9.4; lutz version 0.3.2
# Version: 2025-01-31


######################################
# REQUIRED INPUTS
# tc_input = List of tropical cyclone information extracted from isimip_tc_identifier function
# file_source = Path of data repository for raw netCDF's of TCs from ISIMIP3a (including resampled population rastes from GPWv4)
# pop_weights = TRUE or FALSE, TRUE uses population weights using GPWv4, FALSE takes mean/average of values 
# time_res = "Daily" or "Weekly"
# format_weekly = If weekly, "ISO" or "EPI", ISO is ISO 8601 while EPI is Epiweek of US-CDC


######################################
# OUTPUT(S) OF THE FUNCTION
# Dataframe of time series of extracted tropical cyclone variables
# NOTE: population weights derived from Gridded Population of the World, Version 4 (GPWv4): Population Count Adjusted to Match 2015 Revision of UN WPP Country Totals, Revision 11)



######################################
# FUNCTION 
isimip_tc_extractor <- function(tc_input, file_source, pop_weights, time_res, format_weekly=NULL) {
  
  # check R packages
  reqlib <- c("terra","ncdf4","lubridate","lutz")
  if (sum(reqlib%in%rownames(installed.packages()))!=length(reqlib)) {
    stop(paste0("install ",paste0(reqlib,collapse=" ")," R packages!"))
  }
  
  # check tc_input if correct
  if (sum(names(tc_input)%in%c("cells","trop_id"))!=2) {
    stop("tc_input is incorrect!")
  }
  if (sum(colnames(tc_input[["cells"]])%in%c("ID","locations","cell","coverage","pop2000","pop2005","pop2010","pop2015","pop2020","tzone"))!=10) {
    stop("tc_input is incorrect!")
  }
  if (sum(colnames(tc_input[["trop_id"]])%in%c("tcid","year","month","nature","basin"))!=5) {
    stop("tc_input is incorrect!")
  }
  
  # get attributes
  variable <- attr(tc_input,which="variable")
  model <- attr(tc_input,which="model")
  time_zone <- attr(tc_input,which="time_zone")
  
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
  if (!all(time_zone%in%OlsonNames())) {
    if (!is.logical(time_zone)) {
      stop("time_zone is incorrect!")
    } 
  } 
  cat("\nselected variable:",tolower(variable),"\n","selected model:",modname,"\n","time zone is",time_zone,"\n\n")
  
  
  # check if file source exists
  if (substr(file_source,nchar(file_source),nchar(file_source))=="/") {
    file_source1 <- file_source
  } else {
    file_source1 <- paste0(file_source,"/")
  }
  if (!dir.exists(file_source1)) {
    stop("file_source does not exist or incorrect!")
  }
  
  # loop extraction at daily level
  if (variable=="wind") {
    unitvar <- "1-minute sustained wind speeds"
  } else if (variable=="rain") {
    unitvar <- "millimeters of rainfall"
  }
  
  if (nrow(tc_input[["trop_id"]])>0) {
    message("\nextracting ",unitvar," from ",nrow(tc_input[["trop_id"]])," tropical cyclones...\n") # send notice/message
  } else {
    stop(paste0("no extractable TC ",variable,"!"))
  }
  
  colnam <- c("locations","date","tcid",variable) # all extracted by daily
  dfsav <- as.data.frame(matrix(NA,nrow=0,ncol=length(colnam),dimnames=list(NULL,colnam)))
  tzval <- unique(tc_input[["cells"]]$tzone)
  
  for (tzn in tzval) { # loop by time zone
    message("\n\n",tzn," time zone: ")
    dfcell <- tc_input[["cells"]][tc_input[["cells"]]$tzone==tzn,]
    loc1 <- split(1:nrow(dfcell),dfcell$locations) # potential issue if location has single grid
    
    if (all(unlist(lapply(loc1,length))==1)) {
      loc3 <- loc1 # locations are all single grids
    } else if (any(unlist(lapply(loc1,length))==1)) { # split further if there are subnational locations with a single grid
      #sploc <- which(unlist(lapply(loc1,length))==1)
      loc2 <- loc1[unlist(lapply(loc1,length))!=1] # multiple grids
      loc3 <- loc1[unlist(lapply(loc1,length))==1] # single grids
    }
    
    for (rw in 1:nrow(tc_input[["trop_id"]])) { # loop by tropical cyclone
      #rw=168
      #rw=342
      tid <- tc_input[["trop_id"]]$tcid[rw]
      #tcfile <- list.files(paste0(file_source1,tid),pattern=name_var,full.names=TRUE)
      #timval <- as.POSIXct(unlist(lapply(strsplit(gsub(paste0("^",file_source1,tid,"/",name_var,"_",tid,"_|","HR.tif$"),"",tcfile),"_"),
      #                                   function(z) paste0(paste(z,collapse=" "),":00"))),"UTC")
      tcfile <- terra::rast(paste0(file_source1,tid,"/",name_var,"_",tid,".nc"))
      timval <- terra::time(tcfile)
      dayval <- as.Date(lubridate::with_tz(timval,tzone=tzn))
      #mat1 <- matrix(NA,nrow=nrow(dfcell),ncol=length(tcfile)) # rows are cells, columns are hours
      mat1 <- matrix(NA,nrow=nrow(dfcell),ncol=length(timval)) # rows are cells, columns are hours
      
      #for (tc in 1:length(tcfile)) {
      for (tc in 1:length(timval)) { # loop by hour
        #r1 <- terra::rast(tcfile[tc])
        r1 <- tcfile[[tc]]
        mat1[,tc] <- unlist(r1[dfcell$cell])
        #cat(tc," ")
      }
      mat1[is.na(mat1)] <- 0
      
      if (variable=="wind") {
        mat2 <- apply(mat1,1,function(z)do.call(cbind,lapply(split(z,dayval),max,na.rm=TRUE))) # choose max winds per cell
      } else if (variable=="rain") {
        mat2 <- apply(mat1,1,function(z)do.call(cbind,lapply(split(z,dayval),sum,na.rm=TRUE))) # get sum rain per cell
      }
      if (is.matrix(mat2)) {
        dimnames(mat2)[[1]] <- as.character(unique(dayval))
      } 
      
      if (pop_weights==TRUE) { # weighted means per location
        yr <- lubridate::year(dayval[1])
        yrpop <- seq(2000,2020,5)
        yrdif <- as.numeric(yr-yrpop)
        if (any(yrdif%in%-4:0)) {
          yrval <- paste0("pop",yrpop[which(yrdif%in%-4:0)])
        } else {
          yrval <- paste0("pop",yrpop[which(yrdif==min(yrdif))])  
        }
        popval <- dfcell[,yrval]
        
        if (is.vector(mat2)) {
          mat1 <- do.call(rbind,lapply(loc1,function(z)weighted.mean(mat2[z],w=popval[z],na.rm=TRUE)))
          dimnames(mat1)[[2]] <- as.character(unique(dayval))
        } else if (any(unlist(lapply(loc1,length))==1)) {
          mat11 <- do.call(rbind,lapply(loc2,function(z)apply(mat2[,z],1,function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
          mat12 <- do.call(rbind,lapply(loc3,function(z)sapply(mat2[,z],function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
          mat1 <- rbind(mat11,mat12)
        } else {
          mat1 <- do.call(rbind,lapply(loc1,function(z)apply(mat2[,z],1,function(n)weighted.mean(n,w=popval[z],na.rm=TRUE))))
        }
        
      } else { # regular mean per location
        if (is.vector(mat2)) {
          mat1 <- do.call(rbind,lapply(loc1,function(z)mean(mat2[z],na.rm=TRUE)))
          dimnames(mat1)[[2]] <- as.character(unique(dayval))
        } else if (any(unlist(lapply(loc1,length))==1)) {
          mat11 <- do.call(rbind,lapply(loc2,function(z)apply(mat2[,z],1,function(n)mean(n,na.rm=TRUE))))
          mat12 <- do.call(rbind,lapply(loc3,function(z)sapply(mat2[,z],function(n)mean(n,na.rm=TRUE))))
          mat1 <- rbind(mat11,mat12)
        } else {
          mat1 <- do.call(rbind,lapply(loc1,function(z)apply(mat2[,z],1,function(n)mean(n,na.rm=TRUE))))
        }
        
      }
      
      df1 <- data.frame(matrix(NA,nrow=0,ncol=length(colnam),dimnames=list(NULL,colnam))) # blank dataframe to pool values
      for (a in 1:dim(mat1)[2]) { # per day processing
        if (sum(mat1[,a],na.rm=TRUE)>0) { # check if non-zero 
          df2 <- data.frame(matrix(NA,nrow=nrow(mat1),ncol=length(colnam),dimnames=list(NULL,colnam)))
          df2$locations <- dimnames(mat1)[[1]]
          df2$tcid <- tid
          df2$date <- as.Date(dimnames(mat1)[[2]][a])
          df2[,variable] <- mat1[,a]
          df1 <- rbind(df1,df2)
        }
        
      }
      df3 <- df1[df1[,variable]>0,]
      dfsav <- rbind(dfsav,df3)
      message(rw," ")
    }
  }
  
  # if weekly
  if (tolower(time_res)=="weekly") {
    message("\n\nWeekly conversion does not include tropical cyclone ID number\n")
    if (tolower(format_weekly)=="iso") {
      dfsav$yrwk <- paste0(lubridate::isoyear(dfsav$date),"-",sprintf("%02d",lubridate::isoweek(dfsav$date)))
    } else if (tolower(format_weekly)=="epi") {
      dfsav$yrwk <- paste0(lubridate::epiyear(dfsav$date),"-",sprintf("%02d",lubridate::epiweek(dfsav$date)))
    }
    if (variable=="wind") {
      wk.agg <- aggregate(wind~locations+yrwk,data=dfsav,FUN=max,na.rm=TRUE) # weekly max winds
    } else if (variable=="rain") {
      wk.agg <- aggregate(rain~locations+yrwk,data=dfsav,FUN=sum,na.rm=TRUE) # weekly rain sums
    }
    
    # reorganize
    wk.agg <- wk.agg[order(wk.agg$locations),]
    rownames(wk.agg) <- 1:nrow(wk.agg)
    
    # return dataframe
    return(wk.agg)
    
  } else if (tolower(time_res)=="daily") {
    # reorganize by location
    dfsav2 <- dfsav[order(dfsav$locations),]
    row.names(dfsav2) <- 1:nrow(dfsav2)
    
    # return value
    return(dfsav2)
  }
  
}



#rm(list=ls());gc()