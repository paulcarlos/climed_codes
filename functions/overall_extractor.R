######################################
# Function: Wrapper to extract variables from ERA5 (temperature and precipitation), ISIMIP3a-TC (winds and precipitation), population counts, & climate zones based on polygon/point(s) 
# Author: Paul LC Chua 
# Custom R functions: era5_cell_selector.R; era5_cell_checker.R; era5_cell_extractor.R; isimip_tc_identifier.R; isimip_tc_cell_checker.R; isimip_tc_extractor.R; generic_cell_extractor.R; generic_cell_checker.R
# R packages: terra version 1.8-5; ncdf4 version 1.23; lubridate version 1.9.4; lutz version 0.3.2
# Version: 2025-01-31


######################################
# REQUIRED INPUTS
# polygon = Path or SpatVector of polygon(s) or point(s) of locations or multiple locations 
# name_location = Column name containing the name of the location or add vector of location names (No Defaults)
# path_functions = repository of the custom R functions
# era5_data = "ERA5" or "ERA5-Land" (No Defaults)
# path_era5_base = Path of data repository for ERA5 base rasters
# path_era5_rast = Path of data repository for hourly ERA5 TIFs
# path_tc_base = Path of data repository for baseline rasters of tropical cyclones
# path_tc_rast = Path of data repository for raw netCDF's of TCs from ISIMIP3a (including resampled population rastes from GPWv4)
# tc_model = "H08" or "ER11", no default (H08 is Holland 2008; ER11 is Emanuel and Rotunno 2011)
# pop_weights = TRUE or FALSE, TRUE uses population weights using GPWv4, FALSE takes mean/average of values
# min_cover = Minimum proportion of coverage of a cell/grid to consider (Default is no minimum)
# no_overlap = Assign a single location to a grid; TRUE or FALSE (Default is FALSE) selecting highest proportion of coverage
# buffer_dist = For ERA5-Land, distance in meters to widen the search for most nearest grids with values because some small islands are NA
# time_zone = TRUE or FALSE, FALSE uses UTC 0:00, brings back Olson time zones
# tz_pop_centroid = 2000, 2005, 2010, 2015, 2020, or NULL, selects time zone of the most populated grid in the chosen year, NULL uses polygon/boundary centroid
# time_res = "Daily" or "Weekly"
# format_week = If weekly, "ISO" or "EPI", ISO is ISO 8601 while EPI is Epiweek of US-CDC


######################################
# OUTPUT(S) OF THE FUNCTION
# List of dataframes from each extraction
# PDFs of cell checker plots for ERA5, TC, population, climate zone data


######################################
# FUNCTION 
overall_extractor <- function(country_iso3, subnat_name,
                              polygon, name_location, era5_data, 
                              path_fxn, path_era5_base, path_era5_rast, path_tc_base, path_tc_rast, path_save,
                              start_year, end_year,
                              tc_model, pop_weights=FALSE, min_cover=0, no_overlap=FALSE, buffer_dist=0,
                              time_zone=FALSE, tz_pop_centroid=NULL, time_res, format_week=NULL) {
  
  # import functions
  fnctn <- c("era5_cell_selector.R","era5_cell_checker.R","era5_cell_extractor.R", 
             "isimip_tc_identifier.R","isimip_tc_cell_checker.R","isimip_tc_extractor.R",
             "generic_cell_extractor.R","generic_cell_checker.R")
  for (fn in fnctn) {
    source(paste0(path_fxn,fn))
  }
  
  # check folder to save
  path_sav1 <- paste0(path_save,country_iso3,"/")
  if (!dir.exists(path_save)) {
    stop("folder to save files do not exist!")
  } else if (!dir.exists(paste0(path_sav1))) {
    dir.create(path_sav1)
  }
  
  
  ####################################
  ######### EXTRACT ERA5 #############
  
  # filenames
  fn_erat2m <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_",tolower(era5_data),"-t2m_",time_res,"_",start_year,"-",end_year,".rds")
  fn_eratp <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_",tolower(era5_data),"-tp_",time_res,"_",start_year,"-",end_year,".rds")
  
  if (any(!file.exists(c(fn_eratp,fn_erat2m)))) {
    # cell selector
    cat("\n\nselecting",paste0(toupper(era5_data)),"cell(s)/grid(s)...\n\n")
    csel_era5 <- era5_cell_selector(polygon = polygon,
                                    name_location = name_location,
                                    data_type = era5_data,
                                    file_source = path_era5_base,
                                    min_cover = min_cover,
                                    no_overlap = no_overlap,
                                    buffer_dist = buffer_dist,
                                    time_zone = time_zone,
                                    tz_pop_centroid = tz_pop_centroid)
    # print plots
    cat("\n\nplotting selected",paste0(toupper(era5_data)),"cell(s)/grid(s)...\n\n")
    pdf(paste0(path_sav1,country_iso3,"_",subnat_name,"_",era5_data,"_cell-checks.pdf"),width=7,height=5)
    era5_cell_checker(cellnum = csel_era5,
                      polygon = polygon)
    dev.off()
    # year and month
    stime <- paste0(start_year,"-01")
    etime <- paste0(end_year,"-12")
  
  } 
  # notify existence of files
  if (file.exists(fn_erat2m)) {
    cat("file for",tolower(era5_data),"temperature already exists!\n")
  } 
  if (file.exists(fn_eratp)) {
    cat("file for",tolower(era5_data),"precipitation already exists!\n")
  }
  
  if (!file.exists(fn_erat2m)) {
    # extract temperature
    cat("\n\nextracting",paste0(toupper(era5_data)),"temperatures...\n\n")
    t2m_era5 <- era5_cell_extractor(cellnum = csel_era5,
                                    file_source = path_era5_rast,
                                    start_time = stime,
                                    end_time = etime,
                                    variable = "Temperature",
                                    pop_weights = pop_weights,
                                    time_res = time_res,
                                    format_weekly = format_weekly)
    saveRDS(t2m_era5,file=fn_erat2m)
    
  }
  
  if (!file.exists(fn_eratp)) {
    # extract total precipitation
    cat("\n\nextracting",paste0(toupper(era5_data)),"total precipitation...\n\n")
    tp_era5 <- era5_cell_extractor(cellnum = csel_era5,
                                   file_source = path_era5_rast,
                                   start_time = stime,
                                   end_time = etime,
                                   variable = "Precipitation",
                                   pop_weights = pop_weights,
                                   time_res = time_res,
                                   format_weekly = format_weekly)
    saveRDS(tp_era5,file=fn_eratp)
  }
  
  
  ##################################
  ######### EXTRACT TC #############
  
  ##### TC WINDS
  fn_tcwnd <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_tc-wind_",start_year,"-",end_year,".rds")
  if (!file.exists(fn_tcwnd)) {
    cat("\n\nselecting TC IDs for winds...\n\n")
    csel_tcw <- isimip_tc_id(polygon = polygon, 
                            name_location = name_location, 
                            file_source = path_tc_base, 
                            variable = "wind", 
                            model = tc_model, 
                            year_start = start_year, 
                            year_end = end_year, 
                            min_cover = min_cover, 
                            no_overlap = no_overlap, 
                            time_zone = time_zone, 
                            tz_pop_centroid = tz_pop_centroid)
    
    if (nrow(csel_tcw[["trop_id"]])>0) {
      cat("\n\nplotting selected TC wind cell(s)/grid(s)...\n\n")
      pdf(paste0(path_sav1,country_iso3,"_",subnat_name,"_tc-wind_cell-checks.pdf"),width=7,height=5)
      isimip_tc_cell_checker(tc_input = csel_tcw, polygon = polygon)
      dev.off()
      
      cat("\n\nextracting TC winds...\n\n")
      tc_wnd <- isimip_tc_extractor(tc_input = csel_tcw, file_source = path_tc_rast, pop_weights = pop_weights, 
                                    time_res = time_res, format_week = format_week)
      saveRDS(tc_wnd,file=fn_tcwnd)
    }
    
  } else {
    cat("file for TC winds already exists!\n")
  }
  
  
  ##### TC RAINS
  fn_tcrain <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_tc-rain_",start_year,"-",end_year,".rds")
  if (!file.exists(fn_tcrain)) {
    cat("\n\nselecting TC IDs for rainfall...\n\n")
    csel_tcr <- isimip_tc_id(polygon = polygon, 
                            name_location = name_location, 
                            file_source = path_tc_base, 
                            variable = "rain", 
                            model = tc_model, 
                            year_start = start_year, 
                            year_end = end_year, 
                            min_cover = min_cover, 
                            no_overlap = no_overlap, 
                            time_zone = time_zone, 
                            tz_pop_centroid = tz_pop_centroid)
    
    if (nrow(csel_tcr[["trop_id"]])>0) {
      cat("\n\nplotting selected TC rain cell(s)/grid(s)...\n\n")
      pdf(paste0(path_sav1,country_iso3,"_",subnat_name,"_tc-rain_cell-checks.pdf"),width=7,height=5)
      isimip_tc_cell_checker(tc_input = csel_tcr, polygon = polygon)
      dev.off()
      
      cat("\n\nextracting TC rainfall...\n\n")
      tc_rain <- isimip_tc_extractor(tc_input = csel_tcr, file_source = path_tc_rast, pop_weights = pop_weights, 
                                     time_res = time_res, format_week = format_week)
      saveRDS(tc_rain,file=fn_tcrain)
    }
    
  } else {
    cat("file for TC rain already exists!\n")
  }
  
  
  ##########################################
  ######### EXTRACT POPULATION #############
  
  # variable name
  nmp1 <- "population density"; nmp2 <- "pop_density" ; nmfun <- "mean"
  #nmp1 <- "population counts"; nmp2 <- "pop_count" ; nmfun <- "sum
  
  # file names
  fn_pop <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_",nmp2,".rds")
  
  if (!file.exists(fn_pop)) {
    
    # files to extract
    cat("\n\nextracting",nmp1,"...\n\n")
    path_pop <- paste0(path_era5_rast,"tif_gpw4_population/05km/")
    files_pop <- list.files(path_pop,pattern=nmp2)
    
    # loop
    for (pf in seq(files_pop)) {
      ras1 <- terra::rast(paste0(path_pop,files_pop[pf]))
      yr <- gsub(paste("gpw4_unwpp-adjusted_",nmp2,"_|.tif"),"",files_pop[pf])
      if (pf==1) {
        df1 <- generic_cell_extractor(polygon = polygon, 
                                      name_location = name_location, 
                                      raster = ras1,
                                      min_cover = min_cover, 
                                      no_overlap = no_overlap)
        df_pop <- aggregate(values~locations,data=df1,FUN=nmfun,na.rm=TRUE)
        colnames(df_pop)[2] <- paste0("pop",yr)
      } else {
        df1 <- generic_cell_extractor(polygon = polygon, 
                                      name_location = name_location, 
                                      raster = ras1,
                                      min_cover = min_cover, 
                                      no_overlap = no_overlap)
        df2 <- aggregate(values~locations,data=df1,FUN=nmfun,na.rm=TRUE)
        colnames(df2)[2] <- paste0("pop",yr)
        df_pop <- merge.data.frame(df_pop,df2,by="locations")
      }
    }
    saveRDS(df_pop,fn_pop)
    
  } else {
    cat("file for",nmp1,"already exists!\n")
  }
  
  #############################################
  ######### EXTRACT CLIMATE ZONES #############
  # filename
  fn_cz <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_climzone.rds")
  
  if (!file.exists(fn_cz)) {
    # climate zone raster
    cat("\n\nextracting climate zones...\n\n")
    czras <- terra::rast(paste0(path_era5_rast,"tif_beck_climate-zone/1991_2020/koppen_geiger_0p1.tif"))
    refcz <- read.csv(paste0(path_era5_rast,"tif_beck_climate-zone/beck_climzone_legend.csv"))
    # population raster
    pras1 <- terra::rast(paste0(path_era5_rast,"tif_gpw4_population/05km/gpw4_unwpp-adjusted_pop_count_2010.tif"))
    pras_ag <- terra::aggregate(pras1,fact=2) # get it closer to resolution of raster
    pras2 <- terra::resample(pras_ag,czras,method="bilinear")
    # dataframe
    df1 <- generic_cell_extractor(polygon = polygon, 
                                  name_location = name_location, 
                                  raster = czras,
                                  min_cover = min_cover, 
                                  no_overlap = no_overlap)
    df1$pop <- unlist(pras2[df1$cell])
    # get the most populated climate zones per area
    loc <- unique(df1$locations)
    df_cz <- data.frame("locations"=loc,"climzone"=NA)
    for (cz in seq(loc)) {
      sdf <- df1[df1$locations==loc[cz],]
      if (sum(sdf$pop,na.rm=TRUE)==0) {
        if (length(unique(sdf$values))==1) {
          df_cz$climzone[df_cz$locations==loc[cz]] <- refcz$climcode[refcz$value==unique(sdf$values)]
        } else {
          df_cz$climzone[df_cz$locations==loc[cz]] <- refcz$climcode[refcz$value==sdf$values[1]]
        }
      } else if (sum(sdf$pop,na.rm=TRUE)>0) {
        ag1 <- aggregate(pop~values,data=sdf,FUN="sum",na.rm=TRUE)
        df_cz$climzone[df_cz$locations==loc[cz]] <- refcz$climcode[refcz$value==ag1$values[ag1$pop==max(ag1$pop)]]
      }
      
    }
    saveRDS(df_cz,fn_cz)
  } else {
    cat("file for climate zones already exists!\n")
  }
  
  
  
  #########################################
  ######### COMBINE IN A LIST #############
  #dlist <- list()
  #dlist[["t2m"]] <- t2m_era5
  #dlist[["tp"]] <- tp_era5
  #dlist[["tcwind"]] <- tc_wnd
  #dlist[["tcrain"]] <- tc_rain
  #dlist[["pop"]] <- df_pop
  #dlist[["climzone"]] <- df_cz
  #attr(dlist,which="country_iso3") <- country_iso3
  #attr(dlist,which="subnat_name") <- subnat_name
  #attr(dlist,which="era5_data") <- era5_data
  #attr(dlist,which="tc_model") <- tc_model
  #attr(dlist,which="min_cover") <- min_cover
  #attr(dlist,which="no_overlap") <- no_overlap
  #attr(dlist,which="start_year") <- start_year
  #attr(dlist,which="end_year") <- end_year
  #attr(dlist,which="buffer_dist") <- buffer_dist
  #attr(dlist,which="pop_weights") <- pop_weights
  #attr(dlist,which="time_zone") <- time_zone
  #attr(dlist,which="tz_pop_centroid") <- tz_pop_centroid
  #attr(dlist,which="time_res") <- time_res
  #attr(dlist,which="format_week") <- format_week
  
  #saveRDS(dlist,paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_list_era5-tc-pop-climzone_",start_year,"-",end_year,".rds"))
}

#rm(list=ls());gc()