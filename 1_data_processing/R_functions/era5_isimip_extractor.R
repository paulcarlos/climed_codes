######################################
# Function: Wrapper to extract variables from ERA5 (temperature and precipitation) and ISIMIP3a-TC (winds and precipitation) 
# Author: Paul LC Chua 
# Custom R functions: era5_cell_selector.R; era5_cell_extractor.R; isimip_tc_identifier.R; isimip_tc_extractor.R
# R packages: terra version 1.8-5; ncdf4 version 1.23; lubridate version 1.9.4; lutz version 0.3.2
# Version: 2025-06-11


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
# Dataframes from each extraction


######################################
# FUNCTION 
era5_isimip_extractor <- function(country_iso3, subnat_name,
                              polygon, name_location, era5_data, 
                              path_fxn, path_era5_base, path_era5_rast, path_tc_base, path_tc_rast, path_save,
                              start_year, end_year,
                              tc_model, pop_weights=FALSE, min_cover=0, no_overlap=FALSE, buffer_dist=0,
                              time_zone=FALSE, tz_pop_centroid=NULL, time_res, format_week=NULL) {
  
  # import functions
  fnctn <- c("era5_cell_selector.R", #"era5_cell_checker.R",
             "era5_cell_extractor.R", 
             "isimip_tc_identifier.R", #"isimip_tc_cell_checker.R",
             "isimip_tc_extractor.R")
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
  fn_erat2m <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_",tolower(era5_data),"-t2m_",start_year,"-",end_year,".rds")
  fn_eratp <- paste0(path_save,country_iso3,"/",country_iso3,"_",subnat_name,"_",tolower(era5_data),"-tp_",start_year,"-",end_year,".rds")
  
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
    #cat("\n\nplotting selected",paste0(toupper(era5_data)),"cell(s)/grid(s)...\n\n")
    #pdf(paste0(path_sav1,country_iso3,"_",subnat_name,"_",era5_data,"_cell-checks.pdf"),width=7,height=5)
    #era5_cell_checker(cellnum = csel_era5,
                      #polygon = polygon)
    #dev.off()
    
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
      #cat("\n\nplotting selected TC wind cell(s)/grid(s)...\n\n")
      #pdf(paste0(path_sav1,country_iso3,"_",subnat_name,"_tc-wind_cell-checks.pdf"),width=7,height=5)
      #isimip_tc_cell_checker(tc_input = csel_tcw, polygon = polygon)
      #dev.off()
      
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
      #cat("\n\nplotting selected TC rain cell(s)/grid(s)...\n\n")
      #pdf(paste0(path_sav1,country_iso3,"_",subnat_name,"_tc-rain_cell-checks.pdf"),width=7,height=5)
      #isimip_tc_cell_checker(tc_input = csel_tcr, polygon = polygon)
      #dev.off()
      
      cat("\n\nextracting TC rainfall...\n\n")
      tc_rain <- isimip_tc_extractor(tc_input = csel_tcr, file_source = path_tc_rast, pop_weights = pop_weights, 
                                     time_res = time_res, format_week = format_week)
      saveRDS(tc_rain,file=fn_tcrain)
    }
    
  } else {
    cat("file for TC rain already exists!\n")
  }
  
}

#rm(list=ls());gc()