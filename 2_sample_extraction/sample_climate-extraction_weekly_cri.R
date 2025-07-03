# changing info / objects by country 
country_iso3 <- "cri" # COSTA RICA
subnat_name <- "prov"
polygon <- terra::vect("cri_prov_shapefile.gpkg")
#names(polygon)
name_location <- "prov"
start_year <- 2014
end_year <- 2021
time_zone <- TRUE
tz_pop_centroid <- 2020 # 2000,2005,2010,2015,2020
format_weekly <- "ISO" #ISO EPI

# fixed details / objects
era5_data <- "ERA5-Land" 
fsourc <- "D"
path_fxn <- "functions/" # FOLDER CONTAINING THE BESPOKE R FUNCTIONS
path_era5_base <- paste0(fsourc,":/rast_era5-land/base_rasters/")
path_era5_rast <- paste0(fsourc,":/") 
path_tc_base <- paste0(fsourc,":/rast_isimip3a_tc/base_rasters/") 
path_tc_rast <- paste0(fsourc,":/rast_isimip3a_tc/tcid/") 
path_save <- "E:/R/dia_outc/data/extracted_vars/" 
tc_model <- "ER11"
pop_weights <- TRUE
min_cover <- 0.1 
no_overlap <- FALSE 
buffer_dist <- 100000 #meters
time_res <- "weekly"


# extract
source("functions/overall_extractor.R")
overall_extractor(country_iso3, 
                  subnat_name,
                  polygon, 
                  name_location, 
                  era5_data, 
                  path_fxn, 
                  path_era5_base, 
                  path_era5_rast, 
                  path_tc_base, 
                  path_tc_rast, 
                  path_save,
                  start_year, 
                  end_year,
                  tc_model, 
                  pop_weights, 
                  min_cover, 
                  no_overlap, 
                  buffer_dist,
                  time_zone, 
                  tz_pop_centroid, 
                  time_res, 
                  format_weekly)
