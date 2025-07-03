# CODES FOR EXTRACTING CLIMATE VARIABLES FOR ClimED
This repository contains the R codes in processing and extracting ERA5-Land 2m temperatures and total precipitation, and ISIMIP tropical cyclone (TC) winds and rain by subcountry unit of each country for ClimED project.
Details of the ClimED project is [here](https://paulcarlos.quarto.pub/climed/).

## 1_data_processing
Contains codes for downloading and processing of ERA5-Land and ISIMP3a data in GRIB or NC files.

## 2_sample_extraction
Contains sample codes in extracting 2m temperature, total precipitation, TC winds and rainfall, population density, and climate zone by Costa Rican province.

## base_rasters
Contains processed files for ERA5-Land.

## functions
Contains bespoke R functions (wrapper of terra, ncdf4, lubridate, and lutz R functions) used for extraction.

## pooled
Contains TIF files aggregating TC winds and rain to indicate presence of data to be extracted.
