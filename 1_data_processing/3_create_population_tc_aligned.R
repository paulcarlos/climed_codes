library(terra)

# load population data
pop <- rast("gpw_v4_population_count_adjusted_rev11_2pt5_min.nc") # ~5km resolution data from NASA/CIESIN
pop <- pop[[1:5]] # select every 5 years
names(pop)

# years
yr <- c("2000"=1,"2005"=2,"2010"=3,"2015"=4,"2020"=5)

# loop
fsav <- "" # location to save the base population raster
for (i in seq(yr)) {
  yval <- names(yr)[i]
  p1 <- pop[[yr[i]]]
  aggr <- aggregate(p1,fact=2,fun="sum")
  writeRaster(aggr,filename=paste0(fsav,"gpw411_pop-count_",yval,"_isimip3a-tc.tif"),overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
}

rm(list=ls());gc()
