library(terra)

# load population data
pop <- rast("gpw_v4_population_count_adjusted_rev11_2pt5_min.nc") # downloaded from NASA/CIESIN
pop
names(pop)

# years
yr <- c("2000"=1,"2005"=2,"2010"=3,"2015"=4,"2020"=5)
p1 <- pop[[1]]
plot(p1)

# load ERA5-Land sample
era <- rast("tp_hourly_era5land_global_1999-12-01-00.grib")
#era <- rast("tp_hourly_era5_global_1999-12-01-00.grib")
era

# aggregate and align grids
p1 <- pop[[5]]
p1
sum(p1[],na.rm=T)
p2 <- project(p1,era,method="sum")
p2 <- terra::aggregate(p1,fact=2,fun="sum",na.rm=T)
p2
sum(p2[],na.rm=T)
ncell(p2)
p3 <- project(p2,era,method="sum")
p3
sum(p3[],na.rm=T)
plot(p3)


# loop extraction
#i=yr[1]
fsav <- "" # path to save files
for (i in yr) {
  p1 <- pop[[i]]
  p2 <- resample(p1,era,method="sum")
  writeRaster(p2,filename=paste0(fsav,"/unwpp_pop_count_",names(yr[yr==i]),".tif"),overwrite=TRUE,
              gdal=c("COMPRESS=ZSTD","TFW=NO"))
  cat(i," ")
}

rm(list=ls());gc()
