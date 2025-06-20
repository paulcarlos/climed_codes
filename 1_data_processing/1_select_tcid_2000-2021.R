library(terra)
library(ncdf4)
library(lubridate)

# import
nc <- nc_open("tracks_obsclim_historical_1950_2021.nc")

# world map; for visualization purposes
wmap <- vect("ne_10m_land.shp")

# get IDs and lon lat
sid <- ncvar_get(nc,varid="sid")
tim <- ncvar_get(nc,varid="time")
lat <- ncvar_get(nc,varid="lat") 
lon <- ncvar_get(nc,varid="lon")
tunit <- as.POSIXct(sub("hours since ","",nc$var$time$units),"UTC")
basin <- ncvar_get(nc,varid="basin")
nature <- ncvar_get(nc,varid="nature")
wind <- ncvar_get(nc,varid="windspatialmax")

# select TCs from December 1999 to December 2021
yseq <- 2000:2021
sid1 <- which(as.integer(substr(sid,1,4))>1998)
df <- data.frame("tcid"=character(0),"year"=character(0),"month"=character(0),
                 "nature"=character(0),"basin"=character(0))
for (i in sid1) {
  #i=sid1[200]
  #i=which(sid=="2000244N13229")
  t1 <- substr(tunit+hours(tim[,i]),1,7)
  w1 <- wind[,i]
  w2 <- which(w1>0)
  #pt <- buffer(vect(gpt,geom=c("lon","lat"),crs=crs("OGC:CRS84")),width=10000)
  #plot(wmap); plot(pt,col=2,add=TRUE)
  #plot(pt,col=2); plot(wmap,add=TRUE)
  if (any(as.integer(substr(t1,1,4))%in%yseq) & any(w1>0)) {
    sdf <- data.frame("tcid"=sid[i],"year"=substr(t1[w2[length(w2)]],1,4),"month"=substr(t1[w2[length(w2)]],6,7),
                      "nature"=paste0(unique(nature[,i]),collapse=";"),"basin"=paste0(unique(basin[,i]),collapse=";"))
    df <- rbind(df,sdf)
  }
}

# save 
saveRDS(df,"tc_dataframe_isimip3a_2000-2021.rds")
