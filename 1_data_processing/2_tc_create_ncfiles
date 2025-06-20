library(terra)
library(ncdf4)
library(lubridate)

# import TC IDs
sid <- readRDS("tc_dataframe_isimip3a_2000-2021.rds")

# items/objects
floc <- "" # path of TC files downloaded from ISIMP data repository
fsav <- "" # path to save NC files per TC ID, converted into global rasters
fhrly <- "" # path to save hourly TIF files 

# select variable and model
#var <- "wind"; mod <- "h08"; mod1 <- "Holland 2008"; nm <- "_obsclim_historical_"; short_vnm <- "m s-1"; long_vnm <- "1-minute sustained wind speed"
var <- "rain"; mod <- "h08"; mod1 <- "Holland 2008"; nm <- "_obsclim_historical_"; short_vnm <- "mm"; long_vnm <- "total rainfall"
#var <- "wind"; mod <- "er11"; mod1 <- "Emanuel & Rotunnno 2011"; nm <- "_obsclim_historical_"; short_vnm <- "m s-1"; long_vnm <- "1-minute sustained wind speed"
#var <- "rain"; mod <- "er11"; mod1 <- "Emanuel & Rotunnno 2011"; nm <- "_obsclim_historical_"; short_vnm <- "mm"; long_vnm <- "total rainfall"

cat(paste0(long_vnm," ",mod1)) 

# blank raster entire world
wras <- rast(ext(), crs=crs("OGC:CRS84"), resolution=0.083333)
wras[] <- 0

# split TC into groups to allow multiruns
#nrow(sid)/4
grp <- split(1:nrow(sid),rep(1:4,each=ceiling(nrow(sid)/4)))

# loop 
tcdf <- sid[grp[[1]],] # CHANGE GROUPING NUMBER

for (i in 1:nrow(tcdf)) {
  id <- tcdf$tcid[i]
  cat("\n",i,": ",id,"\n")

  # open NC data and extract data
  n1 <- nc_open(paste0(floc,var,"/",mod,nm,var,"_",id,".nc")) # import 1 TC file
  lon <- as.vector(ncvar_get(n1,"lon"))
  lat <- as.vector(ncvar_get(n1,"lat"))
  if (any(lon<(-160)) & any(lon>160)) { # split NC files crossing beyond WGS84 CRS
    lon1 <- which(lon<0) # negative 
    lon2 <- which(lon>0) # positive
  }
  hr <- as.vector(ncvar_get(n1,"time")) # hourly data values
  tref <- as.POSIXct(gsub("hours since ","",n1$var$time$units),tz="UTC") # time reference
  hrv <- tref+hours(hr) # create POSIXct hours
  ar1 <- ncvar_get(n1,var)

  # file name per TC
  filnm <- paste0(fsav,id,"/",var,"_",mod,"_",id,".nc") 
  
  if (!file.exists(filnm)) {
    
    if (!is.null(dim(ar1)) & length(dim(ar1))>1) { # CHECK IF NOT EMPTY; ONLY LAND HAVE DATA
    
      if (any(lon<(-160)) & any(lon>160)) { # SPLIT MATRIX/RASTER FOR ACROSS PACIFIC OCEAN FOLLOWING WGS84
        
        for (j in 1:length(hrv)) { # hourly loop
          tim1 <- hrv[j] 
          m1 <- ar1[,,j] 
          if (length(dim(m1))==2) {
            dimnames(m1) <- list(lon,lat)
          } else {
            stop(paste0(id," incorrect lon & lat for ",tim1))
          }
          
          if (sum(m1,na.rm=TRUE)>0) {
            
            fnm1 <- paste0(fhrly,id,"/")
            if (!dir.exists(fnm1)) {
              dir.create(fnm1)
            }
            
            fnm2 <- paste0(fnm1,var,"_",mod,"_",id,"_",substr(tim1,1,10),"_",sprintf("%02d",hour(tim1)),"HR.tif")
            
            if (!file.exists(fnm2)) {
              neglon <- m1[lon1,]
              poslon <- m1[lon2,]
              if (length(dim(neglon))>1 & length(dim(poslon))>1) {
                rneg <- rast(apply(t(neglon),2,rev),crs=crs("OGC:CRS84"),extent=ext(min(lon[lon1]),max(lon[lon1]),min(lat),max(lat)))
                rpos <- rast(apply(t(poslon),2,rev),crs=crs("OGC:CRS84"),extent=ext(min(lon[lon2]),max(lon[lon2]),min(lat),max(lat)))
                resneg <- resample(rneg,wras,method="near")
                respos <- resample(rpos,wras,method="near")
                r2 <- merge(resneg,respos)
                terra::time(r2) <- tim1; names(r2) <- var; terra::units(r2) <- short_vnm
                
                writeRaster(r2, filename=fnm2,overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
                
                rm(rneg,rpos,resneg,respos);gc()
                
              } else if (is.vector(neglon) & is.vector(poslon)) {
                pt1 <- data.frame("lon"=as.numeric(lon[lon1]),"lat"=as.numeric(names(neglon)))
                pt2 <- vect(pt1,geom=c("lon","lat"),crs="OGC:CRS84")
                celn <- data.frame(terra::cells(wras,pt2,weights=T,exact=T,small=T))
                resneg <- wras
                resneg[celn$cell] <- neglon
                
                pt1 <- data.frame("lon"=as.numeric(lon[lon2]),"lat"=as.numeric(names(poslon)))
                pt2 <- vect(pt1,geom=c("lon","lat"),crs="OGC:CRS84")
                celn <- data.frame(terra::cells(wras,pt2,weights=T,exact=T,small=T))
                respos <- wras
                respos[celn$cell] <- poslon
                
                r2 <- merge(resneg,respos)
                terra::time(r2) <- tim1; names(r2) <- var; terra::units(r2) <- short_vnm
                
                writeRaster(r2, filename=fnm2,overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
                
                rm(pt1,pt2,celn,rpos,respos,rneg,resneg);gc()
                
              } else if (is.vector(neglon)) {
                pt1 <- data.frame("lon"=as.numeric(lon[lon1]),"lat"=as.numeric(names(neglon)))
                pt2 <- vect(pt1,geom=c("lon","lat"),crs="OGC:CRS84")
                celn <- data.frame(terra::cells(wras,pt2,weights=T,exact=T,small=T))
                rpos <- rast(apply(t(poslon),2,rev),crs=crs("OGC:CRS84"),extent=ext(min(lon[lon2]),max(lon[lon2]),min(lat),max(lat)))
                respos <- resample(rpos,wras,method="near")
                respos[celn$cell] <- neglon
                r2 <- respos
                terra::time(r2) <- tim1; names(r2) <- var; terra::units(r2) <- short_vnm
                
                writeRaster(r2, filename=fnm2,overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
                
                rm(pt1,pt2,celn,rpos,respos);gc()
                
              } else if (is.vector(poslon)) {
                pt1 <- data.frame("lon"=as.numeric(lon[lon2]),"lat"=as.numeric(names(poslon)))
                pt2 <- vect(pt1,geom=c("lon","lat"),crs="OGC:CRS84")
                celn <- data.frame(terra::cells(wras,pt2,weights=T,exact=T,small=T))
                rneg <- rast(apply(t(neglon),2,rev),crs=crs("OGC:CRS84"),extent=ext(min(lon[lon1]),max(lon[lon1]),min(lat),max(lat)))
                resneg <- resample(rneg,wras,method="near")
                resneg[celn$cell] <- poslon
                r2 <- resneg
                terra::time(r2) <- tim1; names(r2) <- var; terra::units(r2) <- short_vnm
                
                writeRaster(r2, filename=fnm2,overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
                
                rm(pt1,pt2,celn,rneg,resneg);gc()
              }
            
            }
          
          } 
          
          cat(j," ")
        }
        
        fnm3 <- list.files(fnm1,pattern=paste0(var,"_",mod),full.names=TRUE)
        r3 <- terra::rast(fnm3)
        writeCDF(r3, filename=paste0(fsav,id,"/",var,"_",mod,"_",id,".nc"),
                 varname=var, longname=long_vnm, unit=short_vnm, atts=c(paste0("model=",mod1)),
                 overwrite=TRUE, prec="float", compression=5)
        rm(fnm3,r3);gc()
        
      } else { # NO CROSSING OVER WGS 84 CRS
        
        for (j in 1:length(hrv)) { # hourly loop
          #j=5
          tim1 <- hrv[j] 
          m1 <- ar1[,,j] 
          
          if (sum(m1,na.rm=TRUE)>0) {
            
            fnm1 <- paste0(fhrly,id,"/")
            if (!dir.exists(fnm1)) {
              dir.create(fnm1)
            }
            fnm2 <- paste0(fnm1,var,"_",mod,"_",id,"_",substr(tim1,1,10),"_",sprintf("%02d",hour(tim1)),"HR.tif")
            
            if (!file.exists(fnm2)) {
              r1 <- rast(apply(t(m1),2,rev),crs=crs("OGC:CRS84"),extent=ext(min(lon),max(lon),min(lat),max(lat)))
              r2 <- resample(r1,wras,method="near")
              terra::time(r2) <- tim1; names(r2) <- var; terra::units(r2) <- short_vnm
              
              writeRaster(r2, filename=fnm2,overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
              rm(r1,r2);gc()
            }
            
          }
        }
        
        fnm3 <- list.files(fnm1,pattern=paste0(var,"_",mod),full.names=TRUE)
        r3 <- terra::rast(fnm3)
        writeCDF(r3, filename=paste0(fsav,id,"/",var,"_",mod,"_",id,".nc"),
                 varname=var, longname=long_vnm, unit=short_vnm, atts=c(paste0("model=",mod1)),
                 overwrite=TRUE, prec="float", compression=5)
        rm(fnm3,r3);gc()
      }
      
    } else if (length(lon)==1 & sum(ar1,na.rm=TRUE)>0) { # a vector of single grid values
      
      for (j in 1:length(hrv)) {
        #j=5
        m1 <- ar1[j]
        tim1 <- hrv[j]
        
        if (sum(m1,na.rm=TRUE)>0) {
          
          fnm1 <- paste0(fhrly,id,"/")
          if (!dir.exists(fnm1)) {
            dir.create(fnm1)
          }
          fnm2 <- paste0(fnm1,var,"_",mod,"_",id,"_",substr(tim1,1,10),"_",sprintf("%02d",hour(tim1)),"HR.tif")
          if (!file.exists(fnm2)) {
            pt1 <- data.frame("lon"=lon,"lat"=lat)
            pt2 <- vect(pt1,geom=c("lon","lat"),crs="OGC:CRS84")
            celn <- data.frame(terra::cells(wras,pt2,weights=T,exact=T,small=T))
            r2 <- wras
            r2[celn$cell] <- m1
            terra::time(r2) <- tim1; names(r2) <- var; terra::units(r2) <- short_vnm
            
            writeRaster(r2, filename=fnm2,overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
            rm(r2);gc()
          }
          
        } 
        
        cat(j," ")
      }
      
      fnm3 <- list.files(fnm1,pattern=paste0(var,"_",mod),full.names=TRUE)
      r3 <- terra::rast(fnm3)
      writeCDF(r3, filename=paste0(fsav,id,"/",var,"_",mod,"_",id,".nc"),
               varname=var, longname=long_vnm, unit=short_vnm, atts=c(paste0("model=",mod1)),
               overwrite=TRUE, prec="float", compression=5)
      rm(fnm3,r3);gc()
    }
    
    nc_close(n1)
    rm(n1,ar1);gc()
    
  }
  
}
