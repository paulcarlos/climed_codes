library(terra)

# guiding/reference dataframe
ref <- readRDS("0_data/tc_dataframe_isimip3a_2000-2021.rds")

# objects
floc <- "" # path of individual processed TCs into global rasters
fsav <- "" # path to save pooled rasters
basin <- unique(substr(ref$basin,1,2))

# select variable and model
#var <- "wind"; mod <- "er11"
var <- "rain"; mod <- "er11"
#var <- "wind"; mod <- "h08"
#var <- "rain"; mod <- "h08"

# loop aggregate
#i=basin[1]
thr <- 300
ilist <- list()
for (i in basin) { # generate single raster per TC per basin
  cat("\n",i,"\n")
  file1 <- list.files(floc,pattern=paste0(var,"_",mod,"_",i))
  slist <- list()
  if (length(file1)>thr) {
    grp <- split(1:length(file1),rep(1:ceiling(length(file1)/thr),each=thr))
    for (j in seq(grp)) {
      #j=1
      sgrp <- grp[[j]]
      f1 <- file1[sgrp[1]]
      r1 <- rast(paste0(floc,f1))
      #r1[is.na(r1)] <- 0
      #plot(r1)
      v1 <- gsub(paste0(var,"_",mod,"_",i,"_|.tif"),"",f1)
      for (k in 2:length(sgrp)) {
        #k=2
        f2 <- file1[sgrp[k]]
        r2 <- rast(paste0(floc,f2))
        #r2[is.na(r2)] <- 0
        r1 <- sum(r1,r2,na.rm=TRUE) # sum is fastest, choosing max will be too slow
        #cat(k," ")
        v2 <- gsub(paste0(var,"_",mod,"_",i,"_|.tif"),"",f2)
        v1 <- c(v1,v2)
      }
      slist[[j]] <- ref[ref$tcid%in%v1,]
      #plot(r1)
      writeRaster(r1,filename=paste0(fsav,var,"_",mod,"_",i,"_grp",j,".tif"),overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
      cat(j," ")
    }
  } else {
    f1 <- file1[1]
    r1 <- rast(paste0(floc,f1))
    v1 <- gsub(paste0(var,"_",mod,"_",i,"_|.tif"),"",f1)
    for (k in 2:length(file1)) {
      f2 <- file1[k]
      r2 <- rast(paste0(floc,f2))
      r1 <- sum(r1,r2,na.rm=TRUE)
      v2 <- gsub(paste0(var,"_",mod,"_",i,"_|.tif"),"",f2)
      v1 <- c(v1,v2)
    }
    writeRaster(r1,filename=paste0(fsav,var,"_",mod,"_",i,"_grp1.tif"),overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
    slist[[1]] <- ref[ref$tcid%in%v1,]
    cat(1," ")
  }
  ilist[[i]] <- slist
  rm(r1,r2);gc()
}
saveRDS(ilist,paste0(fsav,var,"_",mod,"_list-tcid.rds"))

# pile up all rasters
lf <- list.files(fsav,pattern=paste0("^",var,"_",mod,".*\\.tif$"))
r1 <- rast(paste0(fsav,lf[1]))
for (i in 2:length(lf)) {
  r2 <- rast(paste0(fsav,lf[i]))
  #plot(r2)
  r1 <- sum(r1,r2,na.rm=TRUE)
  #plot(r1)
  cat(i," ")
}
writeRaster(r1,filename=paste0(fsav,var,"_",mod,"_overall.tif"),overwrite=TRUE,gdal=c("COMPRESS=ZSTD","TFW=NO"))
plot(r1)

  
rm(list=ls());gc()

