######################################
# Function: Plotting selected cells and the accompanying polygon(s)/point(s) of administrative border for tropical cyclone ISIMIP3a dataset
# Author: Paul LC Chua 
# R package(s): terra version 1.8-5 
# Version: 2025-01-31


######################################
# REQUIRED INPUTS
# tc_input = List of tropical cyclone information extracted from isimip_tc_identifier function
# polygon = Path or SpatVector of polygon of locations or multiple locations


######################################
# OUTPUT(S) OF THE FUNCTION
# Plot of grid(s) containing proportion of coverage within polygon of administrative border overlapped to the actual borders
# NOTE: The grids with color(s) are grids with value 


######################################
# FUNCTION
isimip_tc_cell_checker <- function(tc_input, polygon) {
  
  # check R packages
  reqlib <- c("terra")
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
  file_source <- attr(tc_input,which="file_source")
  
  # check file source
  if (!dir.exists(file_source)) {
    stop("file_source does not exist or incorrect!")
  }
  
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
  message("\nselected variable:",tolower(variable),"\n","selected model:",modname,"\n\n")
  
  # check polygon
  if (class(polygon)[1]=="SpatVector") {
    poly1 <- terra::project(polygon,y=terra::crs("OGC:CRS84"))
  } else if (file.exists(polygon)) {
    poly1 <- terra::vect(polygon)
    poly1 <- terra::project(poly1,y=terra::crs("OGC:CRS84"))
  } else {
    stop("polygon not a SpatVector or a correct path/file name")
  }
  name_location <- attr(tc_input,which="name_location")
  if (sum(names(poly1)==name_location)==0) {
    stop("name_location is incorrect or does not match with polygon!")
  }
  
  # upload base raster
  ras1 <- terra::rast(paste0(file_source,name_var,"_overall.tif"))
  
  # create alternative empty raster 
  ras2 <- ras1
  ras2[] <- NA
  ras2[!is.na(ras1) & ras1>0] <- -0.1
  
  # get cells dataframe
  cellnum <- tc_input[["cells"]]
  
  # get locations
  loc <- unique(cellnum$locations) # locations
  cat("\nplotting proportion of cell coverages for ",length(loc)," locations... \n \n") # send notice/message
  
  # loop plots
  if (length(loc)>1) {
    par(mfrow=c(2,2),mar=c(3,5,4,2))
  } else {
    par(mar=c(3,5,4,2))
  }
  for (x in loc) {
    poly2 <- poly1[which(poly1[[name_location]]==x)] # get polygon
    extnt <- terra::ext(poly2) # extent of plot
    extnt[c(1,3)] <- floor(extnt[c(1,3)])
    extnt[c(2,4)] <- ceiling(extnt[c(2,4)])
    #ras2[cellnum$cell] <- -0.1
    ras3 <- ras2
    ras3[cellnum$cell[cellnum$locations==x]] <- cellnum$coverage[cellnum$locations==x]
    ras4 <- terra::crop(ras3,extnt)
    terra::plot(ras4,col=c("#B2BEB5","#219EBC","#023047","#FFB703","#FB8500"),
                breaks=c(-0.1,0,0.25,0.5,0.75,1.1),
                main=paste0(x," \n(proportion of cell coverage)"),cex.main=0.9)
    if (terra::geomtype(poly2)%in%c("polygon","polygons")) {
      terra::lines(poly2,col="red",lwd=3)
    } else if (terra::geomtype(poly2)%in%c("point","points")) {
      terra::points(poly2,col="red")
    }
  }
  message("\ngrey colored cells have ",variable," values and white colored cells are empty cells\n")
  
}

#rm(list=ls());gc()