library(raster)
nlcd <- raster("data/nlcd2011SE")
#load from book_data path
nlcd <- raster(paste0(book_data, "nlcd2011SE"))
proj4string(nlcd)
nlcd_proj <- projection(nlcd)
res(nlcd)
str(nlcd)
plot(nlcd)
image(nlcd)
nlcd <- as.factor(nlcd) #convert to factor
levels(nlcd)


states<- readOGR("Labs/states_21basic/states.shp")
states <- spTransform(states, nlcd_proj) #set projection
proj4string(states)
summary(states)
plot(states)
plot(nlcd, add = T)



##### 3 Inverse Distance Weighting########
IDW is a simple, non-statistical interpolation technique that calculates an interpolated value at any given spatial location. We can create a raster to cover our ROI, then interpolate the ozone values at the center of each cell.

The idw() function from package gstat is easy to use.

##3.1 Template Grid
gstat::idw() needs a SpatialGridDataFrame as a template.

This code isn’t very elegant, but it works:
  
  nrows = 50; ncols = 100
ozone_template_raster = raster(nrow = nrows, ncol = ncols, ext = extent(roi_border))
ozone_template_raster[, ] = 0
ozone_grid = as(ozone_template_raster, 'SpatialGridDataFrame')

##3.2 Interpolations - The syntax is pretty easy

ozoneIDW1 = gstat::idw(formula = ozone$Av8top ~ 1, locations = ozone, newdata = ozone_grid, idp = 1)
## [inverse distance weighted interpolation]
plot(ozoneIDW1, axes = F); plot(ozone, add = T); plot(roi_counties, add = T)
