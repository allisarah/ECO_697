########LAB 1
####PART 2######
book_data = "C:/Users/aalli/Documents/ECO 697/data/"
dir.exists(book_data)
file.exists("data/environment_vars.R")

file.exists(paste0(book_data, "nlcd2011SE"))
library(raster)
nlcd = raster(paste0(book_data, "nlcd2011SE"))
print(nlcd)
require(raster)


######## chap.2.3.3 pg 28 of text ########
set.seed(16)#sets random number seed for repeatability
toy <- raster(ncol=6, nrow=6, xmn=1, xmx=6, ymn=1, ymx=6)
values(toy) <- rpois(ncell(toy), lambda=3)
ncell(toy)
plot(toy)
text(toy, digits=2)

ncell(toy)
toy2 <- toy
values(toy2) <- 1:ncell(toy)
plot(toy2)
text(toy2, digits=2)

toy_mean <- aggregate(toy, fact=2, fun=mean) #mean value
toy_maj <- aggregate(toy, fact=3, fun=modal) #majority rule

cellStats(toy, mean)
cellStats(toy, var)

cellStats(toy_mean, mean)
cellStats(toy_mean, var)

toy_dis2 <- disaggregate(toy, fact=2)
toy_dis2_bilinear <- disaggregate(toy, fact=2,method='bilinear')

#decrease the extent
e <- extent(2, 4, 2, 4) #first create new, smaller extent
toy_crop <- crop(toy, e)
plot(toy_crop)
text(toy_crop)

#increase the extent
e <- extent(0, 7, 0, 7) #first create new, larger extent
toy_big <- extend(toy, e)
plot(toy_big)
text(toy_big)

#raster similar to toy
set.seed(16)#sets random number seed for repeatability
emp_temp<- raster(ncol=6, nrow=6, xmn=1, xmx=6, ymn=1, ymx=6)
values(emp_temp) <- rpois(ncell(emp_temp), lambda=3)
ncell(emp_temp)
plot(emp_temp)
text(emp_temp, digits=2)

#raster with lambda = 10
set.seed(16)#sets random number seed for repeatability
rast<- raster(ncol=6, nrow=6, xmn=1, xmx=6, ymn=1, ymx=6)
values(rast) <- rpois(ncell(rast), lambda=10)
ncell(rast)
plot(rast)
text(rast, digits=2)

##template with 7 rows
rast7<- raster(ncol=7, nrow=7, xmn=1, xmx=6, ymn=1, ymx=6)

e <- extent(0, 7, 0, 7) #first create new, larger extent
toy_big <- extend(toy, e)
plot(toy_big)
text(toy_big)

####PART 3 LAB1
library(raster)
library(sp)
library(rgdal)
library(rgeos)

##function for creating raster template
rast_templ = function(nrows, ncols, xmn = 1, xmx = 6, ymn = 1, ymx = 6)
{
  return(raster(
    ncol = ncols, nrow = nrows, 
    xmx = xmx, xmn = xmn, 
    ymx = ymx, ymn = ymn))
}

#par(mfrow = c(2,2), mar = c(.5,.5,.5,.5))
par(mfrow = c(2,2), mar = c(1,1,1,1))
dev.off()

par(mfrow = c(2,2), mar = c(1,1,1,1))
####resample
rast1 <- rast_templ(12, 12)
rast4 <- rast_templ(7, 7)
plot(resample(rast, rast4))

plot(resample(rast, rast_templ(10, 10)))
plot(resample(rast, rast_templ(4, 4)))

rast2<-resample(rast, rast1)
plot(rast2)
plot(rast)


##recover original rast
rast3<-resample(rast2, rast)
plot(rast3)
plot(rast)


#####LAB 2
nlcd1<-raster(paste0(book_data, "nlcd2011SE"))
proj4string(nlcd1)

#nlcd_proj <- projection(nlcd)
res(nlcd1)
str(nlcd1)
plot(nlcd1)
image(nlcd1)
nlcd1 <- as.factor(nlcd1) #convert to factor
levels(nlcd1)


states<- readOGR("Labs/states_21basic/states.shp")
states <- spTransform(states, nlcd_proj) #set projection
proj4string(states)
#states_nlcd<-extent(nlcd)
nlcd_sta<- extend(nlcd, states)
proj4string(nlcd)
summary(states)
plot(nlcd_sta)
plot(nlcd, add = T)


nlcd_ext <- extent(states)
state_nlcd <- crop(nlcd, nlcd_ext)
plot(state_nlcd)

e <- extent(2, 4, 2, 4) #first create new, smaller extent
toy_crop <- crop(toy, e)
plot(toy_crop)
text(toy_crop)


####LAB2 

########QUESTION3#############
library(raster)
#nlcd_ext <- extent(nlcd)
#e <- as(extent(nlcd), "SpatialPolygons")
#proj4string(e) <- proj4string(states)
#plot(states)
#plot(e, add = T)
#e <- as(extent(c(-949552.2, 1007165.9, 1164969, 1856979)), "SpatialPolygons")

##current answer
locator(2)
e <- as(extent(c(-639340.8, 1388964.6, 902482.5, 2047878.4)), "SpatialPolygons")
plot(states)
plot(e, add = T, col = rgb(0, 1, 1, 0.3))


###Q4
##create 0 and 1 values based on the reclassification of forest type (cultivated&non-cultivated)
#reclass1 <- c(rep(0, 13), rep(1, 1), rep(0, 2)) ##cultivated vs non-cultivated
reclass2 <- c(rep(0, 2), rep(1, 4), rep(0, 10))
nlcd.levels1 <- levels(nlcd)[[1]]

#create reclassify matrix: first col: orginal; second: change to reclass values created
reclass.mat1 <- cbind(levels(nlcd)[[1]], reclass2)
reclass.mat1
head(reclass.mat1, 3)

#reclassify
forest1 <- reclassify(nlcd, reclass.mat1)

grainarea1 <- res(forest1)[[1]]^2/10000  ## from meter to ha, 10,000 m^2 = 1 ha^2 

#create empty vector for storing output
f1km <- rep(NA, length = nrow(sites))
f.5km <- rep(NA, length = nrow(sites))
f2km <- rep(NA, length = nrow(sites))
f3km <- rep(NA, length = nrow(sites))
f4km <- rep(NA, length = nrow(sites))
f5km <- rep(NA, length = nrow(sites))

#with for loop (all five buffers: 910s; <=3km: 228s)
for(i in 1:nrow(sites)) {
  f3km[i] <- BufferCover(coords=sites,size=3000,landcover=forest1,grain=grainarea1)
  f4km[i] <- BufferCover(coords=sites,size=4000,landcover=forest1,grain=grainarea1)
  print(i)
}

for(i in 1:nrow(sites)) {
  f1km[i] <- BufferCover(coords=sites,size=1000,landcover=forest1,grain=grainarea1)
  f2km[i] <- BufferCover(coords=sites,size=2000,landcover=forest1,grain=grainarea1)
  f.5km[i] <- BufferCover(coords=sites,size=500,landcover=forest1,grain=grainarea1)
  f5km[i] <- BufferCover(coords=sites,size=500,landcover=forest1,grain=grainarea1)
  print(i)
}

#make a data frame
forest.scale1 <- data.frame(site=sites$site,
                           x=sites$coords_x1, y=sites$coords_x2,
                           f3km=f3km, f4km=f4km)

#plot
plot(f.5km, f1km)
plot(f1km, f2km)
plot(f2km, f3km)
plot(f3km, f4km)
plot(f4km, f5km)

#aggregate forest landcover to grain reflected in sampling
forest200_1 <- aggregate(forest1, fact=7, fun=modal)

#inspect
res(forest200_1) #210x210 res

#plot
plot(forest1)
plot(states, add = T)
plot(sites, pch=20, col="red", add=T)




###LAB3####
###QUESTIONS#####
cactus <- read.csv("cactus.csv")
boundary <- read.csv("cactus_boundaries.csv",header=T)

#create spatstat objects
ppp.window1 <- owin(xrange=c(0, 10),
                   yrange=c(0, 10))

#Simulate a homogeneous point pattern 4 times
sim.pp1 <- rpoispp(lambda=intensity(ppp.cactus), nsim=4, win=ppp.window)
sim.pp2 <- rpoispp(lambda= 10 * intensity(ppp.cactus), nsim=4, win=ppp.window)

#plot
plot(sim.pp1)
plot(sim.pp2)
plot(sim.pp2, rgb(0, 0, 0, 0.2))

###with 20,000 points
sim.pp3 <- rpoispp(lambda= 250 * intensity(ppp.cactus), nsim=1, win=ppp.window)
plot(sim.pp3, col = rgb(0, 0, 1, 0.1), pch = 16)
#sum(sim.pp3$n)

#confidence intervals
csr_env_cactus <- envelope(ppp.cactus, Lest, nsim=99, rank=1, correction="trans", 
                           global=T, savefuns = TRUE)

envelope_sims = data.frame(attr(csr_env_cactus, "simfuns"))
head(envelope_sims)[, 1:5]

#L-function plots are easier to read if we subtract the value of the radius from the L-values.
envelope_sims2 = envelope_sims
envelope_sims2[, -1] = envelope_sims[, -1] - envelope_sims2$r

#plot
matplot(x = envelope_sims2$r, y = envelope_sims2[, -1][, 1:5],
        type = "l", ylab = "L(r)", xlab = "r")

#plot many lines
matplot(x = envelope_sims2$r, y = envelope_sims2[, -1][, 1:99],
        type = "l", ylab = "L(r)", xlab = "r", col = gray(level = 0.1, alpha = 0.2),
        ylim = c(-2,2.5))

lines(Liso$r, Liso$iso- Liso$r, type = "l", col = "red")

#Q4
# Repeat 4.3.5 - Marked point patterns for a spatstat dataset
##############################################
View(waka)
View(anemones)

View(longleaf)
plot(longleaf)
summary(longleaf$marks) #median is 26.15
View(boundary)
#ppp.window1 <- owin(xrange=c(min(longleaf$x), max(longleaf$x)),
                   #yrange=c(min(longleaf$y), max(longleaf$y)))

ppp.longleaf <- ppp(longleaf$x, longleaf$y, window=longleaf$window)

#For  call 'Z', which has the data on marks/diameter
longleaf$marksBS<- as.factor(ifelse(longleaf$marks >26, "big", "small"))

#makes a new object with marks, amenable to analysis in spatstat
ppp.BS <- ppp(longleaf$x, longleaf$y, window=longleaf$window, marks=longleaf$marksBS)


#inspect
summary(ppp.BS)

#plot
plot(ppp.BS)
plot(split(ppp.BS))

#Point pattern only considering marks/diameter
marks.data <- subset(longleaf, marks>26)
ppp.marks <- ppp(marks.data$x, marks.data$y, window=longleaf$window)

Lmarks <- envelope(ppp.marks, Lest, nsim=99, rank=1, i="big", global=F)

#plot with envelope
plot(Lmarks, . - r~r, shade=c("hi", "lo"), legend=F)

#Bivariate BS based on a random-labeling:
Lmulti1 <- envelope(ppp.BS, Lcross, nsim=99, rank=1, i="big", global=F,
                   simulate=expression(rlabel(ppp.BS)))

#plot with envelope
plot(Lmulti1, . - r~r, shade=c("hi", "lo"), legend=F)

#-------------------------------#
#Continuous marks
#-------------------------------#
#####dont know how to calculate area for longleaf
#makes a new object with marks, amenable to analysis in spatstat
ppp.area1 <- ppp(longleaf$x, longleaf$y, window=longleaf$window, marks = longleaf$marks)

#ppp.area1 <- ppp(longleaf$x, longleaf$y, window=longleaf$window, marks=area(longleaf$marks))
#area1<- area(longleaf$marks)

#inspect
summary(ppp.area1)

#plot
plot(ppp.area1)

#mark correlation function
mcf.area1 <- markcorr(ppp.area1)

#plot
plot(mcf.area1, legend=F)

#Monte Carlo simulations for confidence envelope
MCFenv1 <- envelope(ppp.area1, markcorr, nsim=99, correction="iso", global=F)

#plot with envelope
plot(MCFenv1,  shade=c("hi", "lo"), legend=F)



#Q6
#Here is an adaptation of part of the code I used to make the Matern cluster figure in the slides:
#kappa is the intensity (lambda) of the parent points
kappa_1 = 1
# scale is the radius of the offspring clusters
scale_1 = 0.2
# mean number of points per cluster.
# (I think the final count is Poisson-distributed -- have to dig into the code to see.)
mu_1 = 4
# Different sized windows for scale comparison
win_small = owin(xrange = c(0, 10), yrange = c(0, 10))
win_large = owin(xrange = c(0, 100), yrange = c(0, 100))
# plotting parameters
pch_mat = 16
cex_mat = 0.8
mat_pp_small = rMatClust(kappa_1, scale_1, mu_1, win = win_small)
mat_pp_large = rMatClust(kappa_1, scale_1, mu_1, win = win_large)
par(mfrow = c(2, 1))

#xlab and ylab not working??
plot(mat_pp_small$x, mat_pp_small$y, pch = pch_mat, cex = cex_mat,
     xlab = "x-coordinates", ylab = "y-coordinates", 
     main = "Matern Process - 10 by 10 window ")

plot(mat_pp_large$x, mat_pp_large$y, pch = pch_mat, cex = 0.3 * cex_mat, 
     xlab = "x-coordinates", ylab = "y-coordinates",
     main = "Matern Process - 100 by 100 window ")

# aspect ratio = asp = 1 makes the plot square

##plot of matern with l function
ppp.mat <- ppp(mat_pp_large$x, mat_pp_large$y, window=mat_pp_large$window)
ppp.mat1 <- ppp(mat_pp_small$x, mat_pp_small$y, window=mat_pp_small$window)

Lmat <- envelope(ppp.mat, Lest, nsim=99, rank=1, i="big", global=F)
Lmat1 <- envelope(ppp.mat1, Lest, nsim=99, rank=1, i="small", global=F )
Lmat2 <- envelope(ppp.mat, Lest, nsim=99, rank=1, i="small", global=F,
                  funargs = list("rmax" = 20), correction = "iso")
#funargs = list("rmax" = 20), correction = "iso"

#plot with envelope
plot(Lmat, . - r~r, shade=c("hi", "lo"), legend=F, main = "Lmat 10 by 10 window")
plot(Lmat1, . - r~r, shade=c("hi", "lo"), legend=F, main = "Lmat 100 by 100 window")
plot(Lmat2, . - r~r, shade=c("hi", "lo"), legend=F, main = "Lmat 100 by 100 window")
################################


#######LAB 5 ##############
library(rgdal)
#ia_counties = readOGR(here("data/census_2010_county_DP1"))
ia_counties = readOGR("tl_2010_19_county10_Iowa")
#plot water area for Iowa counties
blue_pal = colorRampPalette(c(rgb(0.2, 0.2, 0.5), rgb(0.7, 0.7, 1.0)))
county_colors_1 = blue_pal(nrow(ia_counties))[order(ia_counties$AWATER10)]
plot(ia_counties, col = county_colors_1, main = "Iowa counties - water area")

#plot
#Question - best way to add legend to cloropleth maps
county_colors_2 = heat.colors(nrow(ia_counties))[order(ia_counties$ALAND10)]
plot(ia_counties, col = county_colors_2, main = "Iowa counties - land area")
legend("bottomright", legend = ia_counties$ALAND10, fill = county_colors_2, cex = 0.2)

#marked point example
#read csv file
cactus = read.csv(paste0(book_data, "cactus.csv"))
class(cactus)
#specify coordinates to transfrom dataframe to spatialpointsdataframe
coordinates(cactus) = ~ East + North
class(cactus)
View(cactus@coords)

####worldclim data
dat = raster::getData('worldclim', var='tmin', res=10, lon=c(4, 5), lat=c(45, 46))
View(dat)
dat
class(dat)
plot(subset(dat, 1:3))
plot(subset(dat, 1))

####NetCDF example
library(ncdf4)
dat = raster("sresa1b_ncar_ccsm3-example.nc")
str(dat, 0)
plot(dat)

#####OZONE
cal_county <- readOGR("data/tl_2010_06_Cal_county00/tl_2010_06_county00.shp")
extent(ozone)
extent()
locator(1)
la_county <- crop(cal_county, extent(-118.7, -116, 33.6, 34.8))

plot(la_county,  main = "Ozone concentrations in California")
states<- readOGR("Labs/states_21basic/states.shp")
ozone = read.csv(paste0(book_data, "ozone.csv"))
class(ozone)
summary(ozone)
#specify coordinates to transfrom dataframe to spatialpointsdataframe
coordinates(ozone) = ~ Lon + Lat
class(ozone)
summary(ozone)
View(ozone@coords)

proj4string(ozone) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
summary(ozone)
plot(ozone)

#color trick
ozone_colors = terrain.colors(nrow(ozone), alpha = 1)[order(ozone$Station)]
ozone_colors

points(ozone, col = ozone_colors, cex = (ozone$Av8top)/2, pch = 16)
#legend("bottomright", fill = ozone_colors, legend = order(ozone$Av8top, ozone_colors))

#####USA states
extent(states)
proj4string(states) ###albers equal area
locator(2)

e1 <- as(extent(c(-2374592,  2168972, -100841.4, 3291271.2)), "SpatialPolygons")
proj4string(e1) <-proj4string(states)
plot(e1)
states_48 <- crop(states, e1)
states_48$val <- seq(1:49)
states_48$val

states_48_colors = terrain.colors(nrow(states_48), alpha = 1)[order(states_48$val)]
states_48_colors1 <- heat.colors(nrow(states_48))[order(states_48$val)]
states_48_colors2 <- blue_pal(nrow(states_48))[order(states_48$val)]

#some of the options doesnt color two states
#plot(states_48, col = states_48_colors, cex = (states_48$val)/2,  main = "lower 48 states")
plot(states_48, col = states_48_colors1, main = "Lower 48 states")
#multicolor
plot(states_48, col = order(states_48$val), main = "Lower 48 states")
#works
plot(states_48, col = states_48_colors2, main = "Lower 48 states of USA")

#plot(crop(states, e1))

EPSG <- make_EPSG()

#using lambers conformal conic
states_48a<- states_48

states_48a<-spTransform(states_48, CRS("+proj=lcc +lat_1=30 +lat_2=60 +lat_0=38 +lon_0=126 +datum=WGS84" ))     
proj4string(states_48a)
plot(states_48a, col = states_48_colors2, main = "Lower 48 states of USA")

states_48b<-spTransform(states_48, CRS("+init=epsg:28992"))#CRS("+init=epsg:32760")

states_48b<-spTransform(states_48, CRS("+proj=moll"))#has the 2nd most distortion 
states_48b<-spTransform(states_48, CRS("+proj=wintri"))
states_48b<-spTransform(states_48, CRS("+proj=laea"))#has the most distortion

png(filename = "lower48_distorted.png", width=8, height=6,units = 'in', res=300)
plot(states_48b, col = states_48_colors2, main = "Lower 48 states of USA - distorted", border = "red", lty = 2, lwd = 3)
dev.off()
proj4string(states_48b)
proj4string(states_48)
##Q3
#####My data - PM2.5
rots_atr = read.csv(paste0(book_data, "rots_atr.csv"))
class(rots_atr)
summary(rots_atr)
#specify coordinates to transfrom dataframe to spatialpointsdataframe
coordinates(rots_atr) = ~ Lon + Lat
crs(rots_atr) = "+proj=tmerc +lat_0=4.666666666666667 +lon_0=-1 +k=0.99975 +x_0=274319.5100000002
+y_0=0 +ellps=clrk80 +units=ft +no_defs "
class(rots_atr)
summary(rots_atr)
View(rots_atr@coords)
head(rots_atr@coords)
plot(rots_atr)

#read in accra shapefile
gama<- readOGR("C:/Users/aalli/Desktop/ACCRA/GAMA CORRECTED/Gama_disCODE.shp")
proj4string(gama)
gama_proj <- proj4string(gama)
summary(gama)
plot(gama)
head(gama@data)


proj4string(rots_atr) = proj4string(gama)
summary(rots_atr)
plot(rots_atr, add = T)
plot(rots_atr)
plot(gama, add= TRUE)

#color trick
rots_colors = terrain.colors(nrow(rots_atr), alpha = 1)[order(rots_atr$avgPM)]
rots_colors1 = terrain.colors(nrow(rots_atr), alpha = 1)[order(rots_atr$site_type)]

#rots_colors2 <- blue_pal(nrow(rots_atr))[order(rots_atr$site_type)]

plot(rots_atr, add = TRUE, col = rots_col1, cex = rots.cex, pch = 16, 
     main = "PM2.5 concentrations")

legend("bottomright",legend = rots.cex, fill = rots_colors)

plot(gama)
plot(rots_atr, col = rots.col, cex = 2, pch = 16, 
     main = "site_type", add = T)

proj4string(rots_atr) #= proj4string(gama)
lo_rots <- spTransform(rots_atr, crs(gama))
proj4string(lo_rots)

plot(lo_rots, add = T, pch = 16, col = rots.col, cex = 2)
##project points
rots_coords<- rots_atr1[, c("long", "lat")]
data_rots<-rots_atr1[, 1:8]
crs.rots<-CRS("+init=epsg:4326")  ## to project to latlong WGS84
rots_point<-SpatialPointsDataFrame(coords = rots_coords, data = data_rots,
                                 proj4string = crs.rots)
class(rots_point)
proj4string(rots_point)
rots_utm <- spTransform(rots_point, CRS("+init=epsg:32630")) 
proj4string(rots_utm)
proj4string(gama)
rots.cex <- ifelse(rots_atr$avgPM <15, 1, 
                 ifelse(rots_atr$avgPM >= 15 & rots_atr$avgPM <25, 2,
                        ifelse(rots_atr$avgPM >= 25 & rots_atr$avgPM <35, 3,
                               ifelse(rots_atr$avgPM >35 , 4, NA)))) 

#specify color of bubble
rots.col <- ifelse(rots_atr$site_type== "traffic", "blue", 
                 ifelse(rots_atr$site_type == "low-dens", "seagreen",
                        ifelse(rots_atr$site_type == "high-dens", "red", 
                               ifelse(rots_atr$site_type == "other", "gold1",
                                      ifelse(rots_atr$site_type == "commercial", "purple",NA))))) 

#specify color of bubble
rots_col1 <- ifelse(rots_atr$avgPM <15, "gold1", 
                   ifelse(rots_atr$avgPM >= 15 & rots_atr$avgPM <25, "seagreen",
                          ifelse(rots_atr$avgPM >= 25 & rots_atr$avgPM <35, "blue",
                                 ifelse(rots_atr$avgPM >35 , "red", NA))))

gama_utm <- spTransform(gama, CRS("+init=epsg:32630")) # Re-project GAMA to UTM 
proj4string(gama_utm)

png("PM2.5.png", width = 8, height = 6, units = "in", res = 300)
plot(gama_utm, col = "grey85", border = "transparent",lwd = 0.2, main = "PM2.5 concentrations")
#add points
plot(rots_utm, add=TRUE, pch=16, cex = 2,  col= rots_col1)
locator(1)
legend(x = 829365.9, y = 612836.8, legend=c("<15", "15 - 25", "26-35", ">35"), 
       title = c(expression(bold("PM" [2.5]* " (µg/m" ^3* ")"))), 
       title.adj = 0,
       col = c("gold1", "seagreen", "blue", "red"), 
       pch = c(16, 16, 16, 16), 
       pt.cex = 2, 
       bty="n", cex=1,
       x.intersp = .8,
       y.intersp = .8,
       xpd = TRUE)
dev.off()

png("site type.png", width = 8, height = 6, units = "in", res = 300)
plot(gama_utm, col = "grey85", border = "transparent",lwd = 0.2, main = "Site type")
#add points
plot(rots_utm, add=TRUE, pch=16, cex = 2, col= rots.col)
locator(1)
legend(x = 828137, y = 616907.5, legend=c("traffic", "low density", "high density", "background", "commercial"), 
       title = expression(bold("Site type")), 
       title.adj = 0,
       col = c("blue", "seagreen", "red", "gold1", "purple"), 
       pch = c(16, 16, 16, 16, 16), 
       pt.cex = 2, 
       bty="n", cex=1,
       x.intersp = 0.8,
       y.intersp = 0.8, #text. = 12,
       xpd = TRUE)
dev.off()

temp2<- over(rots_atr, geometry(gama))

###############################
#join rots_utm to gama_utm by location
temp2<- over(rots_utm, gama_utm)
temp2a<- over(gama_utm, rots_utm, fn = mean)


#########################
poly1<-gama[1,]
crop(rots_atr, poly1)
################
require(spatialEco)
pts.poly <- point.in.poly(rots_atr, gama)
head(pts.poly@data)

a.data <- over(rots_atr, gama[,"GAMA_Dis_C"])
a.data<-over(rots_atr, gama@polygons[1:12],  fn = mean)
a.data<-over(gama, rots_atr)
###############
rots_atr = read.csv(paste0(book_data, "rots_atr.csv"))
coordinates(rots_atr) = ~ Lon + Lat
a.data<- over(gama@polygons[1:12], rots_atr@coords, fn = mean)

# Convert to sf-objects
srdf.sf <- st_as_sf(gama)
meuse.sf <- st_as_sf(rots_atr)

# Keep all "meuse.sf", sort by row.names(meuse.sf). Default overlay is "intersects".
meuse_srdf <- st_join(meuse.sf, srdf.sf) 

# Keeps all "srdf.sf", sort by row.names(srdf.sf)
srdf_meuse <- st_join(srdf.sf, meuse.sf)

# Convert back to Spatial*
meuse_srdf <- as(meuse_srdf, "Spatial")
srdf_meuse <- as(srdf_meuse, "Spatial")

plot(gama, col = rots_atr$medPM)
plot(srdf_meuse, col = rots_colors, cex = (rots_atr$avgPM)/20, pch = 16, 
     main = "PM2.5 concentrations")

gama_col <- blue_pal(nrow(gama))[order(rots_atr$avgPM)]
plot(gama, col = gama_col)

###package sf
#geom_poly
gama1<-map_data(gama)
rots_atr1 = read.csv(paste0(book_data, "rots_atr.csv")) 
rots_atr1<- rots_atr1 %>% rename(long = Lon, lat = Lat)
gama2<-merge(gama1, rots_atr1)

gama2<-join(gama1, )
ggplot(gama1)+
  geom_polygon(aes(x=long, y= lat, group = group, fill = region), color = "black", lwd = 1)+
  geom_col(data = rots_atr1, aes(avgPM))
  theme_void()+guides(fill = F)

