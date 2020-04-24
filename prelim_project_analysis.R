
# land use and gama shapefiles --------------------------------------------
library(rgdal)
library(raster)
setwd("C:/Users/aalli/Desktop/ACCRA")
gama_land <- readOGR("GAMA_WB_Landuse category/land_gama.shp")
proj4string(gama_land)
crs(gama_land)
plot(gama_land)
##project to UTM
gama_land_utm <- spTransform(gama_land, CRS("+init=epsg:32630")) # Re-project AMA to UTM 
crs(gama_land_utm)
##set color for land use types
levels(gama_land_utm$gridcode)
gridcolor <- c("lightblue1","greenyellow","steelblue1","gold3")[gama_land_utm$gridcode]
gridcolor
plot(gama_land_utm, col = gridcolor, border = "transparent")

##load sites data
setwd("C:/Users/aalli/Desktop/Log form")
fixed_coords <- readxl::read_excel("fixed_coords.xlsx")
all_sites <- read.csv("C:/Users/aalli/Desktop/Log form/PM_allsites.csv")
crs(all_sites)
#specify coordinates to transfrom dataframe to spatialpointsdataframe
coordinates(all_sites) = ~ Lon + Lat
crs(all_sites)<-CRS("+init=epsg:4326")
all_sites_utm <- spTransform(all_sites,  CRS("+init=epsg:32630"))

plot(all_sites_utm, add = TRUE, pch = sites.cex, col = sites.col)
sites.pch <- ifelse(all_sites_utm$ID %in% c("AD", "Ash", "EL", "N1West", "La", "JT", 
                              "Nima", "Taifa", "UGH","TMW"), 17, 16)
sites.col<- ifelse(all_sites_utm$ID %in% c("AD", "Ash", "EL", "N1West", "La", "JT", 
                  "Nima", "Taifa", "UGH","TMW"), "violetred1", "black")
###plot
png(filename = "study area.png", width=8, height=6,units = 'in', res=300)
plot(gama_land_utm, col = gridcolor, border = "transparent")
plot(all_sites_utm, add = TRUE, pch = sites.pch, col = sites.col)

locator(1)
legend(x = 819591, y = 612399.7, legend=c("Fixed sites", "Rotating sites",
       "Low density residential", "High density residential", "Commercial",
       "Other (Forest, Agricultural, etc.)"), 
       col = c("violetred1", "black", "lightblue1","greenyellow","steelblue1","gold3"), 
       pch = c(17, 16, 46, 46, 46, 46), 
       pt.cex = c(1, 1, 14, 14, 14, 14), 
       bty="n", cex=1,
       x.intersp = .8,
       y.intersp = .8,
       xpd = TRUE)
dev.off()

library(prettymapr)
addscalebar(plotunit = "km", plotepsg = 32630, style = "bar", unitcategory = "metric", 
            pos = "bottomleft", widthhint = 0.1, htin = 0.04 , padin = c(1, 0.05) )

####didnt work
legend(x = 829365.9, y = 612836.8, 
       legend=c(levels(gama_land_utm$gridcode[gama_land_utm$gridcode])), 
       fill = gridcolor[gama_land_utm$gridcode], 
       bty="n", cex=1,
       x.intersp = .8,
       y.intersp = .8,
       xpd = TRUE)

#####plot of PM spatial variation 
png(filename = "pm spatial vary.png", width=8, height=6,units = 'in', res=300)
plot(gama_utm)

pm_sites<- read.csv("C:/Users/aalli/Desktop/Log form/PM_allsites.csv")

pm_sites1<- pm_sites %>% group_by(ID, site_type, Lon, Lat) %>% summarise(mass.conc = mean(concn))

#specify coordinates to transfrom dataframe to spatialpointsdataframe
coordinates(pm_sites1) = ~ Lon + Lat
crs(pm_sites1)<-CRS("+init=epsg:4326")
pm_sites_utm <- spTransform(pm_sites1,  CRS("+init=epsg:32630"))

##bubble size
pm.cex <- ifelse(pm_sites_utm$mass.conc <15, 1.5, 
                 ifelse(pm_sites_utm$mass.conc >= 15 & pm_sites_utm$mass.conc <25, 2.5,
                        ifelse(pm_sites_utm$mass.conc >= 25 & pm_sites_utm$mass.conc <35, 3.5,
                               ifelse(pm_sites_utm$mass.conc >35 , 5.0, NA)))) 
#specify color of bubble
pm.col <- ifelse(pm_sites_utm$site_type== "traffic", "red", 
                ifelse(pm_sites_utm$site_type == "low-dens", "seagreen",
                       ifelse(pm_sites_utm$site_type == "high-dens", "blue", 
                              ifelse(pm_sites_utm$site_type == "other", "gold1",
                                     ifelse(pm_sites_utm$site_type == "commercial", "black", NA)))))
pm.col1<- adjustcolor(pm.col, alpha.f = 0.5)

##plot
png(filename = "pm spatial vary.png", width=8, height=6,units = 'in', res=300)
plot(gama_utm)
plot(pm_sites_utm, add = TRUE, cex = pm.cex, col = pm.col1, pch = 16)
View(pm_sites_utm)
#legend
locator(1)
legend1 <- legend(x = 829365.9, y = 612836.8, 
                  legend=c("Traffic", "High density residential", "Low density residential", "Background", "Commercial"), 
                  title =  expression(bold("Site type")),
                  title.adj = 0,
                  col = c("red", "blue",  "seagreen", "gold1", "black"), 
                  pch = c(16, 16, 16, 16, 16), 
                  pt.cex = c(1, 1, 1, 1), 
                  bty="n", cex=.6, 
                  xpd = TRUE)

legend3 <- legend(x = 812703.7, y = 609744.1,
                  legend=c("<15", "15 - 25", "26-35", ">35"), 
                  title = c(expression(bold("PM" [2.5]* " (µg/m" ^3* ")"))), 
                  title.adj = 0,
                  col = c("gray52", "gray52", "gray52", "gray52"), 
                  pch = c(16, 16, 16, 16), 
                  pt.cex = c(0.7, 1, 2, 2.5), 
                  bty="n", cex=.5,
                  x.intersp = 1.3,
                  y.intersp = 1.8,
                  xpd = TRUE)

dev.off()

#######################################################################################################

#####LUR Project##################
#st_as_sf ("input sp object") function is used to transform sp object to sf object
#as_Spatial("input sf object") function is used to transform sf object to sp object
#then add projection, if lost 
##using st_crs and project functions respectively

library(tidyverse)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
require(gstat)
library(geoR) 
library(maptools)
library(sf)
######load pm and weather data 
pmw.dat <- read.csv("C:/Users/aalli/Desktop/R CODES ANALYSIS/Data visualization/pmw_data.csv")

plot(pmw.dat)
pmw.dat<-pmw.dat %>% drop_na(concn)

#specify coordinates to transfrom dataframe to spatialpointsdataframe
coordinates(pmw.dat) = ~ Lon + Lat
crs(pmw.dat)<-CRS("+init=epsg:4326")
pmw.dat.utm <- spTransform(pmw.dat,  CRS("+init=epsg:32630"))
proj4string(pmw.dat.utm)
proj4string(gama_utm)
View(pmw.dat.utm)
#plot
plot(gama_utm)
plot(pmw.dat.utm, add = T)


###load road shapefile
library(rgdal)
road <- readOGR("C:/Users/aalli/Desktop/ACCRA/data/Miscellaneous/Road networks/GAR-road-OpenStreetMap-April 2018/GAR_Road_OSM.shp", layer = "GAR_Road_OSM")
table(road$fclass)

####filter by fclass == motorway, trunk and primary == major roads
road1 <- road[road$fclass == "motorway" | road$fclass == "trunk" | 
                      road$fclass == "primary", ]

unique(road$fclass)

proj4string(road1)
levels(road1$fclass)
table(road1$fclass)

##reproject to utm
road.utm <- spTransform(road1, CRS("+init=epsg:32630")) 
proj4string(road.utm)
plot(road.utm)
##clip road1 to gama extent
road.utm1<-gIntersection(gama_utm, road.utm, byid = T)# , not sure whtether to intersect by id

##clip road1 to gama extent
library(rgeos)
road.int<-gIntersection(gama_utm, road.utm) ##gama_utm
road.int
plot(gama_utm)
plot(road.int, add = T, col = 6)

##50m buffer around site points
#puts a circular buffer around each individual points, calc distance to road within 50m radius
pmdat_50 <- gBuffer(pmw.dat.utm, width=50, byid=TRUE)
plot(pmdat_50, col = "red")
pmdat_100 <- gBuffer(pmw.dat.utm, width=100, byid=TRUE)
pmdat_200 <- gBuffer(pmw.dat.utm, width=200, byid=TRUE)
pmdat_500 <- gBuffer(pmw.dat.utm, width=500, byid=TRUE)
pmdat_1000 <- gBuffer(pmw.dat.utm, width=1000, byid=TRUE)
pmdat_2000 <- gBuffer(pmw.dat.utm, width=2000, byid=TRUE)
pmdat_4000 <- gBuffer(pmw.dat.utm, width=4000, byid=TRUE)


#############################################################################################################
#######################################USE SF AND RASTER PACKAGES FOR ALL ###################
#############################################################################################################


#####LUR Project##################
#st_as_sf ("input sp object") function is used to transform sp object to sf object
#as_Spatial("input sf object") function is used to transform sf object to sp object
#then add projection, if lost 
##using st_crs and project functions respectively
library(tidyverse)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
require(gstat)
library(geoR) 
library(maptools)
library(sf)





###distance of points to major road with sf package - didnt produce desired result####
road.sf <- st_as_sf(road.int)
dis2rd_sf <- apply(st_distance(pm_dat, road.sf1), 2, min)
head(dis2rd_sf)
View(dis2rd_sf)
plot(pm_dat$dis2rd_gdis, pm_dat$dis2rd_sf)


#plot(road.utm, col = 6, add = T)
#points(pmw.dat.utm)
#plot(pmdat_50, add = T)
#plot(pmdat_100, add = T)
#plot(pmdat_1000, add = T)
#text(pmw.dat.utm, labels=pmw.dat.utm$ID,data=pmw.dat.utm, cex=0.9, font=2, col = "red")

# rasterize road data and calculate distance to road ----------------------

#####rasterize road data
#https://rdrr.io/cran/raster/man/distance.html
#create blank raster
blank_raster <- raster(nrow = 100, ncol = 100, extent(pmw.dat.utm))
#first rasterize the road data
road.rast <- rasterize(road.utm, blank_raster)
plot(road.rast)
#then calculate distance to each cell
dist.road <- distance(road.rast)
plot(distance(road.rast))

#extract distance to points
dist.pts <- raster::extract(dist.road, pmw.dat.utm )
dist.pts
##add to points data
pmw.dat.utm$dis2rd_rast <- as.numeric(dist.pts)

##scatterplot
plot(pmw.dat.utm$dis2rd_rast, pmw.dat.utm$conc.log)
plot(pmw.dat.utm$dis2rd_gdis, pmw.dat.utm$conc.log)
plot(pmw.dat.utm$dis2rd_gdis, pmw.dat.utm$dis2rd_rast)##compare road.rast vs road 4rm gdistance

#next do a regression with distance
reg <- lm(log(concn)~dis2rd_gdis, data = pmw.dat.utm)
summary(reg)

plot(pmw.dat.utm$concn, pmw.dat.utm$dis2rd_gdis)

# road length with rgeos package ------------------------------------------
###########  length of roads within each buffer ##############
#function for calculating length - buffer,intersect, extract
road_50 <- gIntersection(pmdat_50, road.utm1, byid = T)
road.50.lnt<- rgeos::gLength(road_50, byid = T)

road.50 <- rgeos::gLength(gIntersection(pmdat_50, road.utm1, byid = T))

View(road.50.lnt)
road.50.lnt <-data.frame(road.50.lnt)
#remove column.number from row.name
row.names(road.50.lnt) <- sub("[ ].*", "", row.names(road.50.lnt))
#temp$ID<-NULL
road.50.lnt$X = row.names(road.50.lnt)
View(road.50.lnt)
#merge by x- row name
pmw.dat.utm<- merge(pmw.dat.utm, road.50.lnt)
View(pmw.dat.utm)
#pmw.dat.utm$rd_50 <-NULL
###
road_100 <- gIntersection(pmdat_100, road.utm1, byid = T)
road.100.lnt<- rgeos::gLength(road_100, byid = T)
road.100.lnt



###############################################################################################
###########################################################################################################


#####distance of points to major road with rgeos package  ------------------------------------ 
pmw.dat.utm$dist.maj.rd<-apply(gDistance(pmw.dat.utm, road.utm,byid=TRUE),2,min)
View(dis2rd)

##distance to major road - sf package
dis2rd_sf <- apply(st_distance(pm_dat, road.sf), 2, min)

###distance of points to major road with sf package ####
maj.rd$id <-1:nrow(maj.rd)
head(maj.rd$id)
site.near.rd <- st_join(dat.pm1, maj.rd, join = st_nearest_feature)
st_crs(site.near.rd)
plot(site.near.rd$geometry)
dist.maj.rd <- st_distance(dat.pm1, maj.rd[site.near.rd$id, ], by_element = T, which = "Euclidean")

dist.maj.rd <- st_distance(dat.pm.utm, maj.rd, by_element = T, which = "Euclidean")
View(dist.maj.rd)
head(dist.maj.rd)


##convert data to sf objects
#points
pm_dat <- st_as_sf(pmw.dat.utm)
#road
road.sf <-st_as_sf(road.utm)# previous major road
#crs.gama <- "+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

osm_road <- read_sf("C:/Users/aalli/Desktop/ACCRA/data/Miscellaneous/Road networks/GAR-road-OpenStreetMap-April 2018/GAR_Road_OSM.shp", layer = "GAR_Road_OSM")
table(osm_road$fclass)
st_crs(osm_road)
osm_road_utm <- st_transform(osm_road, crs.gama)
st_crs(osm_road_utm)
unique(osm_road_utm$fclass)

maj.rd <- osm_road_utm %>% filter(fclass %in% c("motorway", "trunk", "primary", "secondary", "tertiary"))
med.rd <- osm_road_utm %>% filter(fclass %in% c("motorway_link", "trunk_link", "primary_link", "secondary_link"))
min.rd <- osm_road_utm %>% filter(fclass %in% c("unclassified", "residential", "living_street", "pedestrian"))

####road length for each buffer size - 50,100,200,500,1000,2000,4000

#50m buffer
#intersect road and buffer of points
#rd.50 <- st_intersection(road.sf, pm_50) #83 obs
rd.50 <- st_intersection(road.sf, st_buffer(pm_dat, dist = 50))
plot(rd.50)

#length of each line segment
rd.50 <- rd.50 %>%  mutate(rlt.50 = st_length(rd.50))
class(rd.50$rlt.50)
rd.50$rlt.50 <- as.numeric(rd.50$rlt.50)

#sum of line segments per ID
rd.50 <- rd.50 %>% group_by(ID, Start_date) %>% 
        summarise(length.50 = sum(rlt.50)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.50)

class(rd.50)

#join with points data
pm_dat <- merge(pm_dat, rd.50, all = T)

#dat.rd.merge <- merge(pm_dat, rd.50, all = T)
#View(dat.rd.merge)

#scatterplot
plot(pm_dat$road.50.lnt, pm_dat$length.50)
#plot(dat.rd.merge$road.50.lnt, dat.rd.merge$length.50)


#100m buffer
#intersect road and buffer of points
rd.100 <- st_intersection(road.sf, st_buffer(pm_dat, dist = 100))
plot(rd.100[, 2])

#length of each line segment
rd.100 <- rd.100 %>%  mutate(rlt.100 = st_length(rd.100))
class(rd.100$rlt.100)
rd.100$rlt.100 <- as.numeric(rd.100$rlt.100)

#sum of line segments per ID
rd.100 <- rd.100 %>% group_by(ID, Start_date) %>% 
        summarise(length.100 = sum(rlt.100)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.100)

class(rd.100)

#join with points data
pm_dat <- merge(pm_dat, rd.100, all = T)
View(pm_dat)

#200m buffer
#intersect road and buffer of points
rd.200 <- st_intersection(road.sf, st_buffer(pm_dat, dist = 200))
plot(rd.200[, 2])

#length of each line segment
rd.200 <- rd.200 %>%  mutate(rlt.200 = st_length(rd.200))
class(rd.200$rlt.200)
rd.200$rlt.200 <- as.numeric(rd.200$rlt.200)

#sum of line segments per ID
rd.200 <- rd.200 %>% group_by(ID, Start_date) %>% 
        summarise(length.200 = sum(rlt.200)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.200)

class(rd.200)

#join with points data
pm_dat <- merge(pm_dat, rd.200, all = T)
View(pm_dat)

#500m buffer
#intersect road and buffer of points
rd.500 <- st_intersection(road.sf, st_buffer(pm_dat, dist = 500))
plot(rd.500[, 2])

#length of each line segment
rd.500 <- rd.500 %>%  mutate(rlt.500 = st_length(rd.500))
class(rd.500$rlt.500)
rd.500$rlt.500 <- as.numeric(rd.500$rlt.500)

#sum of line segments per ID
rd.500 <- rd.500 %>% group_by(ID, Start_date) %>% 
        summarise(length.500 = sum(rlt.500)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.500)

class(rd.500)

#join with points data
pm_dat <- merge(pm_dat, rd.500, all = T)
View(pm_dat)

#1000m buffer
#intersect road and buffer of points
rd.1000 <- st_intersection(road.sf, st_buffer(pm_dat, dist = 1000))
#plot(rd.1000[, 2])

#length of each line segment
rd.1000 <- rd.1000 %>%  mutate(rlt.1000 = st_length(rd.1000))
class(rd.1000$rlt.1000)
rd.1000$rlt.1000 <- as.numeric(rd.1000$rlt.1000)

#sum of line segments per ID
rd.1000 <- rd.1000 %>% group_by(ID, Start_date) %>% 
        summarise(length.1000 = sum(rlt.1000)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.1000)

class(rd.1000)

#join with points data
pm_dat <- merge(pm_dat, rd.1000, all = T)
View(pm_dat)


#2000m buffer
#intersect road and buffer of points
rd.2000 <- st_intersection(road.sf, st_buffer(pm_dat, dist = 2000))
#plot(rd.2000[, 2])

#length of each line segment
rd.2000 <- rd.2000 %>%  mutate(rlt.2000 = st_length(rd.2000))
class(rd.2000$rlt.2000)
rd.2000$rlt.2000 <- as.numeric(rd.2000$rlt.2000)

#sum of line segments per ID
rd.2000 <- rd.2000 %>% group_by(ID, Start_date) %>% 
        summarise(length.2000 = sum(rlt.2000)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.2000)

class(rd.2000)

#join with points data
pm_dat <- merge(pm_dat, rd.2000, all = T)
View(pm_dat)

#4000m buffer
#intersect road and buffer of points
rd.4000 <- st_intersection(road.sf, st_buffer(pm_dat, dist = 4000))
#plot(rd.4000[, 2])

#length of each line segment
rd.4000 <- rd.4000 %>%  mutate(rlt.4000 = st_length(rd.4000))
class(rd.4000$rlt.4000)
rd.4000$rlt.4000 <- as.numeric(rd.4000$rlt.4000)

#sum of line segments per ID
rd.4000 <- rd.4000 %>% group_by(ID, Start_date) %>% 
        summarise(length.4000 = sum(rlt.4000)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.4000)

class(rd.4000)

#join with points data
pm_dat <- merge(pm_dat, rd.4000, all = T)
View(pm_dat)


# NDVI --------------------------------------------------------------------
#https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/vegetation-indices-NDVI-in-R/
#https://rspatial.org/rs/rs.pdf
#
#vegetation indices function
vi <- function(img, k, i) {
        bk <- img[[k]]
        bi <- img[[i]]
        vi <- (bk - bi) / (bk + bi)
        return(vi)
}

#load landsat data
library(raster)
setwd("C:/Users/aalli/Desktop/ACCRA")
raslist <- paste0('NDVI/LC08_L1TP_193056_20200219_20200219_01_RT_B', 1:11, ".TIF")
raslist
#stack raster list
landsat <- stack(raslist[c(1:7, 9:11)])#band8 has a larger extent
plot(landsat)
# For Landsat NIR = 5, red = 4.
ndvi <- vi(landsat, 5, 4)
View(ndvi)
plot(ndvi, col = rev(terrain.colors(10)), main = "Landsat-NDVI")
#hist(log(ndvi+1), xlab = "NDVI Index Value", main = "NDVI: Distribution of pixels") #log transformed data was also skewed

##extract ndvi for points (sites)
ndvi.pts <- raster::extract(ndvi, pmw.dat.utm)
View(ndvi.pts)
str(ndvi.pts)

hist(ndvi.pts, xlab = "NDVI Index Value", main = "NDVI: Distribution of pixels")

##add to points data
pmw.dat.utm$ndvi <- as.numeric(ndvi.pts)

####NDVI with points buffer
#50m buffer
#ndvi.50 <- raster::extract(ndvi, pm_50, fun = median)##fun = sum
ndvi.50 <- raster::extract(ndvi, st_buffer(pm_dat, dist = 50), fun = median)##fun = sum
pm_dat$ndvi.50 <- as.numeric(ndvi.50)
View(pm_dat)


#100m buffer
ndvi.100 <- raster::extract(ndvi, st_buffer(pm_dat, dist = 100), fun = median)
pm_dat$ndvi.100 <- as.numeric(ndvi.100)
View(pm_dat)

#200m buffer
ndvi.200 <- raster::extract(ndvi, st_buffer(pm_dat, dist = 200), fun = median)
pm_dat$ndvi.200 <- as.numeric(ndvi.200)
View(pm_dat)

#500m buffer
ndvi.500 <- raster::extract(ndvi, st_buffer(pm_dat, dist = 500), fun = median)##fun = sum
pm_dat$ndvi.500 <- as.numeric(ndvi.500)
View(pm_dat)


#1000m buffer
ndvi.1000 <- raster::extract(ndvi, st_buffer(pm_dat, dist = 1000), fun = median)##fun = sum
pm_dat$ndvi.1000 <- as.numeric(ndvi.1000)
View(pm_dat)

#2000m buffer
ndvi.2000 <- raster::extract(ndvi, st_buffer(pm_dat, dist = 2000), fun = median)##fun = sum
pm_dat$ndvi.2000 <- as.numeric(ndvi.2000)
View(pm_dat)

#4000m buffer
ndvi.4000 <- raster::extract(ndvi, st_buffer(pm_dat, dist = 4000), fun = median)##fun = sum
pm_dat$ndvi.4000 <- as.numeric(ndvi.4000)
View(pm_dat)

#extract( r , y = 1:ncell(r) , buffer = 3000 , fun = sum )

#######land use###########
plot(gama_land_utm, col = gridcolor, border = "gray99")
levels(gama_land_utm$gridcode)
gridcolor <- c("lightblue1","greenyellow","steelblue1","gold3")[gama_land_utm$gridcode]
#1 - low dens, 2- high-dens, 3- comm, 4- Other

landcov<- raster('data/Accra_2014_final_classification.img')
#head(landcov@data@attributes)

landcov1 <- crop(landcov, extent(gama_utm)) ##crop to gama to make extract faster
plot(pm_dat, col = 6, add = T)
plot(landcov1)
landcov1
dim(landcov1)##number of celss(x,y)
proj4string(landcov1) = proj4string(pmw.dat.utm) ##make projection similar
plot(landcov1)
plot(pmw.dat.utm, add = T)

##check each land use category
length(landcov1[landcov1== 1])
length(landcov1[landcov1== 2])
length(landcov1[landcov1== 3])
length(landcov1[landcov1== 4])

#extract raster attributes related to buffers of points data
##50m buffer
#luc1 <- raster::extract(landcov1, st_as_sf(dat.rd.merge), buffer = 50)
#temp <- raster::extract(landcov1, st_buffer(pm_dat, dist = 50))
luc1 <- raster::extract(landcov1, st_buffer(pm_dat, dist = 50))
View(luc1)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.50 <- t(sapply(luc1, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.50)
class(ldc.50)
colnames(ldc.50)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
names(ldc.50)
##add to points data
pm_dat$ld.area.50 <- ldc.50[,"low-dens"]
pm_dat$hd.area.50 <- ldc.50[,"high-dens"]
pm_dat$com.area.50 <- ldc.50[,"commercial"]
pm_dat$other.area.50 <- ldc.50[,"other"]
str(pm_dat)

##100m buffer
#luc2 <- raster::extract(landcov1, pm_100)
luc2 <- raster::extract(landcov1, st_buffer(pm_dat, dist = 100))
View(luc1)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.100 <- t(sapply(luc2, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.100) 
colnames(ldc.100)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
pm_dat$ld.area.100 <- ldc.100[,"low-dens"]
pm_dat$hd.area.100 <- ldc.100[,"high-dens"]
pm_dat$com.area.100 <- ldc.100[,"commercial"]
pm_dat$other.area.100 <- ldc.100[,"other"]
str(pm_dat)
View(pm_dat)

##200m buffer
luc3 <- raster::extract(landcov1, st_buffer(pm_dat, dist = 200))
View(luc1)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.200 <- t(sapply(luc3, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.200) 
colnames(ldc.200)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
pm_dat$ld.area.200 <- ldc.200[,"low-dens"]
pm_dat$hd.area.200 <- ldc.200[,"high-dens"]
pm_dat$com.area.200 <- ldc.200[,"commercial"]
pm_dat$other.area.200 <- ldc.200[,"other"]
str(pm_dat)
View(pm_dat)


##500m buffer
luc4 <- raster::extract(landcov1, st_buffer(pm_dat, dist = 500))
View(luc4)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.500 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.500) 
colnames(ldc.500)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
pm_dat$ld.area.500 <- ldc.500[,"low-dens"]
pm_dat$hd.area.500 <- ldc.500[,"high-dens"]
pm_dat$com.area.500 <- ldc.500[,"commercial"]
pm_dat$other.area.500 <- ldc.500[,"other"]
str(pm_dat)
View(pm_dat)

##1000m buffer
luc5 <- raster::extract(landcov1, st_buffer(pm_dat, dist = 1000))
View(luc5)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.1000 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.1000) 
colnames(ldc.1000)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
pm_dat$ld.area.1000 <- ldc.1000[,"low-dens"]
pm_dat$hd.area.1000 <- ldc.1000[,"high-dens"]
pm_dat$com.area.1000 <- ldc.1000[,"commercial"]
pm_dat$other.area.1000 <- ldc.1000[,"other"]
str(pm_dat)
View(pm_dat)

##2000m buffer
luc6 <- raster::extract(landcov1, st_buffer(pm_dat, dist = 2000))
View(luc6)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.2000 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.2000) 
colnames(ldc.2000)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
pm_dat$ld.area.2000 <- ldc.2000[,"low-dens"]
pm_dat$hd.area.2000 <- ldc.2000[,"high-dens"]
pm_dat$com.area.2000 <- ldc.2000[,"commercial"]
pm_dat$other.area.2000 <- ldc.2000[,"other"]
str(pm_dat)
View(pm_dat)

##4000m buffer
luc7 <- raster::extract(landcov1, st_buffer(pm_dat, dist = 4000))
View(luc7)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.4000 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.4000) 
colnames(ldc.4000)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
pm_dat$ld.area.4000 <- ldc.4000[,"low-dens"]
pm_dat$hd.area.4000 <- ldc.4000[,"high-dens"]
pm_dat$com.area.4000 <- ldc.4000[,"commercial"]
pm_dat$other.area.4000 <- ldc.4000[,"other"]
str(pm_dat)
View(pm_dat)

#########population density##############
pop.den<- raster('data/Gridded population rasters/GHS population grid -250m - 2015/pop_GAR_2015_GHS1.tif')
pop.den
plot(pop.den)
proj4string(pop.den)
#project to utm
pop.den.utm<- projectRaster(pop.den, crs = crs(pmw.dat.utm))
proj4string(pop.den.utm)
hist(pop.den.utm, xlab = "Population density", main = "Distribution of Population density")
hist(log(pop.den.utm), xlab = "Population density", main = "Distribution of Population density")
plot(pop.den.utm)
plot(pmw.dat.utm, add = T)
pop.den.utm

#extract population density attributes related to buffers of points data
##50m buffer
temp<-raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 50))
pop.50 <- raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 50), fun = median)
View(pop.50)
pm_dat$pop.50 <- as.numeric(pop.50)

##100m buffer
pop.100 <- raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 100), fun = median)
View(pop.100)
pm_dat$pop.100 <- as.numeric(pop.100)

##200m buffer
pop.200 <- raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 200), fun = median)
View(pop.200)
pm_dat$pop.200 <- as.numeric(pop.200)

##500m buffer
pop.500 <- raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 500), fun = median)
View(pop.500)
pm_dat$pop.500 <- as.numeric(pop.500)

##1000m buffer
pop.1000 <- raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 1000), fun = median)
View(pop.1000)
pm_dat$pop.1000 <- as.numeric(pop.1000)

##2000m buffer
pop.2000 <- raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 2000), fun = median)
View(pop.2000)
pm_dat$pop.2000 <- as.numeric(pop.2000)

##4000m buffer
pop.4000 <- raster::extract(pop.den.utm, st_buffer(pm_dat, dist = 4000), fun = median)
View(pop.4000)
pm_dat$pop.4000 <- as.numeric(pop.4000)


######building density
library(sf)
build.den <- read_sf("data/Gama_buildings_merged_wgs84/Gama_buildings_merged_wgs84.shp")
View(build.den)
plot(st_geometry(build.den))
build.den
plot(build.den)
crs(build.den)
#project to utm
#build.den.utm <- st_transform(build.den, crs = st_crs("+init=epsg:32630"))
build.den.utm <- st_transform(build.den, 32630)
crs(build.den.utm)
View(build.den.utm)
build.den.utm

plot(st_geometry(build.den.utm[1:10,]))
head(build.den.utm)

plot(st_geometry(build.den.utm[1,]))
plot(st_geometry(build.den.utm[2,]))


#intersect building density and buffer - 50
build.count.50 <- lengths(st_intersects(st_buffer(pm_dat, dist = 50), build.den.utm))
pm_dat$build.count.50 <- as.numeric(build.count.50)

#intersect building density and buffer - 100
build.count.100 <- lengths(st_intersects(st_buffer(pm_dat, dist = 100), build.den.utm))
pm_dat$build.count.100 <- as.numeric(build.count.100)

#intersect building density and buffer - 200
build.count.200 <- lengths(st_intersects(st_buffer(pm_dat, dist = 200), build.den.utm))
pm_dat$build.count.200 <- as.numeric(build.count.200)

#intersect building density and buffer - 500
build.count.500 <- lengths(st_intersects(st_buffer(pm_dat, dist = 500), build.den.utm))
pm_dat$build.count.500 <- as.numeric(build.count.500)

#intersect building density and buffer - 1000
build.count.1000 <- lengths(st_intersects(st_buffer(pm_dat, dist = 1000), build.den.utm))
pm_dat$build.count.1000 <- as.numeric(build.count.1000)

#intersect building density and buffer - 2000
build.count.2000 <- lengths(st_intersects(st_buffer(pm_dat, dist = 2000), build.den.utm))
pm_dat$build.count.2000 <- as.numeric(build.count.2000)

#intersect building density and buffer - 4000
build.count.4000 <- lengths(st_intersects(st_buffer(pm_dat, dist = 4000), build.den.utm))
pm_dat$build.count.4000 <- as.numeric(build.count.4000)

#distance to airport
airport <- read_sf("data/gha_airports/GHA_Airports.shp")
crs(airport)
plot(airport$geometry)
airp_utm <- st_transform(airport, crs.gama)
crs(airp_utm)
airp <- st_crop(airp_utm ,pm_dat)
air_dist <- st_distance(pm_dat, airp)
View(air_dist)

plot(pm_dat$geometry)
plot(airp$geometry, add = T, col = 6)

#add to data
pm_dat$dist_airport <- as.numeric(air_dist)
summary(air_dist)
View(pm_dat)


###START HERE###############
#number of bus stops within buffer
transport <- read_sf("data/OSM 2019/ghana-latest-free-2019.shp/gis_osm_transport_free_1.shp")
plot(transport$geometry)
crs(transport)
transport_utm <- st_transform(transport, crs.gama)
crs(transport_utm)
trans_gama <- st_crop(transport_utm ,st_as_sf(gama_utm))
rm(bus.stp)
#temp <- st_crop(transport_utm ,pm_dat)
#plot(temp$geometry)
plot(trans_gama$geometry)
plot(pm_dat, col = 6, add = T)
plot(gama_utm, add = T)
unique(trans_gama$fclass)
bus.stp <- trans_gama %>% filter(fclass %in% "bus_stop")
View(bus.stp)
unique(bus.stp$fclass)

#50m buffer
bus.count <- lengths(st_intersects(st_buffer(pm_dat, dist = 50), bus.stp))
View(bus.count)
str(bus.count)
#pm_dat$bus.count.50 <- as.numeric(bus.count)
pm_dat$bus.count.50 <- bus.count

#100m buffer
bus.count1 <- lengths(st_intersects(st_buffer(pm_dat, dist = 100), bus.stp))
View(bus.count1)
pm_dat$bus.count.100 <- bus.count1

#200m buffer
bus.count2 <- lengths(st_intersects(st_buffer(pm_dat, dist = 200), bus.stp))
View(bus.count2)
pm_dat$bus.count.200 <- bus.count2

#500m buffer
bus.count3 <- lengths(st_intersects(st_buffer(pm_dat, dist = 500), bus.stp))
View(bus.count3)
pm_dat$bus.count.500 <- bus.count3

#1000m buffer
bus.count4 <- lengths(st_intersects(st_buffer(pm_dat, dist = 1000), bus.stp))
View(bus.count4)
pm_dat$bus.count.1000 <- bus.count4

#2000m buffer
bus.count5 <- lengths(st_intersects(st_buffer(pm_dat, dist = 2000), bus.stp))
View(bus.count5)
pm_dat$bus.count.2000 <- bus.count5

#4000m buffer
bus.count6 <- lengths(st_intersects(st_buffer(pm_dat, dist = 4000), bus.stp))
View(bus.count6)
pm_dat$bus.count.4000 <- bus.count6

##number of bus stations within buffer
unique(trans_gama$fclass)
bus.stat<- trans_gama %>% filter(fclass %in% "bus_station")
View(bus.stat)
plot(bus.stat$geometry)
plot(pm_dat, col = 6, add = T)
plot(gama_utm, add = T)

#50m buffer
bus.stat.count <- lengths(st_intersects(st_buffer(pm_dat, dist = 50), bus.stat))
View(bus.stat.count)
summary(bus.stat.count)
pm_dat$bus.stat.50 <- bus.stat.count

#100m buffer
bus.stat.count1 <- lengths(st_intersects(st_buffer(pm_dat, dist = 100), bus.stat))
View(bus.stat.count1)
summary(bus.stat.count1)
pm_dat$bus.stat.100 <- bus.stat.count1

#200m buffer
bus.stat.count2 <- lengths(st_intersects(st_buffer(pm_dat, dist = 200), bus.stat))
View(bus.stat.count2)
summary(bus.stat.count2)
pm_dat$bus.stat.200 <- bus.stat.count2

#500m buffer
bus.stat.count3 <- lengths(st_intersects(st_buffer(pm_dat, dist = 500), bus.stat))
View(bus.stat.count3)
summary(bus.stat.count3)
pm_dat$bus.stat.500 <- bus.stat.count3

#1000m buffer
bus.stat.count4 <- lengths(st_intersects(st_buffer(pm_dat, dist = 1000), bus.stat))
View(bus.stat.count4)
summary(bus.stat.count4)
pm_dat$bus.stat.1000 <- bus.stat.count4

#2000m buffer
bus.stat.count5 <- lengths(st_intersects(st_buffer(pm_dat, dist = 2000), bus.stat))
View(bus.stat.count5)
summary(bus.stat.count5)
pm_dat$bus.stat.2000 <- bus.stat.count5

#4000m buffer
bus.stat.count6 <- lengths(st_intersects(st_buffer(pm_dat, dist = 4000), bus.stat))
View(bus.stat.count6)
summary(bus.stat.count6)
pm_dat$bus.stat.4000 <- bus.stat.count6


##send data to mike
#writeRaster(road.utm1, "C:/Users/aalli/Desktop/major_roads.tiff")
road.utm2 <- SpatialLinesDataFrame(road.utm1, data.frame(id = 1:length(road.utm1)))
writeOGR(road.utm2, dsn="C:/Users/aalli/Documents/ECO 697/LUR project/data", "major_roads",driver = "ESRI Shapefile")
writeOGR(pmw.dat.utm, dsn="C:/Users/aalli/Documents/ECO 697/LUR project/data", "sites_data",driver = "ESRI Shapefile")



##kriging/semivariance#####
###https://rstudio-pubs-static.s3.amazonaws.com/159648_374fc566d7014172b1825e4ed72a390b.html
###https://cran.r-project.org/web/packages/kriging/kriging.pdf
###https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf
###http://gsp.humboldt.edu/OLM/R/04_01_Variograms.html

require(gstat)
hist(pmw.dat.utm$concn, breaks = 10, xlab="PM2.5 (Ug/m^3)", main="Histogram of PM2.5")
summary(pmw.dat.utm$concn)

#the histogram shows the data is highly skewed, using log-transformed value is recommended.
pmw.dat.utm$conc.log <- log10(pmw.dat.utm$concn)
hist(pmw.dat.utm$conc.log , breaks = 10,xlab="Log10 PM2.5 (Ug/m^3)",main="Histogram of Log-transformed PM2.5 (Ug/m^3)")




# Semivariance for Log-transformed Data -----------------------------------
#the fundamental concept is spatial auto-correlation: a variable can be correlated to 
#itself, with the strength of correlation depending on seperation distance 
#(and possibly direction); this should be evident as a relation between seperation 
#distance and correlation (also expresssed as semi-variance);
#each pair of observation has a semivariance,defined as gamma(xi,xj)=1/2[z(xi)-z(xj)]^2, 
#where x is location, z is concentration.
#compute number of point-pairs
n<-length(pmw.dat.utm$conc.log)
n*(n-1)/2
## [1] 123753
#how the semivariances related to the seperation distance: using empirical variogram 
#defined as average semivariance within some seperation ranges (distances) get variogram 
#(semivariance vs distance) of log-tranformed data
#v<-variogram(conc.log~1, pmw.dat.utm)##regress conc.log with itself
# the 1 here represents the spatial mean of the dependent variable - that is conc.log is only depedent on itself.
v<-variogram(conc.log~1, pmw.dat.utm) 
print(plot(v,plot.numbers=T))

head(v)
#plot.numbers show how many point-pairs in each bin (distance interval); 
#more point-pairs are more reliable. The bin interval can be specified using width in the
#variogram function. If there are not enough points (at least 100 to 150) with a good 
#distribution of distance intervals, it is not possible to model a variogram,and kriging 
#should not be used.
#fitting the curve using vgm function for different models
#sph model
vm <-vgm(psill=1000, model="Sph", range=10000, nugget=1000)
print(plot(v,pl=T,model=vm))

###statistically estimated variogram
fit.v <- fit.variogram(v[-1,], vm)
print(plot(v,pl=T,model=fit.v))
fit.v
##use values from sph row in fit.v
vm<-vgm(psill=0.01631598, model="Sph", range=9180.485, nugget=1000)
vm<-vgm(psill=0.01631598, model="Nug", range=0, nugget=1)
print(plot(v,pl=T,model=vm))


#the sill, which is the semivariance value at the plateau of the variogram, representing the
#distance where there is no correlation between the observations
#range means at what distance between point-pairs that there is no more spatial dependence. 
#Nugget means the semivariance at zero seperation. Total sill means the semivariance at the
#range and beyond. kappa is smoothness parameter for the Matern class of variogram models
#there are many fit models to choose from, using vgm() or print(show.vgms()) to view the 
#avaiable models. Selecting a model form based on knowledge of assumed spatial process. 
#A common model for soil attributes is the spherical model.

#adjusting the fit
vmf<-fit.variogram(v,vm)
print(plot(v,pl=T,model=vmf))


plot(gama_utm)
plot(road.utm1, col = "6", add = T)
points(pm_sites_utm,cex = pm.cex, col = pm.col1, pch = 16)
#points(pmw.dat.utm,cex = pm.cex, col = pm.col1, pch = 16)

####variogram of PM2.5 and distance to road
v1<-variogram(conc.log~dis2rd_gdis, pmw.dat.utm)

# residual variogram w.r.t. a linear trend:
distmat <- as.matrix(dist(pmw.dat.utm@coords))
#maximum distance to consider in correlogram/variogram
maxdist <- 1/2 * max(distmat)


v1<-variogram(concn~Lon+Lat, pmw.dat.utm, cutoff = maxdist)
v1
plot(v1)

# "Exp" (for exponential)
# "Sph" (for spherical)
# "Gau" (for Gaussian)
# "Mat" (for Matern)
vm1 <-vgm(psill=1000, model="Sph", range=10000, nugget=1000)
plot(v1, vm1)
print(plot(v1,pl=T,model=vm1))
##fitted model - gives values for psill and range
fit.v1 <- fit.variogram(v1[-1,], vm1)
fit.v1
plot(fit.v1, fit.v1)
print(plot(v1,pl=T,model=fit.v1, main = "semivariogram of coordinates"))
summary(fit.v1)

##use values from fit.v1
vm2 <-vgm(psill=569.2319, model="Sph", range=9110.77, nugget=1000)
plot(v1, vm2)

##excluding the point-pairs at zero, the range of correlation is btwn 2278 and 9111, 
##to keep focus on deriving local estimates, choose the lower bound of the range (~2000) as the max. dist. for buffer



# other codes -------------------------------------------------------------


#setwd("C:/Users/aalli/Desktop/ACCRA/data/Miscellaneous/Road networks/Other")
#road <- readOGR("GAR-road-WorldBank-2009-MajorRoadsOnly/Ghana_Roads.shp")
#file.exists("GAR-road-WorldBank-2009-MajorRoadsOnly/Ghana_Roads.shp")
###################################################################
landcov2 <- ratify(landcov1)
landcov2
levels(landcov2)
View(landcov2)

##subset by land cover types
low_dens <- landcov2[landcov2@data@attributes[[1]] == 1, ]

landcov3 <- extract(landcov2, pmdat_50)
landcov3

low_dens <- mask(landcov2, landcov2 == 1, maskvalue = FALSE)

lowd_50 <- extract(low_dens, pmdat_50)
table(lowd_50[])
View(lowd_50)
#tapply(area(landcov3), landcov3[], sum)
View(landcov3)
temp3<-  aggregate(landcov3, fact=2, fun=function(vals, na.rm) {
        sum(vals==1, na.rm=na.rm)
})

cov_pct <- lapply(unique(r), function(land_class) {
        aggregate(r, fact=2, fun=function(vals, na.rm) {
                sum(vals==land_class, na.rm=na.rm)/length(vals)
        })
})

lowd_50 <- extract(landcov2[landcov2@data@attributes[[1]]$ID == 1], pmdat_50)
raster::subset(landcov2, landcov2@data@attributes[[1]] == 1)

temp <-landcov2[landcov2 ==1]
temp <- landcov2[landcov2@data@attributes[[1]]$ID == 1]

low_dens <- mask(landcov2, landcov2$ID == 1, maskvalue = F)

low_dens <- landcov[landcov==1]
low_dens

temp <- readGDAL('data/gama_landclass_wb.tif')
temp1 <- raster(temp)
temp1
temp1[temp1, temp1@data@values == 1]
#temp <- raster('data/gama_landclass_wb.tif')
View(temp)
str(temp)
temp

lowd_50 <- extract(low_dens, pmdat_50)
temp <- area(lowd_50)
aggregate(getValues(area(lowd_50, weights=FALSE)), by=list(getValues(lowd_50)), sum)
############
summ.land <- lapply(landcov2, function(x){
        prop.table(table(x))})
summ.land



# library(tiff)
# readTIFF("accra_landcover.tif")
# library(png)
# img <- readTIFF(system.file("img", "gama_landclass_wb.tif", package="tiff"))
# img<-readPNG(system.file("img", "Accra_2014_final_classification.img", package="png"))
# ## convert it to a raster, interpolate =F to select only sample of pixels of img
# img.r <- as.raster(img,interpolate=F)
# str(img.r)
# writeRaster(landcov, "C:/Users/aalli/Desktop/landcov.tiff")
#writeOGR(pmdat_50, dsn="C:/Users/aalli/Desktop", "pmdat_50.shp",driver = "ESRI Shapefile")
# landcov1<- raster('C:/Users/aalli/Desktop/landcov.tif')


# spatial overlay
join = st_join(pm_50, int, by = "Strt_dt")

join = st_join(rd1.utm, int)
View(join)
##add len to points
join2 <- st_join(pm_50, join, by = "Strt_dt")
View(join2)
# use the ID of the polygon for the aggregation
out = group_by(join, Id.x) %>%
        summarize(length = sum(len))

tapply(st_length(ints), ints$osm_id, sum)
# #unlist(raslist)
# extent(raslist) <- c(720585, 947415, 524085, 755715)
# extent(raslist[8]) ==extent(720585, 947415, 524085, 755715)
# 
# raslist1 <- unlist(raslist)
# extent(raslist1[8] ,720585, 947415, 524085, 755715)

# band8<-stack(raslist[8])
# band9 <- stack(raslist[9])
# extent(band9)
# extent(band8)
# ##make extent of band8 equal to other bands
# extent(band8)<-extent(720585, 947415, 524085, 755715)
# extent(band8)
# extent(band9)

# ras8 <- raster(vals=values(band8),ext=extent(band9),crs=crs(band9),
#               nrows=dim(band9)[1],ncols=dim(band9)[2])
#pmw.dat.utm$road.50.lnt <-
#pmw.dat.utm$road_50 <- t(gDistance(pmdat_50, road.utm1, byid = T))
#pmw.dat.utm$road_100 <- t(gDistance(pmdat_100, road.utm1, byid = T))
#pmw.dat.utm$road_1000 <- t(gDistance(pmdat_1000, road.utm1, byid = T))
#dist.road.rast1 <-gridDistance(road.rast, pmw.dat.utm@coords)
#dist.road.rast1 <-distanceFromPoints(road.rast, pmw.dat.utm@coords)

#road.utm2 <- SpatialLinesDataFrame(road.utm1, data.frame(id = 1:length(road.utm1))) ##to write shapefile convert to SLDF
#writeOGR(pmdat_50, dsn="C:/Users/aalli/Desktop", "pmdat_50.shp",driver = "ESRI Shapefile")
#read data
#road
#rd<- read_sf("C:/Users/aalli/Desktop/ACCRA/data/Miscellaneous/Road networks/GAR-road-OpenStreetMap-April 2018/GAR_Road_OSM.shp")
#rd1 <- rd[rd$fclass == "motorway" | rd$fclass == "trunk" | rd$fclass == "primary", ]
#st_crs(rd1)
#st_crs(pm_dat)
#rd1.utm <- st_transform(rd1,  crs.gama)#crs("+init=epsg:32630")
##points
#pm_dat <- read_sf("C:/Users/aalli/Desktop/pmdat_50.shp")

# road length with sf package ---------------------------------------------

####road length uisng sf package################
#https://gis.stackexchange.com/questions/280760/intersecting-lines-and-polygons-and-calculating-line-length-in-r
library(sf)
library(maptools)
class(road.utm1)

##convert data to sf objects
#points
pm_dat <- st_as_sf(pmw.dat.utm)
#road
rd1 <-st_as_sf(road.int)
#crs.gama <- "+init=epsg:32630 +proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
st_crs(pm_dat)
st_crs(rd1)

## 50m buffer for points
pm_50 <- st_buffer(pm_dat, dist = 50)
plot(pm_50)

#intersect road and buffer of points
#50m buffer
#rd.int <- st_intersection(rd1, pm_50)

rd.int <- st_intersection(rd1, pm_50) 
#lenght of each line segment
rd.int <- rd.int %>%  mutate(rd.length = st_length(rd.int))
class(rd.int$rd.length)
rd.int$rd.length <- as.numeric(rd.int$rd.length)
#sum of line segments per ID
rd.int <- rd.int %>% group_by(ID, Start_date) %>% 
        summarise(length.50 = sum(rd.length)) %>% ungroup() %>% 
        as.data.frame() %>% select(ID, Start_date, length.50)

class(rd.int)

#join with points data
pm_dat <- merge(pm_dat, rd.int, all = T)

plot(pm_dat)

##scatter plot of road.50.lnt(rgeos) and length.50 (sf)
plot(pm_dat$road.50.lnt, pm_dat$length.50)

str(temp$rlt.50)

View(rd.int)

##sum distances per ID
library(tidyverse)
#using sf package
rd.int1 <- rd.int %>% select(ID, length.50, Strt_dt) %>% group_by(ID, Strt_dt) %>% summarise(sum_len = sum(length.50))
View(rd.int1)

# rd_len_rgeos <- pmw.dat.utm %>% select(ID, road.50.lnt)
# rd_len_rgeos <- as.data.frame (pmw.dat.utm)[, c("Start_date","ID","road.50.lnt")]
# View(rd_len_rgeos)

##road length using rgeos package
rd_len_rgeos1 <- pm_dat %>% drop_na(road.50.lnt) %>% group_by(ID, Start_date) %>% 
        summarise(sum_len = sum(road.50.lnt))
View(rd_len_rgeos1)
#####distance to road with sf package
#pmw.dat.utm$dis2rd_gdis <- as.numeric(t(gDistance(pmw.dat.utm, road.int, byid = T )))# byid = T 
#View(dis2rd_gdis)
#View(pmw.dat.utm)
##distance to road 
road.sf <- st_as_sf(road.int)
road.sf1 <-st_as_sf(road.utm)

pm_dat$dist_maj_rd <- as.numeric(st_distance(pm_dat, road.sf))

View(road.sf1)
sites_near_road <- st_join(pm_dat, road.sf1, join = st_nearest_feature)
rd_dist_temp <- st_distance(pm_dat, road.sf1[sites_near_road$osm_id, ], by_element = TRUE)
View(rd_dist_temp)
#use sf package for road length
pm_dat <- st_as_sf(pmw.dat.utm) 
names(pm_dat)
#pm_dat[, c(1:14, )]

#intersect building density and buffer - 200

# build.50 = st_intersection(build.den.utm, st_buffer(pm_dat, dist = 50))
# View(build.50)
# #count of buliding within buffer
# build.count.50 <- build.50 %>% group_by(ID, Start_date) %>% summarise(build.50 = n()) %>% 
#         ungroup %>% as.data.frame() %>% select(ID, Start_date, build.50) %>% 
#         rename(Start_date = Start_date)
# View(build.count.50)
# # #merge
# dat.rd.merge <- merge(dat.rd.merge, build.count.50, by = c("ID", "Start_date"), all = T)
# which(is.na(dat.rd.merge$build.50))

build.100 = st_intersection(build.den.utm, pm_100)
View(build.100)
#count of buliding within buffer
build.count.100 <- build.100 %>% group_by(ID, Start_date) %>% summarise(build.100 = n()) %>% 
        ungroup %>% as.data.frame() %>% select(ID, Start_date, build.100)
View(build.count.100)
# #merge
dat.rd.merge <- merge(dat.rd.merge, build.count.100, by = c("ID", "Start_date"), all = T)
which(is.na(dat.rd.merge$build.100))

build.200 = st_intersection(build.den.utm, pm_200)
View(build.200)
#count of buliding within buffer
build.count.200 <- build.200 %>% group_by(ID, Start_date) %>% summarise(build.200 = n()) %>% 
        ungroup %>% as.data.frame() %>% select(ID, Start_date, build.200)
View(build.count.200)
# #merge
dat.rd.merge <- merge(dat.rd.merge, build.count.200, by = c("ID", "Start_date"), all = T)
which(is.na(dat.rd.merge$build.200))


#intersect building density and buffer
build.500 = st_intersection(build.den.utm, pm_500) ##error message - shapefile intersects with itself
View(build.500)
#count of buliding within buffer
build.count.500 <- build.500 %>% group_by(ID, Start_date) %>% summarise(build.500 = n()) %>% 
        ungroup %>% as.data.frame() %>% select(ID, Start_date, build.500)
View(build.count.500)
# #merge
dat.rd.merge <- merge(dat.rd.merge, build.count.500, by = c("ID", "Start_date"), all = T)
which(is.na(dat.rd.merge$build.500))

#intersect building density and buffer
build.1000 = st_intersection(build.den.utm, pm_1000) ##error message - shapefile intersects with itself
View(build.1000)
#count of buliding within buffer
build.count.1000 <- build.1000 %>% group_by(ID, Start_date) %>% summarise(build.1000 = n()) %>% 
        ungroup %>% as.data.frame() %>% select(ID, Start_date, build.1000)
View(build.count.1000)
# #merge
dat.rd.merge <- merge(dat.rd.merge, build.count.1000, by = c("ID", "Start_date"), all = T)
which(is.na(dat.rd.merge$build.1000))

###repeat for other buffer sizes
#use to verify that merge works first!
#which(is.na(dat.rd.merge$build.50))
dat.rd.merge1 <- merge(dat.rd.merge, build.count.50, by = c("ID", "Start_date"), all = T)
dat.rd.merge1 %>% is.na(build.50)

#change data to sf object
dat.rd.merge1 <- st_as_sf(dat.rd.merge)
class(dat.rd.merge1)
plot(dat.rd.merge1$geometry)
#temp <- aggregate(id ~ name, data = sam.int, FUN = length)


##buffer sizes
pm_50 <- st_buffer(pm_dat, dist = 50)
pm_100 <-  st_buffer(pm_dat, dist = 100)
pm_200 <-  st_buffer(pm_dat, dist = 200)
pm_500 <-  st_buffer(pm_dat, dist = 500)
pm_1000 <-  st_buffer(pm_dat, dist = 1000)
pm_2000 <-  st_buffer(pm_dat, dist = 2000)
pm_4000 <-  st_buffer(pm_dat, dist = 4000)


# pm_50 <-  st_as_sf(pmdat_50)
# pm_100 <-  st_as_sf(pmdat_100)
# pm_200 <-  st_as_sf(pmdat_200)
# pm_500 <-  st_as_sf(pmdat_500)
# pm_1000 <-  st_as_sf(pmdat_1000)
# pm_2000 <-  st_as_sf(pmdat_2000)
# pm_4000 <-  st_as_sf(pmdat_4000)

class(pm_50)
class(pm_1000)
plot(pm_50[, 2])
plot(pm_1000[, 2])