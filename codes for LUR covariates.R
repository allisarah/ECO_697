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

#######################################USE SF AND RASTER PACKAGES FOR ALL ###################

#pm data
dat.pm <- read.csv("C:/Users/aalli/Desktop/R CODES ANALYSIS/Data visualization/pmw_data.csv")
dat.pm <- dat.pm %>% drop_na(concn)

#specify coordinates to transfrom dataframe to spatialpointsdataframe
coordinates(dat.pm) = ~ Lon + Lat
crs(dat.pm)<-CRS("+init=epsg:4326")
dat.pm.utm <- spTransform(dat.pm,  CRS("+init=epsg:32630"))

dat.pm1 <- st_as_sf(dat.pm.utm)
extent(dat.pm1)
class(dat.pm1)
#dat.pm.utm <- st_as_sf(dat.pm, coords = c("Lon", "Lat"), crs = 32360) - returns a very samll extent

#road shapefile
osm_road <- read_sf("C:/Users/aalli/Desktop/ACCRA/data/Miscellaneous/Road networks/GAR-road-OpenStreetMap-April 2018/GAR_Road_OSM.shp", layer = "GAR_Road_OSM")
table(osm_road$fclass)
st_crs(osm_road)
osm_road_utm <- st_transform(osm_road, crs.gama)
st_crs(osm_road_utm)
unique(osm_road_utm$fclass)
plot(osm_road$geometry)

maj.rd <- osm_road_utm %>% filter(fclass %in% c("motorway", "trunk", "primary", "secondary", "tertiary"))
med.rd <- osm_road_utm %>% filter(fclass %in% c("motorway_link", "trunk_link", "primary_link", "secondary_link"))
min.rd <- osm_road_utm %>% filter(fclass %in% c("unclassified", "residential", "living_street", "pedestrian"))


###distance to road - st_dist - sf package
#distance to nearest road
osm_road_utm$id <-1:nrow(osm_road_utm)
head(osm_road_utm$id)
site.near.roads <- st_join(dat.pm1, osm_road_utm, join = st_nearest_feature)
st_crs(site.near.roads)
plot(site.near.roads$geometry)
dist.roads <- st_distance(dat.pm1, osm_road_utm[site.near.roads$id, ], by_element = T, which = "Euclidean")
View(dist.roads)
head(dist.roads)
dat.pm1$dist.roads <- as.numeric(dist.roads)

##inverse distance to major road
dat.pm1$inv.roads<- 1/(dat.pm1$dist.roads)

##inverse distance to major road squared
dat.pm1$invsq.roads <- (dat.pm1$inv.roads)^2



###distance of points to nearest major road with sf package 
maj.rd$id <-1:nrow(maj.rd)
head(maj.rd$id)
site.near.rd <- st_join(dat.pm1, maj.rd, join = st_nearest_feature)
st_crs(site.near.rd)
plot(site.near.rd$geometry)
dist.maj.rd <- st_distance(dat.pm1, maj.rd[site.near.rd$id, ], by_element = T, which = "Euclidean")
View(dist.maj.rd)
dat.pm1$dist.maj.rd <- as.numeric(dist.maj.rd)

##inverse distance to major road
dat.pm1$inv.maj.rd <- 1/(dat.pm1$dist.maj.rd)

##inverse distance to major road squared
dat.pm1$invsq.maj.rd <- (dat.pm1$inv.maj.rd)^2

library(mapview)
mapview(dat.pm1[3, ]) + maj.rd[site.near.rd$id[3], ]

###distance of points to medium road with sf package 
med.rd$id <-1:nrow(med.rd)
head(med.rd$id)
site.near.medrd <- st_join(dat.pm1, med.rd, join = st_nearest_feature)
st_crs(site.near.medrd)
plot(site.near.medrd$geometry)
dist.med.rd <- st_distance(dat.pm1, med.rd[site.near.medrd$id, ], by_element = T, which = "Euclidean")
View(dist.med.rd)
head(dist.med.rd)
dat.pm1$dist.med.rd <- as.numeric(dist.med.rd)

##inverse distance to medium road
dat.pm1$inv.med.rd <- 1/(dat.pm1$dist.med.rd)

##inverse distance to medium road squared
dat.pm1$invsq.med.rd <- (dat.pm1$inv.med.rd)^2

###distance of points to minor road with sf package 
min.rd$id <-1:nrow(min.rd)
head(min.rd$id)
site.near.minrd <- st_join(dat.pm1, min.rd, join = st_nearest_feature)
st_crs(site.near.minrd)
plot(site.near.minrd$geometry)
dist.min.rd <- st_distance(dat.pm1, min.rd[site.near.minrd$id, ], by_element = T, which = "Euclidean")
View(dist.min.rd)
head(dist.min.rd)
dat.pm1$dist.min.rd <- as.numeric(dist.min.rd)

##inverse distance to minor road
dat.pm1$inv.min.rd <- 1/(dat.pm1$dist.min.rd)

##inverse distance to minor road squared
dat.pm1$invsq.min.rd <- (dat.pm1$inv.min.rd)^2

#cor(log(dat.pm1$concn), log(dat.pm1$invsq.maj.rd)) - 0.16
#summary(lm(log(concn)~invsq.maj.rd, data = dat.pm1))
#cor(log(dat.pm1$concn), log(dat.pm1$dist.med.rd))


#total road length - all roads
#50m buffer
#intersect roads and buffer of points
all.roads.50 <- st_intersection(osm_road_utm, st_buffer(dat.pm1, dist = 50))
head(all.roads.50)
plot(all.roads.50$geometry, col = 6)
#length of each line segment
roads.50 <- all.roads.50 %>% mutate(length.50 = st_length(all.roads.50)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(roads.50.length = sum(length.50)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, roads.50.length)
head(roads.50)
View(roads.50)

roads.50$roads.50.length <- as.numeric(roads.50$roads.50.length)

#join with points data
dat.pm1 <- merge(dat.pm1, roads.50, all = T)

#100m buffer
#intersect roads and buffer of points
all.roads.100 <- st_intersection(osm_road_utm, st_buffer(dat.pm1, dist = 100))
head(all.roads.100)
plot(all.roads.100$geometry, col = 6)
#length of each line segment
roads.100 <- all.roads.100 %>% mutate(length.100 = st_length(all.roads.100)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(roads.100.length = sum(length.100)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, roads.100.length)
head(roads.100)
View(roads.100)

roads.100$roads.100.length <- as.numeric(roads.100$roads.100.length)

#join with points data
dat.pm1 <- merge(dat.pm1, roads.100, all = T)

#200m buffer
#intersect roads and buffer of points
all.roads.200 <- st_intersection(osm_road_utm, st_buffer(dat.pm1, dist = 200))
head(all.roads.200)
plot(all.roads.200$geometry, col = 6)
#length of each line segment
roads.200 <- all.roads.200 %>% mutate(length.200 = st_length(all.roads.200)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(roads.200.length = sum(length.200)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, roads.200.length)
head(roads.200)
View(roads.200)

roads.200$roads.200.length <- as.numeric(roads.200$roads.200.length)

#join with points data
dat.pm1 <- merge(dat.pm1, roads.200, all = T)


#500m buffer
#intersect roads and buffer of points
all.roads.500 <- st_intersection(osm_road_utm, st_buffer(dat.pm1, dist = 500))
head(all.roads.500)
plot(all.roads.500$geometry, col = 6)
#length of each line segment
roads.500 <- all.roads.500 %>% mutate(length.500 = st_length(all.roads.500)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(roads.500.length = sum(length.500)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, roads.500.length)
head(roads.500)
View(roads.500)

roads.500$roads.500.length <- as.numeric(roads.500$roads.500.length)

#join with points data
dat.pm1 <- merge(dat.pm1, roads.500, all = T)


#1000m buffer
#intersect roads and buffer of points
all.roads.1000 <- st_intersection(osm_road_utm, st_buffer(dat.pm1, dist = 1000))
head(all.roads.1000)
plot(all.roads.1000$geometry, col = 6)
#length of each line segment
roads.1000 <- all.roads.1000 %>% mutate(length.1000 = st_length(all.roads.1000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(roads.1000.length = sum(length.1000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, roads.1000.length)
head(roads.1000)
View(roads.1000)

roads.1000$roads.1000.length <- as.numeric(roads.1000$roads.1000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, roads.1000, all = T)


#2000m buffer
#intersect roads and buffer of points
all.roads.2000 <- st_intersection(osm_road_utm, st_buffer(dat.pm1, dist = 2000))
head(all.roads.2000)
plot(all.roads.2000$geometry, col = 6)
#length of each line segment
roads.2000 <- all.roads.2000 %>% mutate(length.2000 = st_length(all.roads.2000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(roads.2000.length = sum(length.2000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, roads.2000.length)
head(roads.2000)
View(roads.2000)

roads.2000$roads.2000.length <- as.numeric(roads.2000$roads.2000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, roads.2000, all = T)


#4000m buffer
#intersect roads and buffer of points
all.roads.4000 <- st_intersection(osm_road_utm, st_buffer(dat.pm1, dist = 4000))
head(all.roads.4000)
plot(all.roads.4000$geometry, col = 6)
#length of each line segment
roads.4000 <- all.roads.4000 %>% mutate(length.4000 = st_length(all.roads.4000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(roads.4000.length = sum(length.4000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, roads.4000.length)
head(roads.4000)
View(roads.4000)

roads.4000$roads.4000.length <- as.numeric(roads.4000$roads.4000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, roads.4000, all = T)


####road length for each buffer size - 50,100,200,500,1000,2000,4000
#major road
#50m buffer
#intersect road and buffer of points
maj.rd.50 <- st_intersection(maj.rd, st_buffer(dat.pm1, dist = 50))
head(maj.rd.50)
plot(maj.rd.50$geometry, col = 6)
#length of each line segment
maj.50 <- maj.rd.50 %>% mutate(length.50 = st_length(maj.rd.50)) %>% 
          group_by(ID, Start_date) %>% 
          summarise(maj.50.length = sum(length.50)) %>% ungroup() %>% 
          as.data.frame() %>% select(ID, Start_date, maj.50.length)
head(maj.50)
View(maj.50)

maj.50$maj.50.length <- as.numeric(maj.50$maj.50.length)

#join with points data
dat.pm1 <- merge(dat.pm1, maj.50, all = T)

#100m buffer
#intersect road and buffer of points
maj.rd.100 <- st_intersection(maj.rd, st_buffer(dat.pm1, dist = 100))
head(maj.rd.100)
plot(maj.rd.100$geometry, col = 6)
#length of each line segment
maj.100 <- maj.rd.100 %>% mutate(length.100 = st_length(maj.rd.100)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(maj.100.length = sum(length.100)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, maj.100.length)
head(maj.100)
View(maj.100)

maj.100$maj.100.length <- as.numeric(maj.100$maj.100.length)

#join with points data
dat.pm1 <- merge(dat.pm1, maj.100, all = T)

#200m buffer
#intersect road and buffer of points
maj.rd.200 <- st_intersection(maj.rd, st_buffer(dat.pm1, dist = 200))
head(maj.rd.200)
plot(maj.rd.200$geometry, col = 6)
#length of each line segment
maj.200 <- maj.rd.200 %>% mutate(length.200 = st_length(maj.rd.200)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(maj.200.length = sum(length.200)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, maj.200.length)
head(maj.200)
View(maj.200)

maj.200$maj.200.length <- as.numeric(maj.200$maj.200.length)

#join with points data
dat.pm1 <- merge(dat.pm1, maj.200, all = T)

#500m buffer
#intersect road and buffer of points
maj.rd.500 <- st_intersection(maj.rd, st_buffer(dat.pm1, dist = 500))
head(maj.rd.500)
plot(maj.rd.500$geometry, col = 6)
#length of each line segment
maj.500 <- maj.rd.500 %>% mutate(length.500 = st_length(maj.rd.500)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(maj.500.length = sum(length.500)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, maj.500.length)
head(maj.500)
View(maj.500)

maj.500$maj.500.length <- as.numeric(maj.500$maj.500.length)

#join with points data
dat.pm1 <- merge(dat.pm1, maj.500, all = T)

#1000m buffer
#intersect road and buffer of points
maj.rd.1000 <- st_intersection(maj.rd, st_buffer(dat.pm1, dist = 1000))
head(maj.rd.1000)
plot(maj.rd.1000$geometry, col = 6)
#length of each line segment
maj.1000 <- maj.rd.1000 %>% mutate(length.1000 = st_length(maj.rd.1000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(maj.1000.length = sum(length.1000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, maj.1000.length)
head(maj.1000)
View(maj.1000)

maj.1000$maj.1000.length <- as.numeric(maj.1000$maj.1000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, maj.1000, all = T)

#2000m buffer
#intersect road and buffer of points
maj.rd.2000 <- st_intersection(maj.rd, st_buffer(dat.pm1, dist = 2000))
head(maj.rd.2000)
plot(maj.rd.2000$geometry, col = 6)
#length of each line segment
maj.2000 <- maj.rd.2000 %>% mutate(length.2000 = st_length(maj.rd.2000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(maj.2000.length = sum(length.2000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, maj.2000.length)
head(maj.2000)
View(maj.2000)

maj.2000$maj.2000.length <- as.numeric(maj.2000$maj.2000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, maj.2000, all = T)

#4000m buffer
#intersect road and buffer of points
maj.rd.4000 <- st_intersection(maj.rd, st_buffer(dat.pm1, dist = 4000))
head(maj.rd.4000)
plot(maj.rd.4000$geometry, col = 6)
#length of each line segment
maj.4000 <- maj.rd.4000 %>% mutate(length.4000 = st_length(maj.rd.4000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(maj.4000.length = sum(length.4000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, maj.4000.length)
head(maj.4000)
View(maj.4000)

maj.4000$maj.4000.length <- as.numeric(maj.4000$maj.4000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, maj.4000, all = T)


#medium/secondary road
#50m buffer
#intersect road and buffer of points
med.rd.50 <- st_intersection(med.rd, st_buffer(dat.pm1, dist = 50))
head(med.rd.50)
plot(med.rd.50$geometry, col = 6)
#length of each line segment
med.50 <- med.rd.50 %>% mutate(length.50 = st_length(med.rd.50)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(med.50.length = sum(length.50)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, med.50.length)
head(med.50)
View(med.50)

med.50$med.50.length <- as.numeric(med.50$med.50.length)

#join with points data
dat.pm1 <- merge(dat.pm1, med.50, all = T)

#100m buffer
#intersect road and buffer of points
med.rd.100 <- st_intersection(med.rd, st_buffer(dat.pm1, dist = 100))
head(med.rd.100)
plot(med.rd.100$geometry, col = 6)
#length of each line segment
med.100 <- med.rd.100 %>% mutate(length.100 = st_length(med.rd.100)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(med.100.length = sum(length.100)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, med.100.length)
head(med.100)
View(med.100)

med.100$med.100.length <- as.numeric(med.100$med.100.length)

#join with points data
dat.pm1 <- merge(dat.pm1, med.100, all = T)

#200m buffer
#intersect road and buffer of points
med.rd.200 <- st_intersection(med.rd, st_buffer(dat.pm1, dist = 200))
head(med.rd.200)
plot(med.rd.200$geometry, col = 6)
#length of each line segment
med.200 <- med.rd.200 %>% mutate(length.200 = st_length(med.rd.200)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(med.200.length = sum(length.200)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, med.200.length)
head(med.200)
View(med.200)

med.200$med.200.length <- as.numeric(med.200$med.200.length)

#join with points data
dat.pm1 <- merge(dat.pm1, med.200, all = T)

#500m buffer
#intersect road and buffer of points
med.rd.500 <- st_intersection(med.rd, st_buffer(dat.pm1, dist = 500))
head(med.rd.500)
plot(med.rd.500$geometry, col = 6)
#length of each line segment
med.500 <- med.rd.500 %>% mutate(length.500 = st_length(med.rd.500)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(med.500.length = sum(length.500)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, med.500.length)
head(med.500)
View(med.500)

med.500$med.500.length <- as.numeric(med.500$med.500.length)

#join with points data
dat.pm1 <- merge(dat.pm1, med.500, all = T)

#1000m buffer
#intersect road and buffer of points
med.rd.1000 <- st_intersection(med.rd, st_buffer(dat.pm1, dist = 1000))
head(med.rd.1000)
plot(med.rd.1000$geometry, col = 6)
#length of each line segment
med.1000 <- med.rd.1000 %>% mutate(length.1000 = st_length(med.rd.1000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(med.1000.length = sum(length.1000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, med.1000.length)
head(med.1000)
View(med.1000)

med.1000$med.1000.length <- as.numeric(med.1000$med.1000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, med.1000, all = T)

#2000m buffer
#intersect road and buffer of points
med.rd.2000 <- st_intersection(med.rd, st_buffer(dat.pm1, dist = 2000))
head(med.rd.2000)
plot(med.rd.2000$geometry, col = 6)
#length of each line segment
med.2000 <- med.rd.2000 %>% mutate(length.2000 = st_length(med.rd.2000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(med.2000.length = sum(length.2000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, med.2000.length)
head(med.2000)
View(med.2000)

med.2000$med.2000.length <- as.numeric(med.2000$med.2000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, med.2000, all = T)

#4000m buffer
#intersect road and buffer of points
med.rd.4000 <- st_intersection(med.rd, st_buffer(dat.pm1, dist = 4000))
head(med.rd.4000)
plot(med.rd.4000$geometry, col = 6)
#length of each line segment
med.4000 <- med.rd.4000 %>% mutate(length.4000 = st_length(med.rd.4000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(med.4000.length = sum(length.4000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, med.4000.length)
head(med.4000)
View(med.4000)

med.4000$med.4000.length <- as.numeric(med.4000$med.4000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, med.4000, all = T)


#minor road
#50m buffer
#intersect road and buffer of points
min.rd.50 <- st_intersection(min.rd, st_buffer(dat.pm1, dist = 50))
head(min.rd.50)
plot(min.rd.50$geometry, col = 6)
#length of each line segment
min.50 <- min.rd.50 %>% mutate(length.50 = st_length(min.rd.50)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(min.50.length = sum(length.50)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, min.50.length)
head(min.50)
View(min.50)

min.50$min.50.length <- as.numeric(min.50$min.50.length)

#join with points data
dat.pm1 <- merge(dat.pm1, min.50, all = T)

#100m buffer
#intersect road and buffer of points
min.rd.100 <- st_intersection(min.rd, st_buffer(dat.pm1, dist = 100))
head(min.rd.100)
plot(min.rd.100$geometry, col = 6)
#length of each line segment
min.100 <- min.rd.100 %>% mutate(length.100 = st_length(min.rd.100)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(min.100.length = sum(length.100)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, min.100.length)
head(min.100)
View(min.100)

min.100$min.100.length <- as.numeric(min.100$min.100.length)

#join with points data
dat.pm1 <- merge(dat.pm1, min.100, all = T)

#200m buffer
#intersect road and buffer of points
min.rd.200 <- st_intersection(min.rd, st_buffer(dat.pm1, dist = 200))
head(min.rd.200)
plot(min.rd.200$geometry, col = 6)
#length of each line segment
min.200 <- min.rd.200 %>% mutate(length.200 = st_length(min.rd.200)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(min.200.length = sum(length.200)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, min.200.length)
head(min.200)
View(min.200)

min.200$min.200.length <- as.numeric(min.200$min.200.length)

#join with points data
dat.pm1 <- merge(dat.pm1, min.200, all = T)

#500m buffer
#intersect road and buffer of points
min.rd.500 <- st_intersection(min.rd, st_buffer(dat.pm1, dist = 500))
head(min.rd.500)
plot(min.rd.500$geometry, col = 6)
#length of each line segment
min.500 <- min.rd.500 %>% mutate(length.500 = st_length(min.rd.500)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(min.500.length = sum(length.500)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, min.500.length)
head(min.500)
View(min.500)

min.500$min.500.length <- as.numeric(min.500$min.500.length)

#join with points data
dat.pm1 <- merge(dat.pm1, min.500, all = T)

#1000m buffer
#intersect road and buffer of points
min.rd.1000 <- st_intersection(min.rd, st_buffer(dat.pm1, dist = 1000))
head(min.rd.1000)
plot(min.rd.1000$geometry, col = 6)
#length of each line segment
min.1000 <- min.rd.1000 %>% mutate(length.1000 = st_length(min.rd.1000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(min.1000.length = sum(length.1000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, min.1000.length)
head(min.1000)
View(min.1000)

min.1000$min.1000.length <- as.numeric(min.1000$min.1000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, min.1000, all = T)

#2000m buffer
#intersect road and buffer of points
min.rd.2000 <- st_intersection(min.rd, st_buffer(dat.pm1, dist = 2000))
head(min.rd.2000)
plot(min.rd.2000$geometry, col = 6)
#length of each line segment
min.2000 <- min.rd.2000 %>% mutate(length.2000 = st_length(min.rd.2000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(min.2000.length = sum(length.2000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, min.2000.length)
head(min.2000)
View(min.2000)

min.2000$min.2000.length <- as.numeric(min.2000$min.2000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, min.2000, all = T)

#4000m buffer
#intersect road and buffer of points
min.rd.4000 <- st_intersection(min.rd, st_buffer(dat.pm1, dist = 4000))
head(min.rd.4000)
plot(min.rd.4000$geometry, col = 6)
#length of each line segment
min.4000 <- min.rd.4000 %>% mutate(length.4000 = st_length(min.rd.4000)) %>% 
  group_by(ID, Start_date) %>% 
  summarise(min.4000.length = sum(length.4000)) %>% ungroup() %>% 
  as.data.frame() %>% select(ID, Start_date, min.4000.length)
head(min.4000)
View(min.4000)

min.4000$min.4000.length <- as.numeric(min.4000$min.4000.length)

#join with points data
dat.pm1 <- merge(dat.pm1, min.4000, all = T)


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

####NDVI with points buffer
#50m buffer
ndvi.50 <- raster::extract(ndvi, st_buffer(dat.pm1, dist = 50), fun = mean)
dat.pm1$ndvi.50 <- as.numeric(ndvi.50)
View(dat.pm1)


#100m buffer
ndvi.100 <- raster::extract(ndvi, st_buffer(dat.pm1, dist = 100), fun = mean)
dat.pm1$ndvi.100 <- as.numeric(ndvi.100)
View(dat.pm1)

#200m buffer
ndvi.200 <- raster::extract(ndvi, st_buffer(dat.pm1, dist = 200), fun = mean)
dat.pm1$ndvi.200 <- as.numeric(ndvi.200)
View(dat.pm1)

#500m buffer
ndvi.500 <- raster::extract(ndvi, st_buffer(dat.pm1, dist = 500), fun = mean)
dat.pm1$ndvi.500 <- as.numeric(ndvi.500)
View(dat.pm1)


#1000m buffer
ndvi.1000 <- raster::extract(ndvi, st_buffer(dat.pm1, dist = 1000), fun = mean)
dat.pm1$ndvi.1000 <- as.numeric(ndvi.1000)
View(dat.pm1)

#2000m buffer
ndvi.2000 <- raster::extract(ndvi, st_buffer(dat.pm1, dist = 2000), fun = mean)
dat.pm1$ndvi.2000 <- as.numeric(ndvi.2000)
View(dat.pm1)

#4000m buffer
beginCluster(n=2)
ndvi.4000 <- raster::extract(ndvi, st_buffer(dat.pm1, dist = 4000), fun = mean)
dat.pm1$ndvi.4000 <- as.numeric(ndvi.4000)
View(dat.pm1)
endCluster()
#extract( r , y = 1:ncell(r) , buffer = 3000 , fun = sum )

#######land use###########
landcov<- raster('data/Accra_2014_final_classification.img')
#head(landcov@data@attributes)
#temp <- crop(landcov, dat.pm1) ##check
landcov1 <- crop(landcov, extent(gama_utm)) ##crop to gama to make extract faster
plot(dat.pm1, col = 6, add = T)
plot(landcov1)
landcov1
dim(landcov1)##number of celss(x,y)
#proj4string(landcov1) = proj4string(dat.pm1) ##make projection similar

##check each land use category
length(landcov1[landcov1== 1])
length(landcov1[landcov1== 2])
length(landcov1[landcov1== 3])
length(landcov1[landcov1== 4])

#extract raster attributes related to buffers of points data
##50m buffer
luc1 <- raster::extract(landcov1, st_buffer(dat.pm1, dist = 50))
View(luc1)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.50 <- t(sapply(luc1, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.50)
head(ldc.50)
colnames(ldc.50)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
names(ldc.50)
##add to points data
dat.pm1$ld.area.50 <- ldc.50[,"low-dens"]
dat.pm1$hd.area.50 <- ldc.50[,"high-dens"]
dat.pm1$com.area.50 <- ldc.50[,"commercial"]
dat.pm1$other.area.50 <- ldc.50[,"other"]
str(dat.pm1)

#add proportion of land use category per site to points data
ldc.prop.50 <- prop.table(ldc.50, 1)
head(ldc.prop.50)
##add to points data
dat.pm1$ld.prop.50 <- ldc.prop.50[,"low-dens"]
dat.pm1$hd.prop.50 <- ldc.prop.50[,"high-dens"]
dat.pm1$com.prop.50 <- ldc.prop.50[,"commercial"]
dat.pm1$other.prop.50 <- ldc.prop.50[,"other"]

##100m buffer
luc2 <- raster::extract(landcov1, st_buffer(dat.pm1, dist = 100))
View(luc1)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.100 <- t(sapply(luc2, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.100) 
colnames(ldc.100)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
dat.pm1$ld.area.100 <- ldc.100[,"low-dens"]
dat.pm1$hd.area.100 <- ldc.100[,"high-dens"]
dat.pm1$com.area.100 <- ldc.100[,"commercial"]
dat.pm1$other.area.100 <- ldc.100[,"other"]
str(dat.pm1)
View(dat.pm1)

#add proportion of land use category per site to points data
ldc.prop.100 <- prop.table(ldc.100, 1)
head(ldc.prop.100)
##add to points data
dat.pm1$ld.prop.100 <- ldc.prop.100[,"low-dens"]
dat.pm1$hd.prop.100 <- ldc.prop.100[,"high-dens"]
dat.pm1$com.prop.100 <- ldc.prop.100[,"commercial"]
dat.pm1$other.prop.100 <- ldc.prop.100[,"other"]

##200m buffer
luc3 <- raster::extract(landcov1, st_buffer(dat.pm1, dist = 200))
View(luc1)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.200 <- t(sapply(luc3, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.200) 
colnames(ldc.200)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
dat.pm1$ld.area.200 <- ldc.200[,"low-dens"]
dat.pm1$hd.area.200 <- ldc.200[,"high-dens"]
dat.pm1$com.area.200 <- ldc.200[,"commercial"]
dat.pm1$other.area.200 <- ldc.200[,"other"]
str(dat.pm1)
View(dat.pm1)

#add proportion of land use category per site to points data
ldc.prop.200 <- prop.table(ldc.200, 1)
head(ldc.prop.200)
##add to points data
dat.pm1$ld.prop.200 <- ldc.prop.200[,"low-dens"]
dat.pm1$hd.prop.200 <- ldc.prop.200[,"high-dens"]
dat.pm1$com.prop.200 <- ldc.prop.200[,"commercial"]
dat.pm1$other.prop.200 <- ldc.prop.200[,"other"]

##500m buffer
luc4 <- raster::extract(landcov1, st_buffer(dat.pm1, dist = 500))
View(luc4)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.500 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.500) 
colnames(ldc.500)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
dat.pm1$ld.area.500 <- ldc.500[,"low-dens"]
dat.pm1$hd.area.500 <- ldc.500[,"high-dens"]
dat.pm1$com.area.500 <- ldc.500[,"commercial"]
dat.pm1$other.area.500 <- ldc.500[,"other"]
str(dat.pm1)
View(dat.pm1)

#add proportion of land use category per site to points data
ldc.prop.500 <- prop.table(ldc.500, 1)
head(ldc.prop.500)
##add to points data
dat.pm1$ld.prop.500 <- ldc.prop.500[,"low-dens"]
dat.pm1$hd.prop.500 <- ldc.prop.500[,"high-dens"]
dat.pm1$com.prop.500 <- ldc.prop.500[,"commercial"]
dat.pm1$other.prop.500 <- ldc.prop.500[,"other"]

##1000m buffer
luc5 <- raster::extract(landcov1, st_buffer(dat.pm1, dist = 1000))
View(luc5)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.1000 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.1000) 
colnames(ldc.1000)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
dat.pm1$ld.area.1000 <- ldc.1000[,"low-dens"]
dat.pm1$hd.area.1000 <- ldc.1000[,"high-dens"]
dat.pm1$com.area.1000 <- ldc.1000[,"commercial"]
dat.pm1$other.area.1000 <- ldc.1000[,"other"]
str(dat.pm1)
View(dat.pm1)

#add proportion of land use category per site to points data
ldc.prop.1000 <- prop.table(ldc.1000, 1)
head(ldc.prop.1000)
##add to points data
dat.pm1$ld.prop.1000 <- ldc.prop.1000[,"low-dens"]
dat.pm1$hd.prop.1000 <- ldc.prop.1000[,"high-dens"]
dat.pm1$com.prop.1000 <- ldc.prop.1000[,"commercial"]
dat.pm1$other.prop.1000 <- ldc.prop.1000[,"other"]

##2000m buffer
luc6 <- raster::extract(landcov1, st_buffer(dat.pm1, dist = 2000))
View(luc6)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.2000 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.2000) 
colnames(ldc.2000)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
dat.pm1$ld.area.2000 <- ldc.2000[,"low-dens"]
dat.pm1$hd.area.2000 <- ldc.2000[,"high-dens"]
dat.pm1$com.area.2000 <- ldc.2000[,"commercial"]
dat.pm1$other.area.2000 <- ldc.2000[,"other"]
str(dat.pm1)
View(dat.pm1)

#add proportion of land use category per site to points data
ldc.prop.2000 <- prop.table(ldc.2000, 1)
head(ldc.prop.2000)
##add to points data
dat.pm1$ld.prop.2000 <- ldc.prop.2000[,"low-dens"]
dat.pm1$hd.prop.2000 <- ldc.prop.2000[,"high-dens"]
dat.pm1$com.prop.2000 <- ldc.prop.2000[,"commercial"]
dat.pm1$other.prop.2000 <- ldc.prop.2000[,"other"]

##4000m buffer
luc7 <- raster::extract(landcov1, st_buffer(dat.pm1, dist = 4000))
View(luc7)
#convert list to table and multiply by area (resolution is 20 *20) ==20m*20m == 400m2
ldc.4000 <- t(sapply(luc4, FUN = function(x) table(factor(x, levels = 1:4))))* 400
View(ldc.4000) 
colnames(ldc.4000)[1:4] <- c("low-dens", "high-dens", "commercial", "other")
##add to points data
dat.pm1$ld.area.4000 <- ldc.4000[,"low-dens"]
dat.pm1$hd.area.4000 <- ldc.4000[,"high-dens"]
dat.pm1$com.area.4000 <- ldc.4000[,"commercial"]
dat.pm1$other.area.4000 <- ldc.4000[,"other"]
str(dat.pm1)
View(dat.pm1)

#add proportion of land use category per site to points data
ldc.prop.4000 <- prop.table(ldc.4000, 1)
head(ldc.prop.4000)
##add to points data
dat.pm1$ld.prop.4000 <- ldc.prop.4000[,"low-dens"]
dat.pm1$hd.prop.4000 <- ldc.prop.4000[,"high-dens"]
dat.pm1$com.prop.4000 <- ldc.prop.4000[,"commercial"]
dat.pm1$other.prop.4000 <- ldc.prop.4000[,"other"]

#########population density##############
pop.den<- raster('data/Gridded population rasters/GHS population grid -250m - 2015/pop_GAR_2015_GHS1.tif')
pop.den
plot(pop.den)
proj4string(pop.den)
#project to utm
pop.den.utm<- projectRaster(pop.den, crs = crs(pmw.dat.utm))
crs(dat.pm1)
crs(pop.den.utm)
proj4string(pop.den.utm)
hist(pop.den.utm, xlab = "Population density", main = "Distribution of Population density")
hist(log(pop.den.utm), xlab = "Population density", main = "Distribution of Population density")
plot(pop.den.utm)
plot(dat.pm1, add = T)
pop.den.utm

##crop to gama
pop.den1 <- crop(pop.den.utm, extent(gama_utm)) ##crop to gama to make extract faster
plot(pop.den1)
plot(dat.pm1, add = T)
#extract population density attributes related to buffers of points data
##50m buffer
pop.50 <- raster::extract(pop.den1, st_buffer(dat.pm1, dist = 50), fun = mean)
View(pop.50)
dat.pm1$pop.50 <- as.numeric(pop.50)


##100m buffer
pop.100 <- raster::extract(pop.den1, st_buffer(dat.pm1, dist = 100), fun = mean)
View(pop.100)
dat.pm1$pop.100 <- as.numeric(pop.100)

##200m buffer
pop.200 <- raster::extract(pop.den1, st_buffer(dat.pm1, dist = 200), fun = mean)
View(pop.200)
dat.pm1$pop.200 <- as.numeric(pop.200)

##500m buffer
pop.500 <- raster::extract(pop.den1, st_buffer(dat.pm1, dist = 500), fun = mean)
View(pop.500)
dat.pm1$pop.500 <- as.numeric(pop.500)

##1000m buffer
pop.1000 <- raster::extract(pop.den1, st_buffer(dat.pm1, dist = 1000), fun = mean)
View(pop.1000)
dat.pm1$pop.1000 <- as.numeric(pop.1000)

##2000m buffer
pop.2000 <- raster::extract(pop.den1, st_buffer(dat.pm1, dist = 2000), fun = mean)
View(pop.2000)
dat.pm1$pop.2000 <- as.numeric(pop.2000)

##4000m buffer
pop.4000 <- raster::extract(pop.den1, st_buffer(dat.pm1, dist = 4000), fun = mean)
View(pop.4000)
dat.pm1$pop.4000 <- as.numeric(pop.4000)


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
#plot(st_intersects(st_buffer(pm_dat, dist = 50), build.den.utm))
build.count.50 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 50), build.den.utm))
dat.pm1$build.count.50 <- as.numeric(build.count.50)

#intersect building density and buffer - 100
build.count.100 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 100), build.den.utm))
head(build.count.100)
dat.pm1$build.count.100 <- as.numeric(build.count.100)

#intersect building density and buffer - 200
build.count.200 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 200), build.den.utm))
dat.pm1$build.count.200 <- as.numeric(build.count.200)

#intersect building density and buffer - 500
build.count.500 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 500), build.den.utm))
dat.pm1$build.count.500 <- as.numeric(build.count.500)

#intersect building density and buffer - 1000
build.count.1000 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 1000), build.den.utm))
dat.pm1$build.count.1000 <- as.numeric(build.count.1000)

#intersect building density and buffer - 2000
build.count.2000 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 2000), build.den.utm))
dat.pm1$build.count.2000 <- as.numeric(build.count.2000)

#intersect building density and buffer - 4000
build.count.4000 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 4000), build.den.utm))
dat.pm1$build.count.4000 <- as.numeric(build.count.4000)

#distance to airport
airport <- read_sf("data/gha_airports/GHA_Airports.shp")
crs(airport)
plot(airport$geometry)
airp_utm <- st_transform(airport, crs.gama)
crs(airp_utm)
airp <- st_crop(airp_utm ,pm_dat)
air_dist <- st_distance(dat.pm1, airp)
View(air_dist)

plot(pm_dat$geometry)
plot(airp$geometry, add = T, col = 6)

#add to data
dat.pm1$dist_airport <- as.numeric(air_dist)
summary(air_dist)
View(dat.pm1)

#number of bus stops within buffer
transport <- read_sf("data/OSM 2019/ghana-latest-free-2019.shp/gis_osm_transport_free_1.shp")
plot(transport$geometry)
crs(transport)
transport_utm <- st_transform(transport, crs.gama)
crs(transport_utm)
trans_gama <- st_crop(transport_utm ,st_as_sf(gama_utm))

#temp <- st_crop(transport_utm ,pm_dat)
#plot(temp$geometry)
plot(trans_gama$geometry)
plot(dat.pm1, col = 6, add = T)
plot(gama_utm, add = T)
unique(trans_gama$fclass)
bus.stp <- trans_gama %>% filter(fclass %in% "bus_stop")
class(bus.stp)
View(bus.stp)
unique(bus.stp$fclass)

#50m buffer
bus.count <- lengths(st_intersects(st_buffer(dat.pm1, dist = 50), bus.stp))
View(bus.count)
str(bus.count)
#pm_dat$bus.count.50 <- as.numeric(bus.count)
dat.pm1$bus.count.50 <- bus.count
table(dat.pm1$bus.count.50) ##373 sites had no bus stop


#100m buffer
bus.count1 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 100), bus.stp))
View(bus.count1)
dat.pm1$bus.count.100 <- bus.count1
table(dat.pm1$bus.count.100) ##284 sites had no bus stop

#200m buffer
bus.count2 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 200), bus.stp))
View(bus.count2)
dat.pm1$bus.count.200 <- bus.count2
table(dat.pm1$bus.count.200) ##227 sites had no bus stop

#500m buffer
bus.count3 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 500), bus.stp))
View(bus.count3)
dat.pm1$bus.count.500 <- bus.count3
table(dat.pm1$bus.count.500) ##123 sites had no bus stop

#1000m buffer
bus.count4 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 1000), bus.stp))
View(bus.count4)
dat.pm1$bus.count.1000 <- bus.count4
table(dat.pm1$bus.count.1000) ##104 sites had no bus stop

#2000m buffer
bus.count5 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 2000), bus.stp))
View(bus.count5)
dat.pm1$bus.count.2000 <- bus.count5
table(dat.pm1$bus.count.2000) ##6 sites had no bus stop

#4000m buffer
bus.count6 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 4000), bus.stp))
View(bus.count6)
dat.pm1$bus.count.4000 <- bus.count6
table(dat.pm1$bus.count.4000) ##all sites had bus stop

##distance to nearest bus stop#####
bus.stp$id <-1:nrow(bus.stp)
head(bus.stp$id)
site.near.busstp <- st_join(dat.pm1, bus.stp, join = st_nearest_feature)
st_crs(site.near.busstp)
plot(site.near.busstp$geometry)
dist.bus.stp <- st_distance(dat.pm1, bus.stp[site.near.busstp$id, ], by_element = T, which = "Euclidean")
View(dist.bus.stp)
dat.pm1$dist.bus.stp <- as.numeric(dist.bus.stp)

##inverse distance to nearest bus stop#####
dat.pm1$inv.dist.busstp <- 1/(dat.pm1$dist.bus.stp)

##inverse squared distance to nearest bus stop#####
dat.pm1$invsq.dist.busstp <- (1/(dat.pm1$dist.bus.stp))^2

##number of bus stations within buffer
unique(trans_gama$fclass)
bus.stat<- trans_gama %>% filter(fclass %in% "bus_station")
View(bus.stat)
plot(bus.stat$geometry)
plot(dat.pm1, col = 6, add = T)
plot(gama_utm, add = T)

#50m buffer
bus.stat.count <- lengths(st_intersects(st_buffer(dat.pm1, dist = 50), bus.stat))
View(bus.stat.count)
summary(bus.stat.count) ##no bus station - drop covariate
dat.pm1$bus.stat.50 <- bus.stat.count

#100m buffer
bus.stat.count1 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 100), bus.stat))
View(bus.stat.count1)
summary(bus.stat.count1) ##no bus station - drop covariate
dat.pm1$bus.stat.100 <- bus.stat.count1

#200m buffer
bus.stat.count2 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 200), bus.stat))
View(bus.stat.count2)
summary(bus.stat.count2)
dat.pm1$bus.stat.200 <- bus.stat.count2
table(dat.pm1$bus.stat.200) ##417 sites had no bus station

#500m buffer
bus.stat.count3 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 500), bus.stat))
View(bus.stat.count3)
summary(bus.stat.count3)
dat.pm1$bus.stat.500 <- bus.stat.count3
table(dat.pm1$bus.stat.500) ##409 sites had no bus station

#1000m buffer
bus.stat.count4 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 1000), bus.stat))
View(bus.stat.count4)
summary(bus.stat.count4)
dat.pm1$bus.stat.1000 <- bus.stat.count4
table(dat.pm1$bus.stat.1000) ##265 sites had no bus station

#2000m buffer
bus.stat.count5 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 2000), bus.stat))
View(bus.stat.count5)
summary(bus.stat.count5)
dat.pm1$bus.stat.2000 <- bus.stat.count5
table(dat.pm1$bus.stat.2000) ##192 sites had no bus station

#4000m buffer
bus.stat.count6 <- lengths(st_intersects(st_buffer(dat.pm1, dist = 4000), bus.stat))
View(bus.stat.count6)
summary(bus.stat.count6)
dat.pm1$bus.stat.4000 <- bus.stat.count6
table(dat.pm1$bus.stat.4000) ##49 sites had no bus station

##distance to nearest bus station#####
bus.stat$id <-1:nrow(bus.stat)
head(bus.stat$id)
site.near.bus.stat <- st_join(dat.pm1, bus.stat, join = st_nearest_feature)
st_crs(site.near.bus.stat)
plot(site.near.bus.stat$geometry)
dist.bus.stat <- st_distance(dat.pm1, bus.stat[site.near.bus.stat$id, ], by_element = T, which = "Euclidean")
View(dist.bus.stat)
dat.pm1$dist.bus.stat <- as.numeric(dist.bus.stat)

##inverse distance to nearest bus station#####
dat.pm1$inv.dist.busstat <- 1/(dat.pm1$dist.bus.stat)

##inverse squared distance to nearest bus stop#####
dat.pm1$invsq.dist.busstat <- (1/(dat.pm1$dist.bus.stat))^2

####distance to railway#######
unique(transport_utm$fclass)
rail <- transport_utm %>% filter(fclass == "railway_station")
st_crs(rail)
st_crs(dat.pm1)
plot(dat.pm1$geometry)
plot(rail$geometry, add = T, col = 6)


rail$id <-1:nrow(rail)
head(rail$id)
site.near.rail <- st_join(dat.pm1, rail, join = st_nearest_feature)
st_crs(site.near.rail)
plot(site.near.rail$geometry)
dist.rail <- st_distance(dat.pm1, rail[site.near.rail$id, ], by_element = T, which = "Euclidean")
View(dist.rail)
head(dist.rail)
dat.pm1$dist.rail <- as.numeric(dist.rail)

##inverse distance to nearest bus station#####
dat.pm1$inv.dist.rail <- 1/(dat.pm1$dist.rail)

##inverse squared distance to nearest bus stop#####
dat.pm1$invsq.dist.rail <- (1/(dat.pm1$dist.rail))^2


####distance to taxi station#######
taxi <- transport_utm %>% filter(fclass == "taxi")

taxi$id <-1:nrow(taxi)
head(taxi$id)
site.near.taxi <- st_join(dat.pm1, taxi, join = st_nearest_feature)
st_crs(site.near.taxi)
plot(site.near.taxi$geometry)
dist.taxi <- st_distance(dat.pm1, taxi[site.near.taxi$id, ], by_element = T, which = "Euclidean")
View(dist.taxi)
head(dist.taxi)
dat.pm1$dist.taxi <- as.numeric(dist.taxi)

##inverse distance to nearest bus station#####
dat.pm1$inv.dist.taxi <- 1/(dat.pm1$dist.taxi)

##inverse squared distance to nearest bus stop#####
dat.pm1$invsq.dist.taxi <- (1/(dat.pm1$dist.taxi))^2

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

