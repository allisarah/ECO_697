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



