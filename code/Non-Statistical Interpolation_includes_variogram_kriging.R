###########Non-statistical Interpolation: CA Temperatures########

if (!require("rspatial")) devtools::install_github('rspatial/rspatial')
require(rspatial)
d <- sp_data('precipitation')
class(d)
head(d)
#Compute annual precipitation
d$prec <- rowSums(d[, c(6:17)])
plot(
  sort(d$prec), 
  main = "CA annual precipitation",
  ylab = 'Annual precipitation (mm)', xlab = 'Stations',
  las = 1 )

#Now make a quick and dirty (and ugly) map:
  
  require(sp)
dsp <- SpatialPoints(d[,4:3], proj4string = CRS("+proj=longlat +datum=NAD83"))
dsp <- SpatialPointsDataFrame(dsp, d)
CA <- sp_data("counties")

# define groups for mapping
cuts <- c(0,200,300,500,1000,3000)

# set up a palette of interpolated colors
blues <- colorRampPalette(c('yellow', 'orange', 'blue', 'dark blue'))
pols <- list("sp.polygons", CA, fill = "lightgray")
spplot(dsp, 'prec', cuts = cuts, col.regions = blues(5), sp.layout = pols, pch = 20, cex = 2)


#Transform longitude/latitude to planar coordinates, using the commonly used coordinate reference system for California 
#(“Teale Albers”) to assure that our interpolation results will align with other data sets we have.

TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0
          +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83
          +units=m +ellps=GRS80 +towgs84=0,0,0")
require(rgdal)
dta <- spTransform(dsp, TA) #spdf - points projected
cata <- spTransform(CA, TA) #california shp projected

#3.2.1 NULL model - We are going to interpolate (estimate for unsampled locations) the precipitation values. 
#The simplest way would be to take the mean of all observations. We can consider that a “Null-model” that we can compare other approaches to. We’ll use the Root Mean Square Error (RMSE) as evaluation statistic.

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm = TRUE))
}
#What are some advantages of using the RMSE instead of the Sum of Squared Errors?
  
#Get the RMSE for the Null-model

# this variable name could lead to ambiguity.
null.mod <- RMSE(mean(dsp$prec), dsp$prec)
null.mod

###########3.2.2 Proximity Polygon Interpolation #############
##Proximity polygons can be used to interpolate categorical variables. 
#Another term for this is nearest neighbour interpolation.
require(dismo)
v <- voronoi(dta)
v
plot(v)

### confine the polygons to the California borders
require(rgeos)
cata
ca <- aggregate(cata)
ca
vca <- raster::intersect(v, ca)
spplot(vca, 'prec', col.regions = rev(get_col_regions()))

plot(cata, main = "cata")
plot(ca, main = "ca")

#These are polygons. We can ‘rasterize’ the results like this.

r <- raster(cata, res = 10000)
vr <- rasterize(vca, r, 'prec') 
plot(vr)

#We can use a cross validation to compare our null model that the best interpolation is to set all locations 
#to the statewide mean precipitation. Cross-validation is a simulation method in which we create a model 
#using a portion of the data, i.e. the training set, and then calculate the discrepancies between the fitted 
#and observed values for the remaining portion of the data, i.e. the testing set.

set.seed(5132015)
kf <- kfold(nrow(dta))

rmse <- rep(NA, 5)
for (k in 1:5) {
  test <- dta[kf ==  k, ]
  train <- dta[kf !=  k, ]
  v <- voronoi(train)
  p <- raster::extract(v, test)
  rmse[k] <- RMSE(test$prec, p$prec)
}
rmse

mean(rmse)
## [1] 196.7708
1 - (mean(rmse) / null.mod)


#3.2.3 Nearest neighbour interpolation
#Here we do nearest neighbor interpolation considering multiple (5) neighbors.
#We can use the gstat package for this. First we fit a model. ~1 means “intercept only”. 
#In the case of spatial data, that would be only ‘x’ and ‘y’ coordinates are used. 
#We set the maximum number of points to 5, and the “inverse distance power” idp to zero, 
#such that all five neighbors are equally weighted
require(gstat)
gs <- gstat(
  formula = prec~1,
  locations = dta, 
  nmax = 5, 
  set = list(idp = 0))
nn <- interpolate(r, gs)
plot(nn)

nnmsk <- mask(nn, vr)
plot(nnmsk)

#NOTE: gstat() creates an object of class gstat.

#Cross validate the result. Note that we can use the predict method to get predictions for the 
#locations of the test points.

rmsenn <- rep(NA, 5)
for (k in 1:5) {
  test <- dta[kf ==  k, ]
  train <- dta[kf !=  k, ]
  gscv <- gstat(
    formula = prec~1, 
    locations = train, 
    nmax = 5, 
    set = list(idp = 0))
  p <- predict(gscv, test)$var1.pred
  rmsenn[k] <- RMSE(test$prec, p)
}

rmsenn
mean(rmsenn)

1 - (mean(rmsenn) / null.mod)


#3.2.4 Inverse distance weighted
#A more commonly used method is “inverse distance weighted” interpolation. The only difference with 
#the nearest neighbor approach is that points that are further away get less weight in predicting a value a location.

require(gstat)
gs <- gstat(formula = prec~1, locations = dta)
idw <- interpolate(r, gs)
idwr <- mask(idw, vr)
plot(idwr)


#Question 4: IDW generated rasters tend to have a noticeable artefact. What is that?
  
#Cross validate. We can predict to the locations of the test points

rmse <- rep(NA, 5)
for (k in 1:5) {
  test <- dta[kf ==  k, ]
  train <- dta[kf !=  k, ]
  gs <- gstat(formula = prec~1, locations = train)
  p <- predict(gs, test)
  rmse[k] <- RMSE(test$prec, p$var1.pred)
}
rmse
mean(rmse)
1 - (mean(rmse) / null.mod)

#Question 5: Inspect the arguments used for and make a map of the IDW model below. What other name could you give to this method (IDW with these parameters)? Why?
gs2 <- gstat(formula = prec~1, locations = dta, nmax = 1, set = list(idp = 1))



#3.3.1 Data preparation
#We use the airqual dataset to interpolate ozone levels for California (averages for 1980-2009). 
#Use the variable OZDLYAV (unit is parts per billion). Original data source.

#Read the data file. To get easier numbers to read, I multiply OZDLYAV with 1000

require(rspatial)
x <- sp_data("airqual")
x$OZDLYAV <- x$OZDLYAV * 1000
#Create a SpatialPointsDataFrame and transform to Teale Albers. Note the units = km, 
#which was needed to fit the variogram.

require(sp)
require(rgdal)
coordinates(x) <- ~LONGITUDE + LATITUDE
proj4string(x) <- CRS('+proj=longlat +datum=NAD83')
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +ellps=GRS80")
aq <- spTransform(x, TA)
#Create an template raster to interpolate to. E.g., given a SpatialPolygonsDataFrame of California, ‘ca’. 
#Coerce that to a ‘SpatialGrid’ object (a different representation of the same idea)

cageo <- sp_data('counties.rds')
ca <- spTransform(cageo, TA)
r <- raster(ca)
res(r) <- 10  # 10 km if your CRS's units are in km
g <- as(r, 'SpatialGrid')
#3.3.2 Fit a variogram
#Use gstat to create an emperical variogram ‘v’

require(gstat)
gs <- gstat(formula = OZDLYAV~1, locations = aq)
v <- variogram(gs, width = 20)
head(v)
plot(v)

#Now, fit a model variogram

fve <- fit.variogram(v, vgm(85, "Exp", 75, 20))
fve

plot(variogramLine(fve, 400), type = 'l', ylim = c(0,120))
points(v[,2:3], pch = 20, col = 'red')

#Try a different type (spherical in stead of exponential)

fvs <- fit.variogram(v, vgm(85, "Sph", 75, 20))
fvs

plot(variogramLine(fvs, 400), type = 'l', ylim = c(0,120) ,col = 'blue', lwd = 2)
points(v[,2:3], pch = 20, col = 'red')

#Both look pretty good in this case. Another way to plot the variogram and the model

plot(v, fve)


#3.3.3 Ordinary kriging
#Use variogram fve in a kriging interpolation

k <- gstat(formula = OZDLYAV~1, locations = aq, model = fve)
# predicted values
kp <- predict(k, g)
spplot(kp)

#variance
ok <- brick(kp)
ok <- mask(ok, ca)
names(ok) <- c('prediction', 'variance')
plot(ok)

#3.3.4 Compare with other methods
#Let’s use gstat again to do IDW interpolation. The basic approach first.

require(gstat)
idm <- gstat(formula = OZDLYAV~1, locations = aq)
idp <- interpolate(r, idm)
idp <- mask(idp, ca)
plot(idp)


#We can find good values for the idw parameters (distance decay and number of neighbors) through optimization. 
#For simplicity’s sake I do not do that k times here. The optim function may be a bit hard to grasp at first. 
#But the essence is simple. You provide a function that returns a value that you want to minimize (or maximize) 
#given a number of unknown parameters. Your provide initial values for these parameters, and optim then searches 
#for the optimal values (for which the function returns the lowest number).

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm = TRUE))
}

f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula = OZDLYAV~1, locations = train, nmax = nmx, set = list(idp = idp))
  p <- predict(m, newdata = test, debug.level = 0)$var1.pred
  RMSE(test$OZDLYAV, p)
}
set.seed(20150518)
i <- sample(nrow(aq), 0.2 * nrow(aq))
tst <- aq[i,]
trn <- aq[-i,]
opt <- optim(c(8, .5), f1, test = tst, train = trn)
opt


#Our optimal IDW model

m <- gstat(formula = OZDLYAV~1, locations = aq, nmax = opt$par[1], set = list(idp = opt$par[2]))
idw <- interpolate(r, m)

idw <- mask(idw, ca)
plot(idw)


#A thin plate spline model

#NOTE: We haven’t talked about these in class. It is just being used here for a comparison to the other interpolation methods.

require(fields)
m <- Tps(coordinates(aq), aq$OZDLYAV)
tps <- interpolate(r, m)
tps <- mask(tps, idw)
plot(tps)

#3.3.5 Cross-validate
#Cross-validate the three methods (IDW, Ordinary kriging, TPS) and add RMSE weighted ensemble model.

require(dismo)

nfolds <- 5
k <- kfold(aq, nfolds)

ensrmse <- tpsrmse <- krigrmse <- idwrmse <- rep(NA, 5)

for (i in 1:nfolds) {
  test <- aq[k != i,]
  train <- aq[k == i,]
  m <- gstat(formula = OZDLYAV~1, locations = train, nmax = opt$par[1], set = list(idp = opt$par[2]))
  p1 <- predict(m, newdata = test, debug.level = 0)$var1.pred
  idwrmse[i] <-  RMSE(test$OZDLYAV, p1)
  
  m <- gstat(formula = OZDLYAV~1, locations = train, model = fve)
  p2 <- predict(m, newdata = test, debug.level = 0)$var1.pred
  krigrmse[i] <-  RMSE(test$OZDLYAV, p2)
  
  m <- Tps(coordinates(train), train$OZDLYAV)
  p3 <- predict(m, coordinates(test))
  tpsrmse[i] <-  RMSE(test$OZDLYAV, p3)
  
  w <- c(idwrmse[i], krigrmse[i], tpsrmse[i])
  weights <- w / sum(w)
  ensemble <- p1 * weights[1] + p2 * weights[2] + p3 * weights[3]
  ensrmse[i] <-  RMSE(test$OZDLYAV, ensemble)
}
rmi <- mean(idwrmse)
rmk <- mean(krigrmse)
rmt <- mean(tpsrmse)
rms <- c(rmi, rmt, rmk)
rms

rme <- mean(ensrmse)
rme


#Question 6: Which method performed best?
#We can use the rmse scores to make a weighted ensemble. Let’s look at the maps

weights <- ( rms / sum(rms) )
s <- stack(idw, ok[[1]], tps)
ensemble <- sum(s * weights)
#And compare maps.

s <- stack(idw, ok[[1]], tps, ensemble)
names(s) <- c('IDW', 'OK', 'TPS', 'Ensemble')
plot(s)
