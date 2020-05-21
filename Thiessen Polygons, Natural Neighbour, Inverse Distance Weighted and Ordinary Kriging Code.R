# Packages

library(rspatial)
library(sp)
library(raster)
library(rgdal)
library(readr)
library(ggplot2)
library(dismo)
library(gstat)

# Data preparation

setwd("~/Uni/Raster data for Charlie/Arsenic Train and Test Data Shapefiles") # Adjust to suit your computer


BD_arsenic_150m_n10_covariates_29Mar20 <- read_csv("BD_arsenic_150m_n10_covariates_29Mar20_outlier_removed.csv")

head(BD_arsenic_150m_n10_covariates_29Mar20)

As_DF <- data.frame(BD_arsenic_150m_n10_covariates_29Mar20)

summary(As_DF)

str(As_DF)

View(As_DF)


#Look at study data

ggplot(As_DF, aes(y = Arsenic, x = seq(1, length(Arsenic))), color="steelblue")+
  geom_point(na.rm=T)


# View as spatial Map

spBD <- SpatialPoints(As_DF[,4:3], proj4string=CRS('+proj=longlat +ellps=WGS84') )

spBD <- SpatialPointsDataFrame(dspBD, As_DF)

BD_basemap <- readOGR("BD_country_polygon.shp")

BD_basemap1 <- readOGR("bgd_admbnda_adm1_bbs_20180410.shp")


plot(BD_basemap1)


# define groups for mapping

cutsAs <- c(0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,1200) 


#Plot Map

spplot(spBD, 'Arsenic', cuts = cutsAs, col.regions=rev(get_col_regions()), sp.layout = BD_basemap, pch=20, cex=0.1)



# Transform CRS to planar for analysis 


BD_CRS <- CRS("+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs")


tspBD <- spTransform(spBD, BD_CRS)


tBD <- spTransform(BD_basemap1, BD_CRS) 


# Null Model RMSE

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

null <- RMSE(mean(spBD$Arsenic), spBD$Arsenic)

null


#1# Theisen Polygons


v <- voronoi(tspBD)


plot(v)


taBD <- aggregate(tBD)


vBD <- intersect(v, taBD)


## Plot


spplot(vBD, 'Arsenic', col.regions=rev(get_col_regions()),at = seq(0, 800, 10))


### Rasterise


r <- raster(taBD, res=2500) # res in unit of CRS  


vr <- rasterize(vBD, r, 'Arsenic')


plot(vr)


#### Five fold cross validation

set.seed(5132015)


kf <- kfold(nrow(tspBD))


rmseTP <- rep(NA, 5)


for (k in 1:5) {
  test <- tspBD[kf == k, ]
  train <- tspBD[kf != k, ]
  v <- voronoi(train)
  p <- extract(v, test)
  rmseTP[k] <- RMSE(test$Arsenic, p$Arsenic)
}

rmseTP

mean(rmseTP)

1 - (mean(rmseTP) / null)


#2# Natural Neighbour 


gsNN <- gstat(formula=Arsenic~1, locations=tspBD, nmax=12, set=list(idp = 0))


nn <- interpolate(r, gsNN)


nnmsk <- mask(nn, vr) 


## Plot


spplot(nnmsk, col.regions=rev(get_col_regions()),at = seq(0, 800, 10))



### 5 fold cross validation


rmseNN <- rep(NA, 5)

for (k in 1:5) {
  test <- tspBD[kf == k, ]
  train <- tspBD[kf != k, ]
  gscv <- gstat(formula=Arsenic~1, locations=train, nmax=10, set=list(idp = 0))
  p <- predict(gscv, test)$var1.pred
  rmseNN[k] <- RMSE(test$Arsenic, p)
}

rmseNN

mean(rmseNN)

1 - (mean(rmseNN) / null)



#3# Inverse distane weighting 


gsIDW <- gstat(formula=Arsenic~1, locations=tspBD, set=list(idp = 2.3))


idw <- interpolate(r, gsIDW)


## Plot

idwr <- mask(idw, vr)


spplot(idwr, col.regions=rev(get_col_regions()), at = seq(0, 800, 10))


### 5 fold Cross Validation


rmseIDW <- rep(NA, 5)


for (k in 1:5) {
  test <- tspBD[kf == k, ]
  train <- tspBD[kf != k, ]
  gs <- gstat(formula=Arsenic~1, locations=train, set=list(idp = 2.3))
  p <- predict(gs, test)
  rmseIDW[k] <- RMSE(test$Arsenic, p$var1.pred)
}


rmseIDW

mean(rmseIDW)

1 - (mean(rmseIDW) / null)


#4# Kriging


## Spatial Points Data Frame


x <- As_DF


coordinates(x) <- ~Lon + Lat


proj4string(x) <- CRS('+proj=longlat +ellps=WGS84')


BD_CRS <- CRS("+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs")


okBD <- spTransform(x, BD_CRS)


okBD <- okBD[-zerodist(okBD)[,1],] # remove dublicate locations 


### create template raster to interpolate


bd <- spTransform(BD_basemap1, BD_CRS)


r <- raster(bd)


res(r) <- 2500  # in CRS's units 


grid <- as(r, 'SpatialGrid')


#### Variogram


gs <- gstat(formula=Arsenic~1, locations=okBD)


v <- variogram(gs, width=250)


head(v)


plot(v)


##### exponetial


fve <- fit.variogram(v, vgm(8500,"Exp",7500, 2000))


fve


###### Semivariogram Plot


plot(v, fve, pch=16, alpha=0.18, col = "red")


####### Ordinary Kriging

ok <- gstat(formula=Arsenic~1, locations=okBD, model=fve)


######## Predicted values

okp <- predict(ok, grid)


######### Plot


spplot(okp)


bok <- brick(okp)


kpmasked <- mask(bok, bd)


spplot(kpmasked$var1.pred, at = seq(-8, 810, 10), col.regions=rev(get_col_regions()))


spplot(kpmasked$var1.var, col.regions=rev(get_col_regions()))


########## Evaluation

nfolds <- 5

k <- kfold(okBD, nfolds)

krigrmse <- rep(NA, 5)

for (i in 1:nfolds) {
  
  test <- okBD[k!=i,]
  train <- okBD[k==i,]
  m <- gstat(formula=Arsenic~1, locations=train, model=fve)
  p2 <- predict(m, newdata=test, debug.level=0)$var1.pred
  krigrmse[i] <-  RMSE(test$Arsenic, p2)
  
}

krigrmse

mean(krigrmse)

1 - (mean(krigrmse) / null)


## Evaluation of all 

#Thiessen Polygon

rmseTP

mean(rmseTP)

1 - (mean(rmseTP) / null)

#Natural Neighbour

rmseNN

mean(rmseNN)

1 - (mean(rmseNN) / null)

# IDW

rmseIDW

mean(rmseIDW)

1 - (mean(rmseIDW) / null)

#Krige

krigrmse

mean(krigrmse)

1 - (mean(krigrmse) / null)  


