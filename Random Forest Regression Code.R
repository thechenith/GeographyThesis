
#Packages
require("raster")
require("rgdal")
require("sp")
require("maptools")
require("randomForest")
require("caret")
require("e1071")
require("readr")


setwd("~/Uni/Raster data for Charlie/Arsenic Train and Test Data Shapefiles")


#Import Data

covariatebrick <- brick("BD_hydro_variables_raster.tif")

#Rename the band names based on their spectra:

names(covariatebrick) <- c("DEM","RAIN","GWD","SAND","LOAM","HCOND","USC","RCHG","RTRND","IRRIG")


# Import the training shapefile:

As_training_points <- readOGR(dsn=".",layer="As_Train_data")


# Repeat the process for the testing dataset:

As_testing_points <-  readOGR(dsn=".",layer="As_Test_data")


# Check out files

covariatebrick

As_training_points

As_testing_points


# View Data class 

class(covariatebrick)

class(As_training_points)

class(As_testing_points)


#Plot Data


BD1 <- readOGR("bgd_admbnda_adm1_bbs_20180410.shp")

covariatebrickmasked <- mask(covariatebrick, BD1)

plot(covariatebrickmasked)


# Extract the multispectral data from under each of the training and testing datasets, and merge it the extracted spectra with the point data.

As_training_points1 <- extract(covariatebrick,
                                               As_training_points,df=TRUE)

As_training_points1.1<- spCbind(As_training_points,
                                                As_training_points1)


# Look at the new attribute table:


as.data.frame(As_training_points1.1)


# Do the same thing with the testing data:


As_testing_points1 <- extract(covariatebrick,
                                              As_testing_points,df=TRUE)

As_testing_points1.1 <- spCbind(As_testing_points,
                                               As_testing_points1)


# Create a randomForest model:


As_randomForest <- randomForest(Arsenic ~ DEM + RAIN + GWD + SAND + LOAM + HCOND + USC + RCHG + RTRND + IRRIG,
                                data=As_training_points1.1,
                                importance=TRUE,
                                ntree = 500,
                                mtry=3,
                              )

As_randomForest



# Random forest Prediction of unclassified data Plot

As_RF_pred <- predict(covariatebrick,
                             As_randomForest)

spplot(As_RF_pred,  col.regions=rev(get_col_regions()),at = seq(0, 800, 10))


BD <- readOGR("BD_country_polygon.shp")


As_RF_pred_mask <- mask(As_RF_pred, BD)


#Final

spplot(As_RF_pred_mask,  col.regions=rev(get_col_regions()),at = seq(0, 800, 10))



# Evaluation

#RMSE

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

RMSE(lm_randomForest$predicted, As_training_points1.1$Arsenic)



As_randomForest_test <- randomForest(Arsenic ~ DEM + USC + RAIN + SAND + LOAM + RCHG + RTRND + HCOND + GWD + IRRIG,
                       
                       data=As_training_points1.1, ntree=500, mtry=3, importance=TRUE)



# Test the prediction using the test data table

As_randomForest_test

As_randomForest_test_pred <- predict(As_randomForest_test, As_testing_points1.1)



As_testing_points1.1$PredArsen <- as.vector(As_randomForest_test_pred)


head(As_testing_points1.1$PredArsen)


summary(As_testing_points1.1)



# Run Linear Regression

summary(lm(As_testing_points1.1$Arsenic ~ As_testing_points1.1$PredArsen))
