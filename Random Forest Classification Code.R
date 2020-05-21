
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

covariate_brick <- brick("BD_hydro_variables_raster.tif")

#Rename the band names based on their spectra:

names(covariate_brick) <- c("DEM","RAIN","GWD","SAND","LOAM","HCOND","USC","RCHG","RTRND","IRRIG")


# Import the training shapefile:

As_train <- readOGR(dsn=".",layer="As_Train_data")


# Repeat the process for the testing dataset:

As_test <-  readOGR(dsn=".",layer="As_Test_data")


# Check out files

covariate_brick

As_train 

As_test


# View Data class 

class(covariate_brick)

class(As_train)

class(As_test)


#Plot Data

plot(covariate_brick)

plot(As_train)


# Extract the multispectral data from under each of the training and testing datasets, and merge it the extracted spectra with the point data.

As_training_points <- extract(covariate_brick, As_train, df=TRUE)

As_training_points<- spCbind(As_train, As_training_points)



# Look at the new attribute table:

as.data.frame(As_training_points)



# Do the same thing with the testing data:

As_testing_points <- extract(covariate_brick, As_test,df=TRUE)

As_testing_points1<- spCbind(As_test, As_testing_points)


# Create a randomForest model:


AsRandomForest <- randomForest(BinaryAs ~ DEM + RAIN + GWD + SAND + LOAM + HCOND + USC + RCHG + RTRND + IRRIG,
                                data=As_training_points,
                                importance=TRUE,
                                ntree = 500,
                                mtry=3) 

AsRandomForest 

importance(AsRandomForest )

varImpPlot(AsRandomForest , main="Variable Importance")


#  Evaulate



# Now calculate the confusionMatrix and accuracy statistics: 

AsRandomForest_test <- predict(AsRandomForest, As_testing_points1)


AsRandomForest_test <- as.factor(AsRandomForest_test)


AsRandomForest_confusionMatrix <- 
  confusionMatrix(
    data=AsRandomForest_test,
    reference=As_testing_points1$BinaryAs, positive='1')


AsRandomForest_confusionMatrix


# Random forest Prediction of unclassified data Plot


As_probability <- predict(covariate_brick,
                          AsRandomForest, type="prob")

plot(As_probability)


BD1 <- readOGR("bgd_admbnda_adm1_bbs_20180410.shp")


As_probability_mask <- mask(As_probability, BD1)



#Final plot

my.col <- colorRampPalette(c("grey","yellow","red","blue","navy"), space="rgb")

spplot(As_probability_mask, col.regions=(my.col))

