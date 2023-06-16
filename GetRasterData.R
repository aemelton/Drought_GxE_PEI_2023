### Script by Walker Morales for Sven Buerki's Spring 2021 VIP course
# Modifications by AEM

# Load required libraries
#install.packages("raster", "rgdal")
library(raster)
library(rgdal)
 
project.folder <- "FILEPATH"
setwd(project.folder)
#read in sample data
Samples_Data <- read.csv("FilePathForLongLatData")[,12:13] # only the 12th and 13th column had longitude and latitude data
#colnames(Samples_Data) <- c("long","lat")
long <- Samples_Data$Longitude
lat <- Samples_Data$Latitude
Samples_Data <- data.frame(long, lat)

# Read in climate data
setwd("Data/EnvironmentalData/")
env.files <- list.files(pattern = ".tif", full.names = TRUE)
envStack <- stack(env.files)
envStack <- setMinMax(envStack)
plot(envStack[[1]])

# Read in separately to trim; skip if already trimmed
AI <- raster("30s_Env_Data/AI_and_PET/ai_et0/ai_et0.tif")
PET <- raster("30s_Env_Data/AI_and_PET/et0_yr/et0_yr.tif")

# Extract climate data at collection points
setwd(project.folder)
df <- raster::extract(x = envStack, y = Samples_Data, cellnumbers=TRUE)
df.AI <- raster::extract(x = AI, y = Samples_Data, cellnumbers=TRUE)
df.PET <- raster::extract(x = PET, y = Samples_Data, cellnumbers=TRUE)
df <- cbind(df, df.AI[,2])
df <- cbind(df, df.PET[,2])
head(df)
colnames(df) <- c(colnames(df[,1:21]), "AI", "PET")

# Get long and lat and bind to clim data
longlat = cbind(df, Samples_Data)

longlat <-  longlat[!duplicated(longlat[,1]),] 
long <- longlat$long
lat <- longlat$lat
longlat <- na.omit(data.frame(long,lat))
write.csv(x = longlat, file = "Outputs/unique_longlat.csv", row.names = F)

df <-  df[!duplicated(df[,1]),] 
   
### Removes the cell numbers
df <- df[,-1]
df <- na.omit(df)
write.csv(x = df, file = "Outputs/env.csv", row.names = F)
