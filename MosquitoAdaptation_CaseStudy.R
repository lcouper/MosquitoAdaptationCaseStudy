#### Mosquito Adaptation Case study #####

# Case study context: dengue transmission by Ae. aegypti in Northern Brazil

library(ncdf4) 
library(PCICt)
library(ncdf4.helpers) 
library(measurements)
library(weathermetrics)

#### Pull in climate data ####
Ta85 = nc_open("~/Downloads/tas_day_HadGEM2-ES_rcp85_r1i1p1_20201201-20251130.nc")
Ta85fut = nc_open("~/Downloads/tas_day_HadGEM2-ES_rcp85_r1i1p1_20751201-20851130.nc")

# extract latitude and longitude and time
lon <- ncvar_get(Ta85, varid = "lon")
lat <- ncvar_get(Ta85, varid = "lat")
time_index2021 = 31:390  # to encompass all of 2021 
time_index2080 = 1471:1830  # to encompass all of 2080

# create time series for the air surface temperature variable
Ta85_time <- nc.get.time.series(Ta85, v = "ta", time.dim.name = "time")
Ta85fut_time <- nc.get.time.series(Ta85fut, v = "ta", time.dim.name = "time")
Ta85_time[c(1:3, length(Ta85_time) - 2:0)]
Ta85fut_time[c(1:3, length(Ta85_time) - 2:0)]
Ta85_time_sub = Ta85_time[time_index2021]
Ta85fut_time_sub = Ta85fut_time[time_index2080]
Ta85_array = ncvar_get(Ta85, "tas")
Ta85fut_array = ncvar_get(Ta85fut, "tas")

x2 = kelvin.to.celsius(Ta85_array[, , time_index2021])
y2 = kelvin.to.celsius(Ta85fut_array[, ,time_index2080])
# 1st dimension = longitude, 2nd = longitude, 3rd = day of year

#### Trait thermal limits & R0 for Ae. aegypti & DENV transmission #####

# Biting rate, a, briere function. 
# units = 1/days
# For aedes aegypti, c = 0.000202, Tmin = 13.35, Tmax = 40.08
BitingRate = function(Temp)
{if (Temp > 13.35 & Temp < 40.08) 
{a = 0.000202 * Temp * (Temp - 13.35) * sqrt((40.08 - Temp))}
  else {a= 0} # if temp entered is below Tmin or above Tmax, return 0
  return(a)} 

# Probability an infected mosquito becomes infectious, b, Briere function
# For aedes aegytpi & DENV, c = 0.000849, Tmin = 17.05, Tmax = 35.83
ProbInfB = function(Temp)
{if (Temp > 17.05 & Temp < 35.83)
{b = 0.000849 * Temp * (Temp - 17.05) * sqrt((35.83 - Temp))}
  else {b = 0}
  return(b)}

# Probability of getting infected from infectious blood meal, c, Briere function
# For aedes aegypti & DENV, c = 0.000491, Tmin = 12.22, Tmax = 37.46
ProbInfC = function(Temp)
{if (Temp > 12.22 & Temp < 37.46)
{c = 0.000491 * Temp * (Temp - 12.22) * sqrt((37.46 - Temp))}
  else {c = 0}
  return(c)}

# Adult mosquito Life span (in days), 1/mu, Quadratic
# For aedes aegypti, c = -0.148, Tmin = 9.16, Tmax = 37.73
Lifespan = function(Temp)
{if (Temp > 9.16 & Temp < 37.73)
{lf = (-0.148) * (Temp - 9.16) * (Temp - 37.73)}
  else {lf = 0}
  return(lf)}

# Parasite development rate (PDR), 1/days, Briere
# For DENV, c = 0.0000665, Tmin = 10.68, Tmax = 45.90
PDR = function(Temp)
{if (Temp > 10.68 & Temp < 45.90)
{pdr = 0.0000665 * Temp * (Temp - 10.68) * sqrt((45.90 - Temp))}
  else {pdr = 0}
  return(pdr)}

# Eggs laid per female per gonotrophic cycle, TFD, #/female, Briere
# For Aedes aegypti, c = 0.00856, Tmin = 14.58, Tmax = 34.61
TFD = function(Temp)
{if (Temp > 14.58 & Temp < 34.61)
{tfd = 0.00856 * Temp * (Temp - 14.58) * sqrt((34.61 - Temp))}
  else {tfd = 0}
  return(tfd)}

# probability egg-to-adult survival, pEA, quadratic
# for Aedes aegypti, c = -0.00599, Tmin = 13.56, Tmax = 38.29
pEA = function(Temp)
{if (Temp > 13.56 & Temp < 38.29)
{pea = (-0.00599) * (Temp - 13.56) * (Temp - 38.29)}
  else {pea = 0}
  return(pea)}

# mosquito egg-to-adult development rate (1/days), MDR, Briere
# for Aedes aegypti, c = 0.0000786, Tmin = 11.36, Tmax = 39.17
MDR = function(Temp)
{if (Temp > 11.36 & Temp < 39.17)
{mdr = 0.0000786 * Temp * (Temp - 11.36) * sqrt((39.17 - Temp))}
  else {mdr = 0}
  return(mdr)}

# R0 model
R0 = function(Temp)
  {if (Temp < 17.05)
  {R0 = NA
return(R0)}
else
  {R0 = sqrt((BitingRate(Temp)^2 * ProbInfB(Temp) * ProbInfC(Temp) * exp(1)^-(1/(Lifespan(Temp) * PDR(Temp))) * 
                (TFD(Temp)*BitingRate(Temp)) * pEA(Temp) * MDR(Temp))/ (1/Lifespan(Temp))^3)
return(R0)}}

#### Calculate R0 for each month, lat & long #####

# Calculate mean monthly temp from these daily temps for 2021 and 2080
# Jan = 1:30, Feb = 31:60, March = 61:90, April = 91:121, ... Dec = 331:360

MonthlyMeanTemp2 = array(data = NA, c(192,145,12))
for (h in 1:12)
  for (i in 1:192)
    for (j in 1:145)
    {MonthlyMeanTemp2[i,j,h] = mean(x2[i,j,(1+((h-1)*30)):(h*30)])}

MonthlyMeanTempFuture2 = array(data = NA, c(192,145,12))
for (h in 1:12)
  for (i in 1:192)
    for (j in 1:145)
    {MonthlyMeanTempFuture2[i,j,h] = mean(y2[i,j,(1+((h-1)*30)):(h*30)])}

MonthlyR02 = array(data = NA, c(192,145,12))
for (h in 1:12)
  for (i in 1:192)
    for (j in 1:145)
    {MonthlyR02[i,j,h] = R0(MonthlyMeanTemp2[i,j,h])}

MonthlyR0future2 = array(data = NA, c(192,145,12))
for (h in 1:12)
  for (i in 1:192)
    for (j in 1:145)
    {MonthlyR0future2[i,j,h] = R0(MonthlyMeanTempFuture2[i,j,h])}

#### Calculate # of months with transmission in 2021 and 2080 #####

TransmissionMonths2 = matrix(nrow = 192,ncol = 145)
for (i in 1:192)
  for (j in 1:145)
  {month = as.numeric(na.omit(MonthlyR02[i,j,]))
  TransmissionMonths2[i,j] = sum(month > 0)}
rownames(TransmissionMonths2) = lon - 360
colnames(TransmissionMonths2) = lat

TransmissionMonthsFuture2 = matrix(nrow = 192,ncol = 145)
for (i in 1:192)
  for (j in 1:145)
  {month = as.numeric(na.omit(MonthlyR0future2[i,j,]))
  TransmissionMonthsFuture2[i,j] = sum(month > 0)}
rownames(TransmissionMonthsFuture2) = lon - 360
colnames(TransmissionMonthsFuture2) = lat

# How much monthly temps exceed 34.61 in the future
MonthlyOver = array(data = NA, c(192,145,12))
for (h in 1:12)
  for (i in 1:192)
    for (j in 1:145)
    {MonthlyOver[i,j,h] = (MonthlyMeanTempFuture2[i,j,h]-34.61)}

# only average the times OVER 34.61
AverageOver = matrix(nrow = 192, ncol = 145)
  for (i in 1:192)
    for (j in 1:145)
    {temp = MonthlyOver[i,j,]
      AverageOver[i,j] = mean(temp[temp > 0])}
rownames(AverageOver) = lon - 360
colnames(AverageOver) = lat


#### Create rasters and plot ####

# used examples from: https://rspatial.org/raster/rosu/Chapter11.html

library(rspatial)
library(raster)
library(sf)
library(rmapshaper)
library(reshape2)
library(RColorBrewer)
library(wesanderson)
library(bnspatial)

# Convert data frames wide to long with columns for lat, long, and value
# current transmission suitability
TransMonth2melt = melt(TransmissionMonths2)
colnames(TransMonth2melt) = c("lon", "lat", "value")
TransMonth2melt = as.matrix(TransMonth2melt)
# create raster from data frame using proper crs
TransMonth2Raster = rasterFromXYZ(TransMonth2melt, crs = "+init=epsg:4326", digits = .4)

# future transmission suitability
TransMonth2melt2 = melt(TransmissionMonthsFuture2)
colnames(TransMonth2melt2) = c("lon", "lat", "value")
TransMonth2melt2 = as.matrix(TransMonth2melt2)
TransMonth2Raster2 = rasterFromXYZ(TransMonth2melt2, crs = "+init=epsg:4326", digits = .4)

# evol. change required to maintain current transmission suitability
AverageOverMelt = melt(AverageOver)
colnames(AverageOverMelt) = c("lon", "lat", "value")
AverageOverMelt = as.matrix(AverageOverMelt)
AverageOverRaster = rasterFromXYZ(AverageOverMelt, crs = "+init=epsg:4326", digits = .4)

# Pull in Brazil shapefile & 'simplify it' for easier/quicker plotting
brazil = getData('GADM', country = 'Brazil', level = 1, type = "sp")
brazil2 = rmapshaper::ms_simplify(brazil)
brazil3 = spTransform(brazil2, try2CRS)
# list names of states
brazil3@data[["NAME_1"]]
# keep only those in the North and Northeast macroregions
brazil4 = brazil3[c(1,4,27,14,3,23,22,10,20,15,26,5,2,17,6,18),]

# Crop raster to Brazil
CurrentSuitabilityCrop <-crop(TransMonth2Raster, brazil4, snap = 'near')
CurrentSuitabilityMask <- mask(CurrentSuitabilityCrop, brazil4)
FutureSuitabilityCrop <-crop(TransMonth2Raster2, brazil4, snap = 'near')
FutureSuitabilityMask <- mask(FutureSuitabilityCrop, brazil4)
EvolChangeCrop <-crop(AverageOverRaster, brazil4, snap = 'near')
EvolChangeMask <- mask(EvolChangeCrop, brazil4)

# Plot
mypalette<-rev(brewer.pal(9,"Spectral"))
mypal = mypalette[2:9]
mypalette= (brewer.pal(9,"Blues"))[c(3,5,7,8)]
cuts = c(4,5,6,7,8,9,10,11,12)
cuts2 = c(0,1,2,3,4)

plot(CurrentSuitabilityMask, main = "Current Transmission Suitability", breaks = cuts, col = mypal, 
     cex.axis = 1.3)
plot(brazil4, fill = "white", color = "black", add=T) 

plot(FutureSuitabilityMask, main = "Future Transmission Suitability" ,breaks = cuts, col = mypal, 
     cex.axis = 1.3)
plot(brazil4, fill = "white", color = "black", add=T)

plot(EvolChangeMask, main = "Degree of Evolutionary Change Required", breaks = cuts2, col = mypalette,
     cex.axis = 1.3)
plot(brazil4, fill = "white", color = "black", add=T)

#### Calculate amount of evoltionary changed required within area of focus ####
ChangeRequired = extractByMask(AverageOverRaster, EvolChangeMask)
mean(ChangeRequired) # 1.573

TransMonthsCurrent = extractByMask(TransMonth2Raster, CurrentSuitabilityMask)
mean(TransMonthsCurrent)

TransMonthsFuture = extractByMask(TransMonth2Raster2, FutureSuitabilityMask)
mean(TransMonthsFuture)





