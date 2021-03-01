## MODIS data calculation
rm(list = ls()) 
# Install and Library package for processing nc data
library(ncdf4)
library(chron)
library(humidity)
# Input data
ncin <- nc_open("path/Ea_1982_2017_CR.nc") #Select path for MODIS data

# Information for variables
print(ncin) 
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
#summary longitude and latitude
lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
time <- ncatt_get(ncin,"time")

# Get data
Ea <- ncvar_get(ncin,"Ea")   #monthly_ET
dim(Ea)
Ea1 <- Ea[,,276]# 265:276 represent 1:12

# Load package for plot
library(fields)
zp <- seq(0,60,10)
## Plot the results
image.plot(lon,lat,Ea1,at <- zp,col=rev(brewer.pal(10,"RdBu")),main="Ea Nov(mm/month)")
cutpts <- seq(0,120,10)
levelplot(Ea1,cuts=7,at <- cutpts, pretty=T, 
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="MEP daily E (W/m2)")
