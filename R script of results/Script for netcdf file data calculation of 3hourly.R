
## RMEP calculations for netcdf(nc)data
# Load packages for nc read
rm(list = ls()) 
library(ncdf4)
library(chron)
library(humidity)

# Read data
ncin <- nc_open("D:/GLDAS-MEP/GLDAS/2014/GLDAS_CLSM025_D.A20140101.020.nc4")
# get variable's information
print(ncin)
# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
summary(lon)
lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
summary(lat)

#Input variable
RnL_array <- ncvar_get(ncin,"Lwnet_tavg")
Rns_array <- ncvar_get(ncin,"Swnet_tavg")
Ts1_array <- ncvar_get(ncin,"AvgSurfT_tavg") #K
Ts_array <- K2C(Ts1_array)  # K to deg
qs_array <- ncvar_get(ncin,"Qair_f_tavg")  # kg/kg   Qair_f_tavg
E_array <- ncvar_get(ncin,"Qle_tavg")   #upward, latent heat,W m-2
H_array <- ncvar_get(ncin,"Qh_tavg")    #upward,sensible heat, W m-2
G_array <- ncvar_get(ncin,"Qg_tavg")    #downward_heat_flux_in_soil, w m-2, positive downward
ET_array <- ncvar_get(ncin,"Evap_tavg")   #kg m-2 s-1

# Matrix converted to vector
RnL <- as.vector(RnL_array);Rns <- as.vector(Rns_array);Ts <- as.vector(Ts_array);qs <- as.vector(qs_array)
Rn <- RnL+Rns
E_vec <- as.vector(E_array);H_vec <- as.vector(H_array);G_vec<- as.vector(G_array)
ET_vec <- as.vector(ET_array)

# Omit NA value
qs <- na.omit(as.vector(qs))
Ts <- na.omit(as.vector(Ts))
RnL <- na.omit(as.vector(RnL))
Rns <- na.omit(as.vector(Rns))
Rn <- RnL+Rns
E <- na.omit(as.vector(E_vec))
H <- na.omit(as.vector(H_vec))
G<- na.omit(as.vector(G_vec))
ET <- na.omit(as.vector(ET_vec))

## MEP calculation
# Load RMEP
library(RMEP)
output<-RMEP(Rn=Rn,RnL=RnL,qs=qs,Ts=Ts,type=1)
simH <-t(output$HMEP); 
simE <-t(output$EMEP);
simQ <- t(output$QMEP);
simET <- t(output$ETMEP) 

## Results analysis
library(hydroGOF)
library(ggplot2)

gof(sim=simH[,1],obs=H)   #H
gof(sim=simE[,1],obs=E)   #E
gof(sim=simQ[,1],obs=G)   #Q
# Energy balance analysis
RN <- H+E+G
lm(RN~Rn)
plot(Rn,H+E+G,xlim = c(-50,250),ylim = c(-50,250))
obsET <- ET*86400   #transfer to the same unit mm/day
gof(sim=simET[,1],obs=obsET)   #ET

## Plots for nc results

## H
# copy lon, lat and time from initial netCDF data set
lon <- lon;lat <- lat;time <- time;tunits <- tunits;nlon <- nlon; nlat <- nlat
fillvalue <- NA
ETsim_array <- array(fillvalue, dim=c(nlon,nlat))
length(simET)
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)
ET_vec <- as.vector(ET_array)
ET_df <- data.frame(cbind(lonlat,ET_vec))

names(ET_df) <- c("lon","lat",paste("Evap_tavg",as.character(sim), sep="_"))
ET_df <- na.omit(ET_df)
length(E_df01$Evap_tavg_sim)
simET <- as.vector(simET)
ET_df$Evap_tavg_sim <- simET    # sim ET is made
head(na.omit(ET_df), 10)

#Replace NA with simulated value
Hsim_array <- array(fillvalue, dim=c(nlon,nlat))
j <- sapply(H_df$lon, function(x) which.min(abs(lon-x)))
k <- sapply(H_df$lat, function(x) which.min(abs(lat-x)))
fillvalue <- NA
ptm <- proc.time()
nobs <- dim(H_df)[1]
Hsim_array[cbind(j,k)] <- as.matrix(H_df[1:nobs,3])   # the final array with simulated value
proc.time() - ptm
Hsim_array
dim(H_array)
#some plots to check creation of arrays
library(fields)
library(lattice)
library(RColorBrewer)
image.plot(lon,lat,Hsim_array,col=rev(brewer.pal(10,"RdBu")),main="3H MEP H (W/m2)")
image.plot(lon,lat,H_array,col=rev(brewer.pal(10,"RdBu")),main="3H GLDAS H (W/m2)")
summary(H)
summary(simH)
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- c(-200,-100,0,100,200,300,400,500)
summary(H)
levelplot(Hsim_array ~ lon * lat, data=grid,  at=cutpts, cuts=10,pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H MEP H (W/m2)")
levelplot(H_array ~ lon * lat, data=grid, at=cutpts, cuts=10,pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H GLDAS H (W/m2)")

## E
lon <- lon;lat <- lat;time <- time;nlon <- nlon; nlat <- nlat
fillvalue <- NA
Esim_array <- array(fillvalue, dim=c(nlon,nlat))
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)
E_vec <- as.vector(E_array)
E_df0 <- data.frame(cbind(lonlat,E_vec))
names(E_df0) <- c("lon","lat",paste("Qle_tavg",as.character(1), sep="_"))
E_df <- na.omit(E_df0)
length(E_df$Qle_tavg_1)
simE <- as.vector(simE)
E_df$Qle_tavg_1 <- simE    # sim ET is made
head(na.omit(E_df), 10)

#Replace NA with simulated value
Esim_array <- array(fillvalue, dim=c(nlon,nlat))
j <- sapply(E_df$lon, function(x) which.min(abs(lon-x)))
k <- sapply(E_df$lat, function(x) which.min(abs(lat-x)))
fillvalue <- NA
ptm <- proc.time()
nobs <- dim(E_df)[1]
Esim_array[cbind(j,k)] <- as.matrix(E_df[1:nobs,3])   # the final array with simulated value
proc.time() - ptm
Esim_array
dim(E_array)
#some plots to check creation of arrays
library(fields)
library(lattice)
library(RColorBrewer)
image.plot(lon,lat,Esim_array,col=rev(brewer.pal(10,"RdBu")),main="3H MEP E (W/m2)")
image.plot(lon,lat,E_array,col=rev(brewer.pal(10,"RdBu")),main="3H GLDAS E (W/m2)")
summary(E)
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- c(-100,0,100,200,300,400,500,600)
levelplot(Esim_array ~ lon * lat, data=grid,  at=cutpts, cuts=10,pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H MEP E (W/m2)")
levelplot(E_array ~ lon * lat, data=grid, at=cutpts, cuts=10,pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H GLDAS E (W/m2)")

## G
lon <- lon;lat <- lat;time <- time;nlon <- nlon; nlat <- nlat
fillvalue <- NA
Gsim_array <- array(fillvalue, dim=c(nlon,nlat))
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)
G_vec <- as.vector(G_array)
G_df0 <- data.frame(cbind(lonlat,G_vec))
names(G_df0) <- c("lon","lat",paste("Qg_tavg",as.character(1), sep="_"))
G_df <- na.omit(G_df0)
length(G_df$Qg_tavg_1)
simG <- as.vector(simQ)
G_df$Qg_tavg_1 <- simG    # sim ET is made
head(na.omit(G_df), 10)

#Replace NA with simulated value
Gsim_array <- array(fillvalue, dim=c(nlon,nlat))
j <- sapply(G_df$lon, function(x) which.min(abs(lon-x)))
k <- sapply(G_df$lat, function(x) which.min(abs(lat-x)))
fillvalue <- NA
ptm <- proc.time()
nobs <- dim(G_df)[1]
Gsim_array[cbind(j,k)] <- as.matrix(G_df[1:nobs,3])   # the final array with simulated value
proc.time() - ptm
Gsim_array
dim(G_array)
#some plots to check creation of arrays
library(fields)
library(lattice)
library(RColorBrewer)
image.plot(lon,lat,Gsim_array,col=rev(brewer.pal(10,"RdBu")),main="3H MEP G (W/m2)")
image.plot(lon,lat,G_array,col=rev(brewer.pal(10,"RdBu")),main="3H GLDAS G (W/m2)")
summary(G)
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- c(-100,-50,0,50,100,150,200,250,300)
levelplot(Gsim_array ~ lon * lat, data=grid,  at=cutpts, cuts=10,pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H MEP G (W/m2)")
levelplot(G_array ~ lon * lat, data=grid, at=cutpts, cuts=10,pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H GLDAS G (W/m2)")

## ET
lon <- lon;lat <- lat;nlon <- dim(lon); nlat <- dim(lat)
fillvalue <- NA
ETsim_array <- array(fillvalue, dim=c(nlon,nlat))
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)
ET_vec <- as.vector(ET_array)
ET_df0 <- data.frame(cbind(lonlat,ET_vec))
names(ET_df0) <- c("lon","lat",paste("Evap_tavg",as.character(1), sep="_"))
ET_df <- na.omit(ET_df0)
length(ET_df$Evap_tavg_1)
simET <- as.vector(simET)
ET_df$Evap_tavg_1 <- simET    # sim ET is made
head(na.omit(ET_df), 10)

#Replace NA with simulated value
ETsim_array <- array(fillvalue, dim=c(nlon,nlat))
j <- sapply(ET_df$lon, function(x) which.min(abs(lon-x)))
k <- sapply(ET_df$lat, function(x) which.min(abs(lat-x)))
fillvalue <- NA
ptm <- proc.time()
nobs <- dim(ET_df)[1]
ETsim_array[cbind(j,k)] <- as.matrix(ET_df[1:nobs,3])   # the final array with simulated value
proc.time() - ptm
ETsim_array
dim(ET_array)
#some plots to check creation of arrays
library(fields)
library(lattice)
library(RColorBrewer)
image.plot(lon,lat,ETsim_array,col=rev(brewer.pal(10,"RdBu")),main="3H MEP ET (mm/day)")
image.plot(lon,lat,ET_array*86400,col=rev(brewer.pal(10,"RdBu")),main="3H GLDAS ET (mm/day)")
summary(na.omit(as.vector(simET)))
summary(ET*86400)
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- c(-2,0,2,4,6,8,10,12,14,16,18)
levelplot(ETsim_array ~ lon * lat, data=grid, at=cutpts, cuts=10, pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H MEP ET (mm/day)")
levelplot(ET_array*86400 ~ lon * lat, at=cutpts, cuts=10,data=grid, pretty=T,
          col.regions=(rev(brewer.pal(10,"RdBu"))),main="3H GLDAS ET (mm/day)")

