## Library packages needed for nc data read in China region grid 
rm(list = ls()) 
library(ncdf4)
library(chron)
library(humidity)

## read data
# set path and filename
ncin <- nc_open("F:/HSun_5D3_02NASA_DAS/GLDAS_Monthly/GLDAS_NOAH025_M.A200412.021.nc4")
print(ncin) #read information of variables

# get longitude and latitude
lon <- ncvar_get(ncin,"lon")[1013:1260]#China
nlon <- dim(lon)
lat <- ncvar_get(ncin,"lat")[253:452]#China
nlat <- dim(lat)

##Input variable
# Rn
RnL_array <- ncvar_get(ncin,"Lwnet_tavg")[1013:1260,253:452]
Rns_array <- ncvar_get(ncin,"Swnet_tavg")[1013:1260,253:452]
Rn_array <- RnL_array+Rns_array
#Ts
Ts1_array <- ncvar_get(ncin,"AvgSurfT_inst")[1013:1260,253:452]#K
Ts_array <- K2C(Ts1_array)  # K to deg
#qs 
qs_array <- ncvar_get(ncin,"Qair_f_inst")[1013:1260,253:452]  # kg/kg   Qair_f_tavg

#E, H, G and ET
E_array <- ncvar_get(ncin,"Qle_tavg")[1013:1260,253:452]   #upward, latent heat,W m-2
H_array <- ncvar_get(ncin,"Qh_tavg") [1013:1260,253:452]   #upward,sensible heat, W m-2
G_array <- ncvar_get(ncin,"Qg_tavg") [1013:1260,253:452]   #downward_heat_flux_in_soil, w m-2, positive downward
ET_array <- ncvar_get(ncin,"Evap_tavg")[1013:1260,253:452]  #kg m-2 s-1
dim(ET_array)

## Matrix convert to vector
Rn <- Rn_array;RnL <- RnL_array;Ts <- Ts_array;qs <- qs_array;
E <- E_array;H <- H_array;G <- G_array;ET <-ET_array 
RnL <- as.vector(RnL_array);Rns <- as.vector(Rns_array);Ts <- as.vector(Ts_array);qs <- as.vector(qs_array)
E_vec <- as.vector(E_array);H_vec <- as.vector(H_array);G_vec<- as.vector(G_array)
ET_vec <- as.vector(ET_array)
summary(ET_vec)

# Omit NA value
qs <- na.omit(as.vector(qs))
Ts <- na.omit(as.vector(Ts))
RnL <- na.omit(as.vector(RnL))
Rns <- na.omit(as.vector(Rns))
Rn <- na.omit(as.vector(Rn))
E <- na.omit(as.vector(E_vec))
H <- na.omit(as.vector(H_vec))
G<- na.omit(as.vector(G_vec))
ET <- na.omit(as.vector(ET_vec))

## RMEP calculation
# Load RMEP
library(RMEP)
output<-RMEP(Rn=Rn,RnL=RnL,qs=qs,Ts=Ts,type=1)
simH <-t(output$HMEP); 
simE <-t(output$EMEP);
simQ <- t(output$QMEP);
simET <- t(output$EMEP*10^(-6)/2.45*86400*31) #mm/momth

## Results analysis
gof(sim=simH[,1],obs=H)   #H
gof(sim=simE[,1],obs=E)   #E
gof(sim=simQ[,1],obs=G)   #Q

#copy lon, lat and time from initial netCDF data set
#China lon lat
lon <- lon[lon>73&lon<135];lat <- lat[lat>3&lat<53];
lon <- lon;lat <- lat;time <- time;nlon <- dim(lon); nlat <- dim(lat)
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
dim(ETsim_array)

## Plot for results
library(lattice)
library(RColorBrewer)
library(fields)

#determine range
image.plot(lon,lat,ETsim_array,col=rev(brewer.pal(10,"RdBu")),main="MEP ET Nov (mm/month)")
image.plot(lon,lat,ETsim_array,col=rainbow(150)[1:60],main="MEP ET Nov (mm/month)")

#Figure
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- seq(-10,100,10)
levelplot(ETsim_array ~ lon * lat, data=grid, at=cutpts,cuts=11,pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu"))),main="MEP ET Aug (mm/month)")
# Grid data for sites
lon[221];lat[158];ETsim_array[221,158];#Changbaishan
lon[73];lat[110];ETsim_array[73,110];#Dangxiong
lon[159];lat[81];ETsim_array[159,81];#Dinghushan
lon[113];lat[138];ETsim_array[113,138];#Haibei
lon[169];lat[95];ETsim_array[169,95] #Qianyanzhou

