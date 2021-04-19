
############# Half-hourly temporal scale for RMEP, Penman-Monteity, Priestley-Taylor model############### 
rm(list=ls())
##Data Input in .csv format

Data <- read.csv("Path/HaibeiHH0305daytime.csv")

#######    Penman-Monteith model  calculation ###########
#Rn, net radiation, original unit:W/m2, needed to be converted to MJ/m2day
Rn <-Data$NETRAD*0.0864      # Unit:MJ/m2day

#T:temperature  air temperature at 2m
T <- Data$TA_F               # daily air temperature C

#VPD, Vapor pressure Deficit,  unit:hpa
esea <- Data$VPD_F/10   #Unit:Kpa

#delta, slope vapor pressure curve, unit:Kpa/C

del <- 2504*exp(17.27*T/(T+237.2))/(T+273.3)^2

#P, atmospheric pressure,Unit:Kpa
P <- Data$PA_F  #Unit:Kpa

#gammma,psychrometric constant, Unit:Kpa/C
gam <- 0.000665*P  #p, Unit:Kpa

#U wind speed at 2m, Unit:m/s
U <- Data$WS_F

#G, soil heat flux, Unit:MJ/m2day, G=0.1Rn for daytime and G=0.5Rn for nighttime
G <- 0.1*Rn

## calculation 
PETPM <- (0.408*del*(Rn-G) + gam*900/(T+273)*U*(esea))/(del+gam*(1+0.34*U))  #Unit: mm/day


###############    Priestley-Taylor model calculation   ################
#P-T model  latent heat of vaporization=2.45 MJ/kg

ETPri <- 1.26*del/(del+gam)*(Rn-G)/2.45  #Priestley-Taylor model,Unit: mm/day

##############  RMEP calculation ########################
##Install and library RMEP package
library(RMEP)

outputPET<-RMEPPET(Rn <- Data$NETRAD,RnL <- Data$NETRAD,Ts <- Data$TS_F_MDS_1,I=600,z=5,type=2) #RMEP model,Unit: mm/day


######## Result analysis and plot ############################# 
library(hydroGOF)  ### result analysis
library(ggplot2)   ### plot
gof(sim=outputPET$PETMEP[1,],obs=PETPM)
lm(outputPET$PETMEP[1,]~PETPM)  #linear regression 

gof(sim=outputPET$PETMEP[1,],obs=ETPri) 
lm(outputPET$PETMEP[1,]~ETPri)

resultHH <- data.frame(outputPET$PETMEP[1,],PETPM,ETPri)  ## data.frame of results
write.csv(resultHH,file = "Path /half hourly result haibei.csv") ## save results, this can be ignored if not save

#########  Plot for PM and RMEP ###
ggplot(data=resultHH,mapping = aes(x=PETPM,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(-5, 20),ylim=c(-5,20))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-M PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Half-hourly)")+
  annotate("text",x=0,y=17,label="R^2=0.93\nRMSE=1.04 mm/day\ny=0.97x-0.07")
########## Plot for P-T and RMEP ###
ggplot(data=resultHH,mapping = aes(x=ETPri,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(-5, 25),ylim=c(-5,25))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-T PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Half-hourly nighttime)")+
  annotate("text",x=2,y=22,label="\nR^2=0.98\nRMSE=2.22 mm/day\ny=0.70x+0.28")


#### Plot for night time at Dinghushan ###
ggplot(data=resultHH,mapping = aes(x=PETPM,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(-2, 7),ylim=c(-2,7))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-M PET (mm/day)",y="RMEP PET (mm/day)",title = "Dinghushan(Half-hourly nighttime)")+
  annotate("text",x=0,y=5.2,label="R^2=0\nRMSE=1.90 mm/day\ny=-0.08x-0.72")

ggplot(data=resultHH,mapping = aes(x=ETPri,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(-2, 6),ylim=c(-2,6))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-T PET (mm/day)",y="RMEP PET (mm/day)",title = "Dinghushan(Half-hourly nighttime)")+
  annotate("text",x=0,y=5,label="\nR^2=0.99\nRMSE=0.47 mm/day\ny=1.62x-0.01")

#### Plot for night time at Haibei ###
ggplot(data=resultHH,mapping = aes(x=PETPM,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(-1, 5),ylim=c(-1,5))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-M PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Half-hourly nighttime)")+
  annotate("text",x=0,y=4,label="R^2=0.15\nRMSE=0.97 mm/day\ny=0.61x-0.69")

ggplot(data=resultHH,mapping = aes(x=ETPri,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(-2, 5),ylim=c(-2,5))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-T PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Half-hourly nighttime)")+
  annotate("text",x=-0.5,y=4,label="\nR^2=0.91\nRMSE=0.48 mm/day\ny=1.68x-0.14")

