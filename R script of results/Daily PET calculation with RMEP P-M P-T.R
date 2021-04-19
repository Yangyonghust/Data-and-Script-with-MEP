
############# Daily temporal scale for RMEP, Penman-Monteity, Priestley-Taylor model############### 
rm(list=ls())
##Data Input in .csv format

Data <- read.csv("Path/Dinghushan0305DD.csv")

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

#G, soil heat flux, Unit:MJ/m2day
G <- 0  #G is neglected in daily data

## calculation 
PETPM <- (0.408*del*(Rn-G) + gam*900/(T+273)*U*(esea))/(del+gam*(1+0.34*U))  #Unit: mm/day


###############    Priestley-Taylor model calculation   ################
#P-T model  latent heat of vaporization=2.45 MJ/kg

ETPri <- 1.26*del/(del+gam)*(Rn-G)/2.45  #Priestley-Taylor model, Unit:mm/day

##############  RMEP calculation ########################
##Install and library RMEP package
library(RMEP)

outputPET<-RMEPPET(Rn <- Data$NETRAD,RnL <- Data$NETRAD,Ts <- Data$TS_F_MDS_1,I=600,z=5,type=2)  #RMEP model, Unit:mm/day

######## Result analysis and plot ############################# 
library(hydroGOF)  ### result analysis
library(ggplot2)   ### plot
gof(sim=outputPET$PETMEP[1,],obs=PETPM)
lm(outputPET$PETMEP[1,]~PETPM)  #linear regression 

gof(sim=outputPET$PETMEP[1,],obs=ETPri) 
lm(outputPET$PETMEP[1,]~ETPri)

resultDD <- data.frame(outputPET$PETMEP[1,],PETPM,ETPri)  ## data.frame of results
write.csv(resultDD,file = "Path /half hourly result haibei.csv") ## save results, this can be ignored if not save

################ Plots between PM and RMEP at Dinghushan ###
gof(sim=outputPET$PETMEP[1,],obs=PETPM)
lm(outputPET$PETMEP[1,]~PETPM)
plot(PETPM,outputPET$PETMEP[1,])
ggplot(data=resultDD,mapping = aes(x=PETPM,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 8),ylim=c(0,8))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-M PET (mm/day)",y="RMEP PET (mm/day)",title = "Dinghushan(Daily)")+
  annotate("text",x=1.5,y=7,label="R^2=0.85\nRMSE=0.96 mm/day\ny=0.76x-0.08")
################ Plots between P-T and RMEP at Dinghushan ###
plot(ETPri,outputPET$PETMEP[1,])
gof(sim=outputPET$PETMEP[1,],obs=ETPri) 
lm(outputPET$PETMEP[1,]~ETPri)
ggplot(data=resultDD,mapping = aes(x=ETPri,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 8),ylim=c(0, 8))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-T PET (mm/day)",y="RMEP PET (mm/day)",title = "Dinghushan(Daily)")+
  annotate("text",x=1.5,y=7,label="\nR^2=1.00\nRMSE=0.65 mm/day\ny=0.77x+0.08")

################ Plots between P-M and RMEP at Haibei ###
gof(sim=outputPET$PETMEP[1,],obs=PETPM)
lm(outputPET$PETMEP[1,]~PETPM)
plot(PETPM,outputPET$PETMEP[1,])
ggplot(data=resultDD,mapping = aes(x=PETPM,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 6.5),ylim=c(0,6.5))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-M PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Daily)")+
  annotate("text",x=1.5,y=5.5,label="R^2=0.91\nRMSE=0.47 mm/day\ny=0.93x-0.16")

################ Plots between P-T and RMEP at Haibei ###
plot(ETPri,outputPET$PETMEP[1,])
gof(sim=outputPET$PETMEP[1,],obs=ETPri) 
lm(outputPET$PETMEP[1,]~ETPri)
ggplot(data=resultDD,mapping = aes(x=ETPri,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 6.5),ylim=c(0, 6.5))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-T PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Daily)")+
  annotate("text",x=1.5,y=5.5,label="\nR^2=0.99\nRMSE=0.60 mm/day\ny=0.73x+0.10")


