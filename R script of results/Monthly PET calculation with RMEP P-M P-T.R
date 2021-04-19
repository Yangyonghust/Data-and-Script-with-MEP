
############# Monthly temporal scale for RMEP, Penman-Monteity, Priestley-Taylor model############### 
rm(list=ls())
##Data Input in .csv format

Data <- read.csv("Path/DinghushanMM0305.csv")

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

#G, soil heat flux, Unit:MJ/m2day, G=0.14*(Tmonth,i - Tmonth,i-1)
G <- Data$G_calculated

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

resultMM <- data.frame(outputPET$PETMEP[1,],PETPM,ETPri)  ## data.frame of results
write.csv(resultMM,file = "Path /half hourly result haibei.csv") ## save results, this can be ignored if not save

################ Plots between PM and RMEP at Dinghushan ###
gof(sim=outputPET$PETMEP[1,],obs=PETPM)
lm(outputPET$PETMEP[1,]~PETPM)
plot(PETPM,outputPET$PETMEP[1,])
ggplot(data=resultmm,mapping = aes(x=PETPM,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 4),ylim=c(0,4))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-M PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Monthly)")+
  annotate("text",x=1,y=3.5,label="R^2=0.98\nRMSE=0.32 mm/day\ny=1.10x-0.43")
################ Plots between P-T and RMEP at Dinghushan ###
plot(ETPri,outputPET$PETMEP[1,])
gof(sim=outputPET$PETMEP[1,],obs=ETPri) 
lm(outputPET$PETMEP[1,]~ETPri)
ggplot(data=resultmm,mapping = aes(x=ETPri,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 4),ylim=c(0, 4))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-T PET (mm/day)",y="RMEP PET (mm/day)",title = "Haibei(Monthly)")+
  annotate("text",x=1,y=3.5,label="\nR^2=0.99\nRMSE=0.47 mm/day\ny=0.78x+0.03")

################ Plots between P-M and RMEP at Haibei ###
gof(sim=outputPET$PETMEP[1,],obs=PETPM)
lm(outputPET$PETMEP[1,]~PETPM)
plot(PETPM,outputPET$PETMEP[1,])
ggplot(data=resultmm,mapping = aes(x=PETPM,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 5),ylim=c(0,5))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-M PET (mm/day)",y="RMEP PET (mm/day)",title = "Dinghushan(Monthly)")+
  annotate("text",x=1,y=4.5,label="R^2=0.82\nRMSE=0.87 mm/day\ny=0.80x-0.21")
################ Plots between P-T and RMEP at Haibei ###
plot(ETPri,outputPET$PETMEP[1,])
gof(sim=outputPET$PETMEP[1,],obs=ETPri) 
lm(outputPET$PETMEP[1,]~ETPri)
ggplot(data=resultmm,mapping = aes(x=ETPri,y=outputPET$PETMEP[1,]))+geom_point(col="blue")+
  geom_smooth(method="lm",formula = y~x,col="red",lty=1)+coord_cartesian(xlim = c(0, 5),ylim=c(0, 5))+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="P-T PET (mm/day)",y="RMEP PET (mm/day)",title = "Dinghushan(Monthly)")+
  annotate("text",x=1,y=4.5,label="\nR^2=0.99\nRMSE=0.56 mm/day\ny=0.78x+0.07")


