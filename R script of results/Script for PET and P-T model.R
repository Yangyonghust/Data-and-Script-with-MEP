## Priestley-Taylor model calculation
rm(list=ls()) 
setwd("")

# Input data
data<-read.csv("EC site.csv")  #select the file of EC site data
data <- data.frame(data)
data <- na.omit(data)  #delete strange value
data[data==-9999] <- NA   #delete strange value
data <- na.omit(data)  #delete NA value
dim(data)

#P-T model
T <-  data$TS_F_MDS_1 #Temperature
ea <- 6.1088*exp(17.27*T/(T+237.3))  # ea is Saturated vapor pressure, unit:KPa
slp <- 4098*ea/(T+237.3)^2  #the slope of vapor pressure curve,KPa/C
H <- 738                   #altitude,m,738m for Changbaishan£¬3200m for Haibei
P <-101.3*(293-0.0065*H)/293    #Pressure at altitude H (m),KPa
gam <- 1.61452*P/2.45  #KPa/C
Rn <- data$NETRAD; G <- 0
ETPri <- 1.26*slp/(slp+gam)*(Rn-G) #Priestley-Taylor model,G can be 0

# Calculate potential evapotranspiration(PET)
# z=5m for Changbaishan,z=10m for Haibei
library(RMEP)
outputPET<-RMEPPET(Rn <- data$NETRAD,RnL <- data$NETRAD,Ts <- data$TS_F_MDS_1,z=5,type=1)

##Analysis of results
library(hydroGOF)
library(ggplot2)
gof(sim=outputPET$PETMEP[1,],obs=ETPri) 
lm(outputPET$PETMEP[1,]~ETPri)
resultPET <- data.frame(outputPET$PETMEP[1,],ETPri)

## Output result
write.csv(resultPET,file = ".csv") # select path

