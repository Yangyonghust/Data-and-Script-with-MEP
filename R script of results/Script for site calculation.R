
##read and process the data
#read the EC data of different sites in csv file 
data<-read.csv("EC site.csv")  #select the file of EC site data
data <- data.frame(data)
data <- na.omit(data)  #delete strange value
data[data==-9999] <- NA   #delete strange value
data <- na.omit(data)  #delete NA value
dim(data)

# installing RMEP package
library(devtools)
devtools::install_github("Yangyonghust/RMEP")

# Library the RMEP package
library(RMEP)

##Prepare for the inputs variables, variable's name is changeable with different data files 
# Calculate specific humidity: qs, in kg/kg
Shum(TA <- data$TA_1_1_1,PA <-data$PA*1000,RH <-data$RH )

# Net radiation: Rn, in W/m2
Rn<-data$NETRAD
RnL<-Rn  #Results is not affected by RnL for type 1 and 2

# Surface temperature: Ts, in Celsius
Ts <- data$PA*1000

## RMEP calculation
output<-RMEP(Rn,RnL,qs,Ts,type=1) #type is changeable according to site's situation 

# Matrix of results
simH<-t(output$HMEP);
simE<-t(output$EMEP);
simQ<-t(output$GMEP);

## Results analysis
library(hydroGOF) # Package to fit model
library(ggplot2)  # Package for plot

gof(sim=simH[,1],obs=data$H) # Analysis for model performance
lm(simH[,1]~data$H) #Linear fit
# Time series plot of H, E and G, annotate is defined according to linear fit result
qplot(data$H,simH[,1])+geom_point(col="blue")+geom_smooth(method="lm",formula = y~x,col="red",lty=1)+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="Hobs",y="HMEP",title = "BR-Sa3")+
  annotate("text",x=-10,y=350,label="y=1.22x+18.10\nR2=0.77")

gof(sim=simE[,1],obs=data$LE) # Analysis for model performance
lm(simE[,1]~data$LE)
qplot(data$LE,simE[,1])+geom_point(col="blue")+geom_smooth(method="lm",formula = y~x,col="red",lty=1)+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="Eobs",y="EMEP",title = "BR-Sa3")+
  annotate("text",x=0,y=600,label="y=0.90x-1.17\nR2=0.84")

gof(simG[,1],obs=data$G) # Analysis for model performance
lm(simG[,1]~data$G)
qplot(data$G,simG[,1])+geom_point(col="blue")+geom_smooth(method="lm",formula = y~x,col="red",lty=1)+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="Gobs",y="GMEP",title = "BR-Sa3")+
  annotate("text",x=-30,y=300,label="y=4.07x+24.11\nR2=0.44")

# Analysis for energy balance
HE <- data$H+data$LE;RnG <- data$NETRAD-data$G;lm(HE~RnG)
qplot(RnG,HE,xlab="Rn-G",ylab="H+E",main="BR-Sa3")+geom_point(col="blue")+
  geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+geom_smooth(method="lm",formula = y~x,col="red",lty=1)+
  annotate("text",x=100,y=800,label="y=0.84x+6.95\n")
qplot(data$H+data$LE+data$G,data$NETRAD,xlab="H+LE+G",ylab="Rn")+ geom_smooth(method = lm)

