##############################################################################################
########### Teporales Data Mining, SS 2016 ####################################################
########### Projekt 2, Wissen aus Hydrologie ##################################################
###############################################################################################

###############################################################################################
#Aufgabe 1
#Inspizieren sie die einzelnen Zeitreihen mit den Methoden der Knowledge Discovery
############################################################################################
FileName <-c('Chirimachay_NO3','Chirimachay_Netradiation.lrn','Chirimachay_Precip.lrn')

## Pfad bitte hier anpassen!!!
Directory <-'C:/#Private/Studium/SoSe 2016/Temporales Data Mining/Uebung/Zettel2'


NO3=ReadLRN(FileName[1],Directory)
Netra=ReadLRN(FileName[2],Directory)
Prec=ReadLRN(FileName[3],Directory)


# Funktion, um die erste bzw.letzte 10 Ertraege zu sehen
printData=function(data){
  print("Erste bzw. letzte 10 Werte des Datensatzes")
  print(data$Data[1:10,])
  temp=nrow(data$Data)
  print(data$Data[(temp-10):temp,])
  
}



# Histogram,PDEplot,QQplot,Boxplot,NaN Anzahl
overviewPlot = function(data,name="data"){
  par(mfrow=c(2,3)) 
  #mat <- matrix(c(1,2,0,3,4,5),2,byrow = T)
  #layout(mat,c(5,5,1),c(1,1))
  
  MinD = nanmin(data); MaxD= nanmax(data);
  
  # Histogramm
  histopt(data,name,ylab='Nr of Data');  
  
  # PDE-Dichte
  #PDEplot(data,title=name, xlab=name);  # Bei mir kann das Plot nicht verkleinert werden
  pdeVal = ParetoDensityEstimation(data,ParetoRadius(data))
  plot(pdeVal$kernels,pdeVal$paretoDensity,type='l',xlab='Data',ylab='PDE',main=name)
  #plot(density(data),main=name);
  
  plot.new()
  
  # QQplot
  qqnormfit(data,ylab=name);gridOn(); 
  
  
  
  # Boxplot
  boxplot(data,xlab=name, main=paste('Range:[ ', num2str(round(MinD,5)), ' ,', num2str(round(MaxD,5)) ,' ]'),axes=FALSE);
  
  
  # Barplot fuer die NaNs
  NaNs = (sum(is.infinite(data)) + sum(is.na(data)))/length(data)
  barplot(NaNs, ylab = "NaNs in %", main = paste(round(NaNs, 4), " %"), xlim = c(0, 3), ylim = c(0, 1))
  if (any(is.nan(data), na.rm = TRUE)) 
    print("NaNs in Data found. This message is only important, if after rounding the percent of NaN is zero in the bar plot.")
  if (any(is.infinite(data), na.rm = TRUE)) 
    warning("Infinite values in Data found.")
  
  par(mfrow=c(1,1))
}

  
############################################################################
#NO3
attach(NO3)

# Messungsfrequenz
printData(NO3)

# Ueberblick von Zeitreihe
year=Data[,3]
indy=which(year==2015,arr.ind = T)

plot(Data[,1],type='l',ylab='NO3',xlab='Zeit',col="blue")
vline(min(indy),Color = 'red',LineWidth = 3)


# Verteilungsbetrachtung
#library('reshape2')
#library('caTools')
#par(mfrow=c(2,3))
#InspectVariable(Data[,1],Header[1])  # Muss zuerst Pakete importieren. 
#par(mfrow=c(1,1))                    #Bei mir kann die Graphiken nicht innerhalb ein Bild legen. Version von R?

overviewPlot(Data[,1],Header[1])

# Eliminierung der Nullen und extrem kleinen Werten 
ohneNull=Data[,1][Data[,1]!=0]

overviewPlot(ohneNull,paste(Header[1],"ohne Null"))

# Transformieren Daten, um Normalverteilung zu annaehren
boxcoxTrans(ohneNull[sample(1:length(ohneNull),5000)])

#trans= slog(ohneNull)*100
#trans=slog(1/ohneNull)
trans=(sqrt(ohneNull)*100)
overviewPlot(trans,paste("Transformierte",Header[1]))


# GMM
temp=AdaptGauss(trans)


Werte1=list(Means = c(7.029219, 17.786148, 29.106267), 
           SDs = c(5.08932, 4.50063, 5.69700), 
           Weights = c( 0.176, 0.348, 0.490))
res1=AdaptGauss(trans,Werte1$Means,Werte1$SDs,Werte1$Weights)
dput(res1)
ab=Bayes4Mixtures(trans,Werte1$Means,Werte1$SDs,Werte1$Weights,PlotIt = T) 


#Verifizierung
QQplotGMM(trans,Werte1$Means,Werte1$SDs,Werte1$Weights)

abc=Chi2testMixtures(trans,Werte1$Means,Werte1$SDs,Werte1$Weights,PlotIt = T)


#Clusterung
dec=BayesDecisionBoundaries(Werte1$Means,Werte1$SDs,Werte1$Weights)
Cls=ClassifyByDecisionBoundaries(trans,dec)

ind=which(Data[,1]!=0)
ind2=which(Data[,1]==0)

plot(ind,Data[ind,1],ylab='NO3',xlab='Zeit',col=Cls,main="Clusterung",ylim=c(-0,0.4))
points(ind2,Data[ind2,1],col='blue')
vline(min(indy),Color = "pink",LineWidth = 3)
legend("topright",c("Hoch","Mittel","Niedrig","Null"),col=c("green","red","black","blue"),pch=1)


detach(NO3)

#############################################################################
#Netradiation
printData(Netra)

##############################################################################
#Precip
printData(Prec)



#######################################################################################
###################################################################################
#Aufgabe 2
#Versuchen Sie jede der Variablen mit einem passenden Ansatz zu modellieren. 
#Orientieren sie sich hierbei an der Vorlesungseinheit zur Hydrologie.
########################################################################################
#1. Komponentmodell

#Residualanalysis
showTrend = function(data,trend,trendName){
  index=1:length(data)
  par(mfrow=c(2,2))
  
  plot(data,type="l",col="blue",main=trendName)
  lines(trend$fitted.value,col="green")
  
  res=trend$residuals
  plot(res,type="l",col="red",xlab="Zeit",ylab="Residual",main=paste("Residuals, EW=", round(mean(res),4)))
  abline(h=0)
  
  qqnormfit(res,ylab="Residuals",main="Normal QQplot");gridOn();
  
  plot(res[-1],res[-length(data)], main=paste("Autokorrelation: ",round(cor(res[-1],res[-length(data)]),4)), xlab="X_t", ylab="X_t+1")
  abline(0,1, col="blue",lwd=3)
 
  
  par(mfrow=c(1,1))
  
}

# Wahl des Polynomgrades
polynomGrad=function(data){
  index=1:length(data)
  trend1=lm(data~index)
  fit1.SSE <- sum(resid(trend1)^2)
  for (i in 2:10){
    trend.poly=lm(data~poly(index,i))
    sse <- sum(resid(trend.poly)^2)
    fit1.SSE <- c(fit1.SSE,sse)
  }
  plot(fit1.SSE,type="l",xlab="Grad",ylab="SSE",main="Screen Kriterium")
}


## Hier data einsetzen
polynomGrad(data)

#OptimalTrend: 3th Polynom
trend=lm(data~poly(index,3))
showTrend(data,trend3,"3th Polynom")

#################################################################################################################
#2. ARMA/ARIMA

# Step0. Prüfe zuerst die Stationarität
#plot(data,type="l",col="blue");
#title(paste("Erwartungswert: ", mean(data)))
      
data=NO3$Data[,1]

Differenz=diff(data)
plot(Differenz,type="l",col="green")
title(paste("Differenz, Erwartungswert: ", mean(Differenz)))


# Step1. Modell Postulat : Zeichne die ACF und PACF
acfpacfPlot=function(sample){
  par(mfrow=c(1,2))
  
  acf(sample)
  pacf(sample)
  
  par(mfrow=c(1,1))
}

acfpacfPlot(Differenz)


# Step2. Model Parameter Bestimmung
requireRpackage('forecast')
fit=Arima(Differenz, order=c(0,0,3), include.mean=F)
fit$coef

# Step3. Modell Validierung
# one-step forecast
plot(Differenz,type='l',col='blue')
points(fitted(fit),col='red',type='l')
res2=Differenz-fitted(fit)

# Residuum aktuelle Modell
plot(Differenz,type='l',col='blue')
points(fit$residuals,col='red',type='l')
res=Differenz-fit$residuals

#Residuumanalysis
par(mfrow=c(2,3))
plot(res2, xlab="time",ylab="residuals",main="MA(3) Residuals with one-step forecast")
plot(res2[-length(res2)],res2[-1],main="Autokorelation") 
qqnormfit(res2,ylab="MA(3) Residuals with one-step forecast");gridOn();
 
plot(res,xlab="time",ylab="residuals",main="MA(3) Residuals with Model fitting")
plot(res[-length(res)],res[-1],main="Autokorelation")
qqnormfit(res,ylab="MA(3) Residuals with Model fitting");gridOn();

par(mfrow=c(1,1))

###################################################################################################################
#3. ARCH/GARCH
#Betse AR Modell Schatzen, Berechne die Regression der Yt-1 auf die Yt (KQ) um die Fehler e_t zu erhalten
#Berechne die KQ Regression von Fehler
#Teste die gegenseitige Signifikanz der Parameter
#Falls einige dieser Koeffizienten signifikant konstruiere ein entsprechendes ARCH Modell

acfpacfPlot(Differenz)
#acfpacfPlot(Differenz^2)

requireRpackage('FinTS')
ArchTest(Differenz)


requireRpackage('tseries')
gar1<-garch(Differenz,order=c(1,0))
gar2<-garch(Differenz,order=c(2,0))
gar3<-garch(Differenz,order=c(3,0))
gar4<-garch(Differenz,order=c(1,1))
gar5<-garch(Differenz,order=c(1,2))
gar6<-garch(Differenz,order=c(2,1))
gar7<-garch(Differenz,order=c(2,2))
aic=c()
aic=c(aic,AIC(gar1))
aic=c(aic,AIC(gar2))
aic=c(aic,AIC(gar3))
aic=c(aic,AIC(gar4))
aic=c(aic,AIC(gar5))
aic=c(aic,AIC(gar6))
aic=c(aic,AIC(gar7))
plot(aic,type="l")
points(aic)
legend("topright",paste(1:7,": GARCH(",c("1,0","2,0","3,0","1,1","1,2","2,1","2,2"),")"),pch=1)

summary(gar4)
#Jarque Bera Test: test, ob Residuals normalverteilt ist.
#Residual soll Chi2 verteilt ist --> heavy tail
#Box-Ljung test: test,ob Autokorrelation gibt.
#Residual soll kein Autokorrelation vorliegen

overviewPlot(gar4$residuals,"Residuen für GARCH(1,1)")

plot(gar4)


###############################################################################################################
#4. Hidden Markov Modell
library('RHmm')

#NO3 ohne Nullen hat 3 Gaussian Mixture.
HMM=HMMFit(trans, nStates=3)
HMM$HMM$initProb
HMM$HMM$transMat
HMM$HMM$distribution
summary(HMM)

vit=viterbi(HMM,trans)
HMMGraphicDiag(vit,HMM,trans)



par(mfrow=c(3,1))
plot(vit$states,type="l",main="Clusterung nach Viterbi")
plot(Cls,type="l",main="Clusterung nach GMM")
plot(vit$states-Cls,type="l",main="Diffenz von Clusterung")
par(mfrow=c(1,1))


############################################################################################################
#5. Fast Foriour Transformation




###############################################################################################################
###############################################################################################################
#Aufgabe 3
#Beschreiben Sie den Wissensgewinn durch diese Modelle. Wie und wodurch wird NO3 beeinflusst?
##############################################################################################################



















