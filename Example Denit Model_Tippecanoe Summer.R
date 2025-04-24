### Script File for modeling single station metabolism and open-channel denitrification based upon MIMS data and Bayesian models using rjags
### written by AJR
### Last updated 9 Oct 2015
### Everything has been updated with new priors, as well as the new modeling approach from Hall et al. 2015
### On 9 Oct 2015 I updated so that we can use three replicate samples from each time point in the jags models 
rm(list=ls())
#First we set our working directory - this needs to be modified to wherever your data are located
setwd("C:/Users/abaga/Desktop/for.r/diel")

#We have two data files - the first has data for the single station metabolism model measured using a sonde (e.g., DO only)
#Call in oxygen, temp, bp, light, and depth data:
oxydata<-read.csv(file='tippe.onestation.sonde.7.20.22.csv', header=T)

#The second file has data for the gas ratio models
#call in MIMS data (N2toAr measured and O2toAr Measured) and diel sampling data (temp, time)
diel.data.replicated<-read.csv(file="Tippe.mims.repdata.7.20.21.csv", header=T)

#We also need to open some packages in order to be able to use certain functions:
library(rjags)
library(chron)

#Combine dates and times to get continuous date/time column
dtime<-chron(dates=as.character(oxydata$date), times=as.character(oxydata$time))

oxydata$dtime<-dtime

#We pulled bp from a weather station, but that needs to be corrected for elevation:
###temp is degC, alt is m, and bpst is in inches of Hg.  Temp is usually relative to a standard, 15 degC.  This is from Colt's book
bpcalc<- function(bpst, alt, temp) {
  
  bpst*25.4*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+temp)))
  
}

#Calculate elevation corrected bp: bp from the weather station were in in Hg, and need to be converted to mm Hg, and altitude needs to be in m but is entered in ft. Temp is degrees C
oxydata$bpcalc<-bpcalc(oxydata$bp*0.0393701, alt=814*0.3084, temp=oxydata$temp)

### Calculate saturation for O2, Ar and N2 outside of the model
### These are really complex equations and are easier to understand outside
### of the JAGS model, rather than having to type out all of these exponents


###First calculate water density at different temps (from Patterson and Morris 1994 Meterologia):
watdens<-function(temp){
  
  t<-temp
  
  A <- 7.0132e-5
  B <- 7.926295e-3 
  C <-  -7.575477e-5 
  D<- 7.314701e-7
  E <-  -3.596363e-9
  to<- 3.9818
  
  dens<- (999.97358- (A*(t-to) + B*(t-to)^2 +C*(t-to)^3 + D*(t-to)^4+E*(t-to)^5) ) -4.873e-3 + 1.708e-4*t - 3.108e-6 * t^2
  dens/1000
}

##We have two different functions for oxygen saturation, both from Garcia and Gordon 1992 L&O, calculates mL gas / dm^3 water; the below equation converts to mg/L. 
osat.sonde<- function(temp, bp) {
  ts<-log((298.15-temp) / (273.15 + temp))
  a0<-2.00907
  a1<-3.22014
  a2<-4.0501
  a3<-4.94457
  a4<- -0.256847
  a5<- 3.88767
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  sato<-(exp(a0 + a1*ts + a2*ts^2 + a3*ts^3 + a4*ts^4+a5*ts^5))*((bp-u)/(760-u))*1.42905
  
}

###alternative formulation from Garcia and Gordon (umol/kg), which is converted to mg/L and corrected for density.  This gives the same values as from Colt and is the one that should be used when dealing with MIMS data.
osat<- function(temp, bp) {
  u<-10^(8.10765-(1750.286/(235+temp)))
  ts<-log((298.15-temp) / (273.15 + temp))
  a0<-5.80871
  a1<-3.20291
  a2<-4.17887
  a3<-5.1006
  a4<- -9.86643e-2
  a5<-3.88767
  sato<-(exp(a0 + a1*ts + a2*ts^2 + a3*ts^3 + a4*ts^4)+a5*ts^5)*((bp-u)/(760-u))
  watdens(temp)*sato*(31.9988/1000)###converts umol/kg to mg/L
  
}

#Compare the two functions to see they give slightly different results
x<-osat.sonde(25,760)
x
osat(25,760)

###Next, based upon solute type (N2 or Ar), river temp, salinity 
###(always 0 for our model), and bp, calculate N2 and Ar expected
###equilibrium concentration based upon Hamme and Emerson 2004

# Equation (1) from Hamme and Emerson calculates gas concentration at equilibrium with atmosphere:

# ln(C) = A0 + A1Ts + A2Ts^2 + A3TS^3 + S (BO + B1Ts + B2Ts^2)
# C = gas concentration at equilibrium with a moist atmosphere at 1-atm (umol/kg), t is the temperature, and S is the salinity

#Equation (2) is the temperature correction term:
# Ts = ln((298.15-t / 273.15 + t))

#We will use these equations to establish saturation for N2 and Ar by writing
#our own function:


#solute = N2 or Ar, temp is in deg C, and barpress is elevation-corrected mm Hg
getSatn<-function(solute, temp, salinity, barpress){
  Ts=log((298.15-temp)/(273.15+temp))
  if(solute=="N2"){
    A0=6.42931
    A1=2.92704
    A2=4.32531
    A3=4.69149
    B0=-0.00744129
    B1=-0.00802566
    B2=-0.0146775
  }
  else{
    A0=2.79150
    A1=3.17609
    A2=4.13116
    A3=4.90379
    B0=-0.00696233
    B1=-0.00766670
    B2=-0.0116888
  }
  ln.C = A0 + A1 * Ts + A2 * Ts^2 + A3 * Ts^3 + salinity*(B0 + B1 * Ts + B2 * Ts^2)
  u<-10^(8.10765-(1750.286/(235+temp)))
  C<-exp(ln.C)*((barpress-u)/(760-u))
  if(solute=="N2"){
  converted<-C*watdens(temp)*(28.014/1000) # converts N2 from umolN2/kg to mgN/L
  }
  else{
    converted<-C*watdens(temp)*(39.948/1000)# converts Ar from umol/kg to mg/L
  }
  return(converted)
}


#Check the function
getSatn("N2", 15,0, barpress=760)/(28.014/1000) # Check - if getSatn is correct it should return 589.0272 umol/L N



#Now the traditional one-station metabolism model
#First calculate oxygen saturation using the sonde function at each time point
#Next we have to set up some vectors to feed data to jags
oxySat<-(osat.sonde(oxydata$temp, oxydata$bpcalc))
#And other data needed for the model
N<-length(oxydata$oxy) # Number of timepoints
z<-mean(oxydata$z) #Average river depth
bp<-mean(oxydata$bpcalc) # bp at each timepoint
light<-oxydata$light #light at each timepoint
temp<-oxydata$temp #temp at each timepoint
y<-oxydata$oxy #measured DO (mg/L) at each time point from sonde

#Make the call to jags:
metabolism<-jags.model('SingleStationSondeModel.bug',#you will need to update this directory with wherever you saved the .bug file containing the jags model script. If you have it in the same directory as the data, you only need to include the file name and extension here.
                       data=list('N'=N, #all of these data are things not defined within the jags script or unknowns
                                 'y'=y,
                                 'z'=z,
                                 'temp'=temp,
                                 'light'=light,
                                 'oxySat'=oxySat),
                       n.chains=4, #Number of different parallel chains to run
                       n.adapt=1000) #Number of iterations for adaptation

#Now run the model to get posterior estimates
metab.results<-coda.samples(metabolism, 
                            c('GPP', 'ER', 'K'),
                            100000)
#Print the posterior mean and credible quantiles
summary(metab.results) 
#Paste Results Below (if you run the model multiple times you will get *slightly* different results, but generally not different up to the hundredths place):

#Iterations = 1001:101000
#Thinning interval = 1 
#Number of chains = 4 
#Sample size per chain = 1e+05 

#1. Empirical mean and standard deviation for each variable,
#plus standard error of the mean:
  
#  Mean     SD  Naive SE Time-series SE
#ER  -12.223 0.24547 0.0003881      0.0034272
#GPP   5.145 0.08908 0.0001409      0.0009098
#K     4.180 0.10541 0.0001667      0.0014333

#2. Quantiles for each variable:
  
#  2.5%     25%     50%     75%   97.5%
#ER  -12.715 -12.387 -12.219 -12.056 -11.754
#GPP   4.972   5.084   5.144   5.204   5.321
#K     3.978   4.108   4.178   4.250   4.391
#Model Convergence Diagnostics - you can visually estimate convergence by looking at the trace of the model runs and the distribution - the gelman plot and diagnostic should be at or very close to 1 to insure convergence
windows()
plot(metab.results)
gelman.diag(metab.results)
gelman.plot(metab.results)

GPP.sonde<-5.14 # Mean of credible interval
ER.sonde<--12.22 # Mean of credible interval
K.sonde<-4.18 # Mean of credible interval #20.66 for prec

#We write this function to set up a plotting routine later
o2modeled<-function(oxy, z, light, temp, GPP, ER, K){
  O2.mass.predicted<-rep(NA, length(oxy))
  O2.mass.predicted[1]<-oxy[1]
  for(i in 2:length(O2.mass.predicted)){
    O2.mass.predicted[i]<-(O2.mass.predicted[i-1]+((GPP/z)*(((light[i]+light[i-1])*0.5)/sum(light)))+ ER*0.006944/z+ (K/(600/(1800.6-(temp[i]*120.1)+(3.7818*temp[i]^2)-(0.047608*temp[i]^3)))^-0.5) * 0.0069444*((oxySat[i-1]+oxySat[i]-O2.mass.predicted[i-1])/2))/(1+((K/(600/(1800.6-(temp[i]*120.1)+(3.7818*temp[i]^2)-(0.047608*temp[i]^3)))^-0.5)*0.006944)/2)
  }
  O2.mass.predicted
}

#Call the function from above
sonde.metabolism.output<-o2modeled(oxy=oxydata$oxy, z=0.71, light=oxydata$light, temp=oxydata$temp, GPP=GPP.sonde, ER=ER.sonde, K=K.sonde)

#Set up results of each parameter into a dataframe to allow posterior density plots (see below)
metab.results.matrix<-as.matrix(metab.results)
metab.results.df<-as.data.frame(metab.results.matrix)
#We will use the specific parameter results for plotting below
metab.results.er<-as.numeric(metab.results.df[,1])
metab.results.gpp<-as.numeric(metab.results.df[,2])
metab.results.k<-as.numeric(metab.results.df[,3])

#That's it for the single station sonde model

###Next we will modify some of the raw data for the diel gas ratio model

#Combine dates and times to get continuous date/time column
dtime<-chron(dates=as.character(diel.data.replicated$Date), times=as.character(diel.data.replicated$Time))

diel.data.replicated$dtime<-dtime

#Make a new column in the data that has the length of time between timepoints - this allows for different times between timepoints, but our data are consistently 1 h apart
timesteps<-rep(0,length(diel.data.replicated$dtime))
for(i in 2:length(timesteps)){
  timesteps[i]=diel.data.replicated$dtime[i]-diel.data.replicated$dtime[i-1]
}

diel.data.replicated$TimeStep<-timesteps

# Calculate elevation corrected bp the same way as for the metabolism model

# here we have bp in the raw data as in and mm Hg, so we don't have to correct to mm Hg
diel.data.replicated$bpcalc<-bpcalc(diel.data.replicated$bp, alt=814*0.3084, temp=diel.data.replicated$Temp)

#We also have to calculate the measured N, Ar, and O2 from the ratios by multiplying the ratios by Ar Sat
#First convert molar ratios to mass ratios for each replicate:
#First for N2:Ar
diel.data.replicated$N2toArMassRatioRep1<-(diel.data.replicated$Rep1N2toArMeasured*(28.014/1000))/(39.948/1000)
diel.data.replicated$N2toArMassRatioRep2<-(diel.data.replicated$Rep2N2toArMeasured*(28.014/1000))/(39.948/1000)
diel.data.replicated$N2toArMassRatioRep3<-(diel.data.replicated$Rep3N2toArMeasured*(28.014/1000))/(39.948/1000)
#And for O2:Ar
diel.data.replicated$O2toArMassRatioRep1<-(diel.data.replicated$Rep1O2toArMeasured*(31.988/1000))/(39.948/1000)
diel.data.replicated$O2toArMassRatioRep2<-(diel.data.replicated$Rep2O2toArMeasured*(31.988/1000))/(39.948/1000)
diel.data.replicated$O2toArMassRatioRep3<-(diel.data.replicated$Rep3O2toArMeasured*(31.988/1000))/(39.948/1000)
#Next calculate N and Ar saturation at each timepoint (O2 saturation will be calculated later)
diel.data.replicated$NSat<-getSatn(solute="N2", temp=diel.data.replicated$Temp, salinity=0, barpress=diel.data.replicated$bpcalc)
diel.data.replicated$ArSat<-getSatn(solute="Ar", temp=diel.data.replicated$Temp, salinity=0, barpress=diel.data.replicated$bpcalc)
#We need initial values for the actual dissolved N, so we multiple the N2:Ar ratio by Ar saturation to get dissolved N2 gas
diel.data.replicated$Rep1NMeasured<-(diel.data.replicated$N2toArMassRatioRep1*diel.data.replicated$ArSat)
diel.data.replicated$Rep2NMeasured<-(diel.data.replicated$N2toArMassRatioRep2*diel.data.replicated$ArSat)
diel.data.replicated$Rep3NMeasured<-(diel.data.replicated$N2toArMassRatioRep3*diel.data.replicated$ArSat)


#Next we run the denitrification MCMC
#First we have to set some parameters:

N<-length(diel.data.replicated$N2toArMassRatioRep1) # number of time points
z<-0.71 # depth in m (from previous knowledge)
timestep<-diel.data.replicated$TimeStep # length of time (days) between samples
NSat<-diel.data.replicated$NSat #N2 (mg/L) expected at 100% saturation
ArSat<-diel.data.replicated$ArSat # Ar (mg/L) expected at 100% saturation
y1<-diel.data.replicated$N2toArMassRatioRep1 #N2/Ar measured on MIMS
y2<-diel.data.replicated$N2toArMassRatioRep2 #N2/Ar measured on MIMS
y3<-diel.data.replicated$N2toArMassRatioRep3 #N2/Ar measured on MIMS
temp<-diel.data.replicated$Temp # Temp at each timepoint
Ninit<-mean(c(diel.data.replicated$Rep1NMeasured[1], diel.data.replicated$Rep2NMeasured[1], diel.data.replicated$Rep3NMeasure[1])) # Average of 3 initial dissolved N2 reps


#The model and all of the units are in grams, meters, and days, output is g/m2/d
#Initialize the model
denit<-jags.model('N2toArMassBalanceModel_7.20.22.tippe.bug',
                  data= list('N' = N,
                             'y1'=y1,
                             'y2'=y2,
                             'y3'=y3,
                             'z'=z,
                             'NSat'=NSat,
                             'ArSat'=ArSat,
                             'Ninit'=Ninit,
                             'temp'=temp,
                             'timestep' = timestep),
                  n.chains=4,
                  n.adapt=1000)

#And run the full model
denit.results<-coda.samples(denit,
                            c('DN', 'K'),
                            100000)

windows()
plot(denit.results)
gelman.diag(denit.results)
gelman.plot(denit.results)
summary(denit.results)
#Paste Results
#Iterations = 1001:101000
#Thinning interval = 1 
#Number of chains = 4 
#Sample size per chain = 1e+05 

#1. Empirical mean and standard deviation for each variable,
#plus standard error of the mean:
  
#  Mean      SD  Naive SE Time-series SE
#DN  0.0227 0.01369 2.165e-05      0.0001103
#K  0.2802 0.16828 2.661e-04      0.0013648

#2. Quantiles for each variable:
  
#  2.5%     25%     50%     75%   97.5%
#DN 0.001624 0.01215 0.02128 0.03148 0.05283
#K  0.024845 0.15302 0.25894 0.38268 0.66253

################################

#Model results are in gN/m2/d

################################

### CONVERT TO mgN/m2/h to compare to microcosms and other literature
DN.mean<-.0227*(1000/24)
DN.upper<-.05283*(1000/24)
DN.lower<-.001624*(1000/24)

DN.mean #.9458
DN.upper #2.2
DN.lower #.06767

DN.model<-.0227
K.dn<-.2802

denit.results.matrix<-as.matrix(denit.results)
denit.results.DN<-as.numeric(denit.results.matrix[,1])*(1000/24)
denit.results.k<-as.numeric(denit.results.matrix[,2])

#Now let's look at O2:Ar for metabolism:
#Calculate O2 Measured:
O2Measured1<-diel.data.replicated$O2toArMassRatioRep1*diel.data.replicated$ArSat
O2Measured2<-diel.data.replicated$O2toArMassRatioRep2*diel.data.replicated$ArSat
O2Measured3<-diel.data.replicated$O2toArMassRatioRep3*diel.data.replicated$ArSat
o2Sat<-(osat(diel.data.replicated$Temp, diel.data.replicated$bpcalc))

N<-length(diel.data.replicated$O2toArMassRatioRep1)
z<-0.71 # depth in m (from previous knowledge)
timestep<-diel.data.replicated$TimeStep # length of time (days) between samples
ArSat<-diel.data.replicated$ArSat # Ar (mg/L) expected at 100% saturation
y1<-diel.data.replicated$O2toArMassRatioRep1
y2<-diel.data.replicated$O2toArMassRatioRep2
y3<-diel.data.replicated$O2toArMassRatioRep3#O2/Ar (mass ratio) measured on MIMS
temp<-diel.data.replicated$Temp
bp<-diel.data.replicated$bpcalc
light<-diel.data.replicated$light
o2init<-mean(c(O2Measured1[1], O2Measured2[1], O2Measured3[1]))

#Initialize the model. Units are g, m, and days, with output as gO2/m2/d for GPP and ER, and /d for K
metab.mims<-jags.model('O2toArMetabModel_7.20.22.tippe.bug',
                  data= list('N' = N,
                             'y1'=y1,
                             'y2'=y2,
                             'y3'=y3,
                             'z'=z,
                             'ArSat'=ArSat,
                             'o2Sat'=o2Sat,
                             'o2init'=o2init,
                             'light'=light,
                             'temp'=temp,
                             'timestep' = timestep),
                  n.chains=4,
                  n.adapt=1000)
#Run the model
metab.mims.results<-coda.samples(metab.mims,
                            c('GPP','ER', 'K'),
                            100000)

windows()
plot(metab.mims.results)
gelman.diag(metab.mims.results)
gelman.plot(metab.mims.results)

summary(metab.mims.results)
#Iterations = 1001:101000
#Thinning interval = 1 
#Number of chains = 4 
#Sample size per chain = 1e+05 

#1. Empirical mean and standard deviation for each variable,
#plus standard error of the mean:
  
#  Mean     SD  Naive SE Time-series SE
#ER  -8.930 0.4990 0.0007891       0.006625
#GPP  3.593 0.1631 0.0002580       0.001390
#K    3.396 0.2348 0.0003712       0.003061

#2. Quantiles for each variable:
  
#  2.5%    25%    50%    75%  97.5%
#ER  -9.946 -9.259 -8.919 -8.588 -7.987
#GPP  3.280  3.483  3.591  3.702  3.920
#K    2.955  3.235  3.391  3.550  3.876
GPP.mims<-3.593
ER.mims<--8.93
K.mims<-3.396

metab.mims.results.matrix<-as.matrix(metab.mims.results) #converts coda list to dataframe
#We will use the specific parameter results for plotting below
metab.mims.results.er<-(as.numeric(metab.mims.results.matrix[,1]))
metab.mims.results.gpp<-(as.numeric(metab.mims.results.matrix[,2]))
metab.mims.results.k<-as.numeric(metab.mims.results.matrix[,3])

o2armodeled<-function(o2Sat, ArSat, light, temp, z, timestep, GPP, ER, K, o2init){
O2modeled<-rep(NA, length(o2Sat))
O2modeled[1]<-o2init
for(i in 2:length(O2modeled)){
  O2modeled[i]<-(O2modeled[i-1]+((GPP/z)*(((light[i]+light[i-1])*0.5)/sum(light)))+ ER*timestep[i]/z+ (K/(600/(1800.6-(temp[i]*120.1)+(3.7818*temp[i]^2)-(0.047608*temp[i]^3)))^-0.5) * timestep[i]*((o2Sat[i-1]+o2Sat[i]-O2modeled[i-1])/2))/(1+((K/(600/(1800.6-(temp[i]*120.1)+(3.7818*temp[i]^2)-(0.047608*temp[i]^3)))^-0.5)*timestep[i])/2)
}

Armodeled.O2<-rep(NA, length(ArSat))
Armodeled.O2[1]<-ArSat[1]
for(i in 2:length(Armodeled.O2)){
  Armodeled.O2[i]<-(Armodeled.O2[i-1]+ (K/(600/(1759.7-(temp[i]*117.37)+(3.6959*temp[i]^2)-(0.046527*temp[i]^3)))^-0.5) * timestep[i]*((ArSat[i-1]+ArSat[i]-Armodeled.O2[i-1])/2))/(1+((K/(600/(1759.7-(temp[i]*117.37)+(3.6959*temp[i]^2)-(0.046527*temp[i]^3)))^-0.5)*timestep[i])/2)
}

O2toArModeled<-O2modeled/Armodeled.O2
}

mims.model.output<-o2armodeled(o2Sat=o2Sat, ArSat=diel.data.replicated$ArSat, light=diel.data.replicated$light, temp=diel.data.replicated$Temp, z=0.71, timestep=diel.data.replicated$TimeStep, GPP=GPP.mims, ER=ER.mims, K=K.mims, o2init=o2init)

#tiff("Figures/O2ModelFit.tiff", height=7, width=3.5, units='in', res=300)
windows(height=7, width=3.5)
par(mfcol=c(2,1),mar=c(1,4,1,1), oma=c(3,0,0,0), family='serif')
plot(oxydata$dtime, oxydata$oxy, pch=16, xaxt='n', yaxt='n', ylim=c(4,8), bty='o', xlab='', ylab='', col='gray70')
lines(oxydata$dtime, sonde.metabolism.output, lwd=2, col='black')
axis(1, c(oxydata$dtime[9], oxydata$dtime[45], oxydata$dtime[81], oxydata$dtime[117]), labels=c("","", "", "" ), cex=1.1)
axis(2, c(4,5,6,7,8), las=2, cex=0.9)
mtext(expression(paste("Dissolved O"[2], " (g m"^-3, ")")), side=2, line=2.25, cex=1.4)
legend("topright", legend=c("Observed", "Modeled"), lty=c(NA, 1), lwd=c(NA, 2), pch=c(16,NA), col=c("gray70", "black"), bty='n')
#text(oxydata$dtime[9], 9.1, "(a)", cex=1.5)

plot(rep(diel.data.replicated$dtime,3), c(diel.data.replicated$O2toArMassRatioRep1, diel.data.replicated$O2toArMassRatioRep2, diel.data.replicated$O2toArMassRatioRep3), pch=16, ylim=c(9,20), bty='o', xaxt='n', yaxt='n', xlab='', ylab='', col='gray70')
axis(1, c(oxydata$dtime[9], oxydata$dtime[45], oxydata$dtime[81], oxydata$dtime[117]), labels=c("12:00","18:00", "00:00", "06:00" ), cex=0.9)
axis(2, c(9,10,11,12,13,14,15,16,17,18,19,20), las=2, cex=0.9)
lines(diel.data.replicated$dtime, mims.model.output, col='black', lwd=2)
mtext(expression(paste({O}[2], ":Ar")), side=2, line=2.6, cex=1.4)
mtext("Time of Day (hh:mm)", side=1, line=2.4, cex=1.4)
#text(oxydata$dtime[9], 15.2,"(b)", cex=1.5)
dev.off()

#Now for DN

narmodeled<-function(ArSat, NSat, temp, z, timestep, DN, K, Ninit){
Armodeled<-rep(NA, length(ArSat))
Armodeled[1]<-ArSat[1]
for(i in 2:length(Armodeled)){
  Armodeled[i]<-(Armodeled[i-1]+ (K/(600/(1759.7-(temp[i]*117.37)+(3.6959*temp[i]^2)-(0.046527*temp[i]^3)))^-0.5) * timestep[i]*((ArSat[i-1]+ArSat[i]-Armodeled[i-1])/2))/(1+((K/(600/(1759.7-(temp[i]*117.37)+(3.6959*temp[i]^2)-(0.046527*temp[i]^3)))^-0.5)*timestep[i])/2)
}


N2withdenit<-rep(NA, length(NSat))
N2withdenit[1]<-Ninit

for(i in 2:length(N2withdenit)){
  N2withdenit[i]<-(N2withdenit[i-1]+ DN*timestep[i]/z+ (K/(600/(1970.7-(temp[i]*131.45)+(4.139*temp[i]^2)-(0.052106*temp[i]^3)))^-0.5) * timestep[i]*((NSat[i-1]+NSat[i]-N2withdenit[i-1])/2))/(1+((K/(600/(1970.7-(temp[i]*131.45)+(4.1390*temp[i]^2)-(0.052106*temp[i]^3)))^-0.5)*timestep[i])/2)
}

N2toAr<-N2withdenit/Armodeled
}

DNmodeled<-narmodeled(ArSat=diel.data.replicated$ArSat, NSat=diel.data.replicated$NSat, temp=diel.data.replicated$Temp, z=0.71, timestep=diel.data.replicated$TimeStep, DN=DN.model, K=K.dn, Ninit=Ninit)
DNmodeledNoDenit<-narmodeled(ArSat=diel.data.replicated$ArSat, NSat=diel.data.replicated$NSat, temp=diel.data.replicated$Temp, z=0.71, timestep=diel.data.replicated$TimeStep, DN=0, K=K.dn, Ninit=NSat[1])

windows(height=3.5, width=3.5)
#tiff(file="Figures/ModelTest.tiff", height=3.5, width=3.5, units='in', res=300)
par(oma=c(0,0,0,0), mar=c(4,4,0.5,0.5), family='serif')
plot(rep(diel.data.replicated$dtime,3), c(diel.data.replicated$N2toArMassRatioRep1,diel.data.replicated$N2toArMassRatioRep2,diel.data.replicated$N2toArMassRatioRep3), pch=16, col='black', ylim=c(27,31), xaxt='n', yaxt='n', bty='o', xlab='', ylab='')
points(diel.data.replicated$dtime, NSat/ArSat, pch=16, col='gray')
lines(diel.data.replicated$dtime, DNmodeled, lty=2, lwd=2)
lines(diel.data.replicated$dtime, DNmodeledNoDenit, lty=2, lwd=2, col='gray')
axis(1, c(diel.data.replicated$dtime[2], diel.data.replicated$dtime[8], diel.data.replicated$dtime[14], diel.data.replicated$dtime[20]), labels=c("12:00","18:00", "00:00", "06:00" ), cex.axis=1)
axis(2, c(27,28,29,30,31), las=2, cex.axis=1)
mtext(expression(paste("N"[2], ":Ar")),side=2,cex=1.6,line=2.6)
mtext("Time of day (hh:mm)", side=1, cex=1.6, line=2.6)
legend("topright", legend=c(expression(paste({N}[2], ":Ar Observed")), expression(paste({N}[2], ":Ar Equilibrium")), "Model Fit (With DN)", "Model Fit (No DN)"), bty='n', pch=c(16,16,NA,NA), lty=c(NA,NA,2,2), lwd=c(NA,NA,2,2),col=c("black", "gray", "black", "gray"), cex=0.9)
rect(xleft=diel.data.replicated$dtime[11], ybottom=26.71, xright=diel.data.replicated$dtime[20], ytop=26.72, col='black')
dev.off()


###Horizontal 5 panel plot of posteriors:
# first set up the densities:
dens.gpp<-density(metab.results.gpp)
dens.gpp.mims<-density(metab.mims.results.gpp)
dens.er<-density(metab.results.er)
dens.er.mims<-density(metab.mims.results.er)
dens.k<-density(metab.results.k)
dens.k.mims<-density(metab.mims.results.k)
dens.k.dn<-density(denit.results.k)
dens.dn<-density(denit.results.DN)



#Horizontal
windows(width=9, height=3)
#tiff(file="Figures/Posteriors.tiff", width=9, heigh=3, units='in', res=300)
par(mfrow=c(1,4), mar=c(5.5,0.5,1,1), oma=c(0,2.5,0,0), family='serif')

plot(dens.gpp$x, dens.gpp$y/max(dens.gpp$y), type='l', bty='o', xaxt='n', yaxt='n', xlim=c(0.38,0.81), ylim=c(0,1),xlab='', ylab='', lwd=2, col='black', lty=2)
lines(dens.gpp.mims$x, dens.gpp.mims$y/max(dens.gpp.mims$y), col='gray60', lwd=2)
#axis(2, c(0,10,20,30), cex.axis=1.5, las=2)
axis(1, c(0.4,0.5,0.6,0.7,0.8), cex.axis=1.5)
#abline(v=mean(dens.gpp$x), lty=3, col='gray80')
#abline(v=mean(dens.gpp.mims$x), lty=3, col='black')
mtext(expression(paste("GPP (g O"[2], " m"^-2, " d"^-1, ")")), side=1, line=3.5, cex=1.3)
mtext("Probability Density", side=2, line=1, cex=1.2)
text(0.38,0.95,"(a)", cex=1.7, pos=4)

plot(dens.er$x,dens.er$y/max(dens.er$y), type='l', bty='o',lty=2, col='black', xlab='',yaxt='n',  ylab='', xaxt='n', xlim=c(-4, -2), lwd=2, ylim=c(0,1))
lines(dens.er.mims$x, dens.er.mims$y/max(dens.er.mims$y), col='gray60', lwd=2)
#axis(2,c(0,1,2,3,4,5), cex.axis=1.5, las=2)
axis(1, c(-4,-3.5,-3,-2.5,-2), cex.axis=1.5)
#abline(v=mean(dens.er$x), lty=3, col='gray80')
mtext(expression(paste("ER (g O"[2], " m"^-2, " d"^-1, ")")), side=1, line=3.5, cex=1.3)
text(-4,0.95,"(b)", cex=1.7,pos=4)

plot(dens.k$x, dens.k$y/max(dens.k$y), type='l', bty='o', xlab='',xlim=c(2.5, 4.5), ylim=c(0,1), xaxt='n', yaxt='n', lty=2, col='black', ylab='', lwd=2)
#axis(2,c(0,1,2,3,4), cex.axis=1.5, las=2)
axis(1, c(2.5,3,3.5,4,4.5), cex.axis=1.5)
lines(dens.k.dn$x, dens.k.dn$y/max(dens.k.dn$y), col='black', lwd=2)
lines(dens.k.mims$x, dens.k.mims$y/max(dens.k.mims$y), col='gray60', lwd=2)
#abline(v=mean(dens.k.dn$x), lty=3, col='gray80')
mtext(expression(paste(italic({K}[600]), "(d"^-1,")")), side=1, line=3.5, cex=1.3)
text(2.5,0.95,"(c)", cex=1.7,pos=4)

plot(dens.dn$x, dens.dn$y/max(dens.dn$y), type='l', bty='o', ylim=c(0,1), xlim=c(6,12), xaxt='n', yaxt='n', xlab='', ylab='', lwd=2)
#axis(2, c(0,0.2,0.4,0.6,0.8,1), cex.axis=1.5, las=2)
axis(1, c(6,8,10,12), cex.axis=1.5)
#abline(v=mean(dens.dn$x), lty=3, col='gray80')
mtext(expression(paste("DN (mg N m"^-2, " h"^-1, ")")), side=1, line=3.5, cex=1.3)
text(6,0.95,"(d)", cex=1.7,pos=4)
legend(9,0.9, legend=c("Base Met.", "Ratio Met.", "Ratio DN"), lty=c(2,1,1), col=c("black", "gray60", "black"), bty='n', lwd=2, cex=1)

dev.off()


##Denit:Resp
d.r<-DN.mean/ER.mims
abs(d.r) #.1
  
tipsondesumdata<-oxydata
tipsondesumdata$sondemetabmodel<-sonde.metabolism.output

tipmimssumdata<-diel.data.replicated
tipmimssumdata$dnmodel<-DNmodeled
tipmimssumdata$nodnmodel<-DNmodeledNoDenit
tipmimssumdata$mimsoxy<-mims.model.output


write.csv(tipsondesumdata, file="tipsondesumdata.csv")
write.csv(tipmimssumdata, file="tipmimssumdata.csv")



#Create own plots in different r code