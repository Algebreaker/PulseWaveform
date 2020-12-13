setwd("xxxxx")

library(tidyverse)                                 
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils) #SplinesUtils is best downloaded directly from Github
library(spectral)
library(seewave)

source("PPG_funcs.R")

########################################################################    

#             Step 1 : Reading in Data and Preprocessing               #

######################################################################## 


#### Bioradio (finger) ####
data <- read.csv("Source.csv", header = T)
#downsampling data and undetrending 
undetrendedData <- data.frame(preproc(dat=data))
undetrended <- undetrendedData$undetrended
samplingRate <- 75

#### Mimic database (sample) ####
#data <- read.csv("~/Desktop/data.csv", header = T)
#data <- data[, c(1, 6)]
#undetrendedDat <- data
#colnames(undetrendedDat)[2] <- "undetrended"
#undetrended <- undetrendedDat$undetrended
#plot(undetrended[1:10000], type = "l")
#undetrended <- undetrended[1:10000, ]
#samplingRate <- 60


#### ISO-3 (finger) #### 
#data <- read.delim(file.choose(), header = T)
#data <- cbind(data, (1:nrow(data)))
#undetrended <- data
#colnames(undetrended)[1] <- "undetrended"
#samplingRate <- 40


#### Intero-battery (ear) ####
#data <- read.csv("~/Desktop/Khalsa_Fletcher_collaboration copy/Intero_Battery/AN Float/data-organized/AT053/T0/behavioral_session/AT053-T0-__BH-R1-PHYS.csv")
#data <- data[, c(1, 4)]
#data <- data[seq(1, nrow(data), 100), ]
#undetrendedDat <- data
#colnames(undetrended)[2] <- "undetrended"
#undetrended <- undetrendedDat$undetrended
# Most ear PPG traces have large periods of nothing followed by large artefact followed by waveforms, remove the first two:
#undetrended <- undetrended[-c(1:(max(which(undetrended == max(undetrended)))+1000)), ]  
#samplingRate <- 20

########################################################################    

#             Step 2 : Create splines, find W, U, V and O              #

######################################################################## 


# Create discrete splines:
sfunction <- splinefun(1:length(undetrended), undetrended[1:length(undetrended)], method = "natural")
deriv1 <- sfunction(seq(1, length(undetrended)), deriv = 1)
spline1 <-  sfunction(seq(1, length(undetrended)), deriv = 0)

# Create polynomial splines: 
splinePoly <- CubicInterpSplineAsPiecePoly(1:length(undetrended), undetrended[1:length(undetrended)], "natural")
deriv1Poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended), deriv1, "natural") 

# Finding inflexion points on splinePoly (points on deriv1 where x = 0):
inflexX <- solve(splinePoly, b = 0, deriv = 1)
inflexY <- predict(splinePoly, inflexX)

# Find W:
w <- find_w(d1p = deriv1Poly, deriv1 = deriv1, sp = splinePoly)

# Find U and V:
uv <- find_u_v(dat = undetrended, wx = w$wX, wy = w$wY, d1 = deriv1, d1p = deriv1Poly, spline = splinePoly, spline_o = spline1, plot=FALSE)

# Find O in order to find the baseline

o <- find_o(wx = w$wX, inx = inflexX, iny = inflexY, d1p = deriv1Poly, sp = splinePoly)
#plot(spline1[1:500], type = "l")
#points(inflexX[o], inflexY[o])


# Find where W lies from U to V on the x-axis, for each wave find the percentage distance to V: 
vDist <- c()
for(i in 1:length(w$wX)){
  vDist[i] <- (w$wX[i] - uv$uX[i]) / (uv$vx[i] - uv$uX[i])*100
}
#plot(w$wX, vDist, type = "l")

# Make a vector of abnormal pdtv (allowing any values between 35 and 65): 
sdpdtv <- sd(vDist)
pdtvWaves <- c(which(vDist > (sdpdtv + median(vDist)) & vDist > 65), which(vDist < (median(vDist) - sdpdtv) & vDist < 35))

# Remove pdtvWaves:
if(length(pdtvWaves) > 0){
  w <- w[-pdtvWaves, ]
  uv <- uv[-pdtvWaves, ]
}

########################################################################    

#             Step 3 : Correct Baseline and refind W, U, V and O       #

######################################################################## 

# Correct Baseline
baseCor <- baseline(inx = inflexX, iny = inflexY, o = o, dat = undetrended, sp = splinePoly, plot=T)

# Redefine discrete splines:
sfunctionBC <- splinefun(1:length(baseCor), baseCor, method = "natural")
deriv1BC <- sfunctionBC(seq(1, length(baseCor)), deriv = 1)
spline1BC <- sfunctionBC(seq(1, length(baseCor)), deriv = 0)

# Redefine polynomial splines: 
splinePolyBC <- CubicInterpSplineAsPiecePoly(1:length(baseCor), baseCor, "natural")
deriv1PolyBC <- CubicInterpSplineAsPiecePoly(1:length(baseCor), deriv1BC, "natural") 

# Refind W (y-values):
w$wY <- predict(splinePolyBC, w$wX)

# Refind U and V (y_values):
uv$uY <- predict(splinePolyBC, uv$uX)
uv$vY <- predict(splinePolyBC, uv$vX)

wuv <- cbind(w, uv)
tmp <- clean_wuv(wuv = wuv, sp = splinePolyBC, inx = inflexX, o = o)
wuv <- tmp[[1]]
ibi <- tmp[[2]]
oDiff <- tmp[[3]]
rm(tmp, w, uv)

########################################################################    

#       Step 4 : Find individual waves and the average wave            #

########################################################################
# Find the average length of a wave (and 15 since we are starting the wave from before O):
waveLen <- round(median(oDiff)+15) 

tmp <- sep_beats(odiff = oDiff, bc = baseCor, samp = samplingRate, wuv = wuv, wvlen = waveLen) 
pulse <- tmp[[2]]
avWave <- tmp[[1]]
rm(tmp)

# Can now stack and plot the mean +/- median wave 

pulse_stacked <- gather(pulse, key = "wave_ID", value = "values", -c("x"))

average <- data.frame(seq((-141/(samplingRate*10)), ((waveLen*10 -9)-142)/(samplingRate*10), by = 1/(samplingRate*10)))  
average <- cbind(average, avWave)
colnames(average)[1] <- "x"

ggplot(data = pulse_stacked, aes(x, values, col = wave_ID), col = "black") +
  scale_color_manual(values = rep("black", ncol(pulse))) +  
  geom_line(size = 1.5, alpha = ((1/length(wuv$wX)*10)-(1/length(wuv$wX)))) + geom_line(data = average, aes(x, avWave), size = 1.125, color = "red") +  # ylim will vary based on source data    
  theme(legend.position = "none") + labs( y= "PPG Output", x = "Time (Seconds)")  # + xlim(c(pulse$x[max(which(is.na(avWave)))], pulse$x[length(avWave)])) + ylim(range(avWave[!is.na(avWave)]*1.5)) 

# Create a polynomial spline for each wave:
polyWave <- list()
for(i in 2:ncol(pulse)){
  polyWave[[i-1]] <-CubicInterpSplineAsPiecePoly(pulse$x, pulse[, i], "natural")
}


########################################################################    

#                     Step 5 : Find O, S, N and D                      #

########################################################################

## Find OSND on the average wave:
# Find the diastolic peak on the average wave to inform OSND finding (also some adjusment of x-values for removal of NA values):
avWave <- avWave[!is.na(avWave)]
# Need to find new W position (0.5) after removing NAs
x_shift1 <- which(abs(avWave-0.5) == min(abs(avWave - 0.5)))
avWavePoly <- CubicInterpSplineAsPiecePoly(1:length(avWave), avWave, "natural")
inflexion_points_av <- solve(avWavePoly, b = 0, deriv = 1)
inflexion_points_av_yval <- predict(avWavePoly, inflexion_points_av)
#plot(avWavePoly)
#points(inflexion_points_av, inflexion_points_av_yval)
# Specify limitations for where the diastolic peak can first be found i.e between 120:230 on x-axis, and below 1 on y-axis:
peaks <- order(inflexion_points_av_yval[which(inflexion_points_av < 215 & inflexion_points_av > 120 & inflexion_points_av_yval < 1)], decreasing = TRUE)
diastolic_peak <- inflexion_points_av[which(inflexion_points_av < 215 & inflexion_points_av > 120 & inflexion_points_av_yval < 1)][peaks[1]]
# diastolic_peak will be NA for class 3 waveforms, in which case set a default value
if(is.na(diastolic_peak) | diastolic_peak < inflexion_points_av[peaks[1]]){
  diastolic_peak <- 10*samplingRate
}
# Find OSND, then extend range as necessary (class 3 + waveforms only)
osnd <- osnd_of_average(avWave, dp = diastolic_peak, diff = 0)
if(diastolic_peak == 5*samplingRate){
  diastolic_peak <- osnd$x[4]*1.2 
}
# artefacts causing sharp inclines can create a step effect which can be falsely detected as a notch
# if N and D are too close together e.g < 1.5, reduce the diastolic peak value
#  most of these artefacts are at the end of the average wave when waves start to drop off
if((osnd$x[4]-osnd$x[3]) < 1.5 & (osnd$x[4]-osnd$x[3]) > 0){
  diastolic_peak <- diastolic_peak*0.95
  osnd <- osnd_of_average(avWave, dp = diastolic_peak, diff = 0)
}


## Find OSND for each individual wave:
# Make a list of OSND for each individual wave in pulse:
osnd_all <- list()
for(i in 2:ncol(pulse)){  #ncol(pulse)
  wave_no_nas <- pulse[, i][!is.na(pulse[, i])]
  x_shift2 <- (which(abs(wave_no_nas - 0.5) == min(abs(wave_no_nas - 0.5))))  # the new 0.5 
  diff <- x_shift1 - x_shift2
  dpa <- diastolic_peak - diff
  osnd_all[[i-1]] <- osnd_of_average(aw = wave_no_nas, dp = dpa, diff = diff)
}

# Can plot all OSND values against the average to see if there are any obvious anomalies:
for(i in 1:length(osnd_all)){
  points(osnd_all[[i]][4, 1], osnd_all[[i]][4, 2])
}


##########################################################################################################



#Use new polynomial splines to find w/u/v/notch values for each waveform 
change____var____name <- find_wuv(p=pulse, col_len = waveLen, p_w = polyWave)

## Find O, S, N, D on the new polynomial splines:
osnd_xy <- find_osnd(p = pulse, p_w = polyWave, col_len = waveLen, wuvn = wuv)
osnd_y <- osnd_xy[1:(length(osnd_xy)/2)]
osnd_x <- osnd_xy[(length(osnd_xy)/2+1):length(osnd_xy)]


#Fit S and D sines using the OSND points

sd_sines <- find_sd_sine(p = pulse, wuvn = wuv, osndx = osnd_x, osndy = osnd_y, pw = polyWave, plot=FALSE)
s_sines <- sd_sines[1:(length(sd_sines)/2)]
d_sines <- sd_sines[(length(sd_sines)/2+1):length(sd_sines)]

#Fit N sines using the S and D sines

n_sines <- fit_n_sine(p=pulse, ss = s_sines, ds = d_sines, osndx = osnd_x, osndy = osnd_y, wuvn = wuv, plot=TRUE)


##Refitting the SND peaks 

avg_period <- 3*(mean(wuv$v_x)-mean(wuv$u_x))
refitted_snd <- refit_peaks(p=pulse, ss= s_sines, ds= d_sines, ns= n_sines, period= avg_period)

##Plotting the refitted SND peaks back on the waveforms 
#for(i in 2:(ncol(pulse2))){
#  plot(1:nrow(pulse2), pulse2[,i], type='l')
#  points(refitted_snd$s_x[i-1], refitted_snd$s_y[i-1], pch = 19, col ='red')
#  points(refitted_snd$n_x[i-1], refitted_snd$n_y[i-1], pch = 19, col = 'blue')
#  points(refitted_snd$d_x[i-1], refitted_snd$d_y[i-1], pch = 19, col = 'yellow')
#}

#Creating a new osnd_x/y list incorporating the reffited SND peaks 
refit_osnd_x <- list()
refit_osnd_y <- list()

for(i in 1:length(osnd_x)){
  refit_osnd_x[[i]] <- c(NA, NA, NA, NA)
  refit_osnd_x[[c(i,1)]] <- osnd_x[[c(i,1)]]
  refit_osnd_x[[c(i,2)]] <- refitted_snd$s_x[i]
  refit_osnd_x[[c(i,3)]] <- refitted_snd$n_x[i]
  refit_osnd_x[[c(i,4)]] <- refitted_snd$d_x[i]
  refit_osnd_y[[i]] <- c(NA, NA, NA, NA)
  refit_osnd_y[[c(i,1)]] <- osnd_y[[c(i,1)]]
  refit_osnd_y[[c(i,2)]] <- refitted_snd$s_y[i]
  refit_osnd_y[[c(i,3)]] <- refitted_snd$n_y[i]
  refit_osnd_y[[c(i,4)]] <- refitted_snd$d_y[i]
}



#Refit sines using the refitted SND values
refit_sd_sines <- find_sd_sine(p = pulse, wuvn = wuv, osndx = refit_osnd_x, osndy = refit_osnd_y, pw = polyWave, plot=FALSE)
refit_s_sines <- refit_sd_sines[1:(length(sd_sines)/2)]
refit_d_sines <- refit_sd_sines[(length(sd_sines)/2+1):length(sd_sines)]

refit_n_sines <- fit_n_sine(p=pulse, ss = s_sines, ds = d_sines, osndx = refit_osnd_x, osndy = refit_osnd_y, wuvn = wuv, plot=FALSE)

#Finding the residual after subtracting all the sines from the original pulse 
sine_resid <- list()
const <- c()
for(i in 2:ncol(pulse)){
  const <- rep(refit_s_sines[[i-1]][[length(refit_s_sines[[i-1]])]],length(refit_s_sines[[i-1]]))
  sine_resid[[i-1]] <- pulse[,i]
  sine_resid[[i-1]] <- (((((sine_resid[[i-1]] - refit_s_sines[[i-1]]) + const) - refit_d_sines[[i-1]]) + const) - refit_n_sines[[i-1]]) + const
}

#Plotting all the sines on the pulse trace, as well as the residual
#for(i in 2:ncol(pulse2)){
#  plot(1:(nrow(pulse2)), pulse2[,i], type = 'l')
#  lines(1:(nrow(pulse2)), refit_s_sines[[i-1]], col='red')
#  lines(1:(nrow(pulse2)), refit_d_sines[[i-1]], col='purple')
#  lines(1:(nrow(pulse2)), refit_n_sines[[i-1]], col='green')
#  lines(1:(nrow(pulse2)), sine_resid[[i-1]], col='yellow')
#}


#Correct all OSNDs so that O = 0
for(i in 1:length(osnd_y)){
  osnd_correction <- osnd_y[[i]]
  osnd_correction <- osnd_correction - osnd_correction[1]
  osnd_y[[i]] <- osnd_correction
  print(osnd_y[[i]])
}




## Features (canonical waveform first):

# Peak to notch time (waveforms need to be normalized on o-o interval first (or can divide by the pulse duration)):
pn_time <- c()
for(i in 1:length(osnd)){
  pn_time[i] <- (x_osnd[[c(i, 3)]] - x_osnd[[c(i, 2)]])/(next_o[i] - x_osnd[[c(i, 1)]])
}

# Peak to peak time (PPT):
PPT <- c()
for(i in 1:length(osnd)){
  PPT[i] <- x_osnd[[c(i, 4)]] - x_osnd[[c(i, 2)]]
}

# Stiffness Index = height / PPT
SI <- c()
for(i in 1:length(osnd)){
  SI[i] <- # Insert height here # / PPT[i]
}

# Notch to peak ratio  = height of notch / height of primary peak (this was referred to as RI in multivariate paper)
np_ratio <- c()
for(i in 1:length(osnd)){
  np_ratio[i] <-  osnd[[c(i, 3)]] / osnd[[c(i, 2)]]
}

# Notch-time ratio = time interval from notch to end of pulse / time interval from notch to beginning of pulse
nt_ratio <- c()
for(i in 1:length(osnd)){
  nt_ratio[i] <- (next_o[i] - x_osnd[[c(i, 3)]]) /  (x_osnd[[c(i, 3)]] - x_osnd[[c(i, 1)]])
}

# Systolic amplitude (aka amplitude):
sa <- c()
for(i in 1:length(osnd)){
  sa[i] <-  osnd[[c(i, 2)]]
}

# Reflectance peak to forward peak ratio (Augmentation index as per Takazawa et al, Reflection index as per Padilla et al) )
AI <- c()
for(i in 1:length(osnd)){
  AI[i] <-  osnd[[c(i, 4)]] / osnd[[c(i, 2)]]
}

# Alternate augmentation index (Rubins et al) ( (x-y)/x ):
aAI <- c()
for(i in 1:length(osnd)){
  aAI[i] <- (osnd[[c(i, 2)]] - osnd[[c(i, 4)]]) / osnd[[c(i, 2)]]
}

# Crest time (time from foot of waveform (o) to peak (S)):
CT <- c()
for(i in 1:length(osnd)){
  CT[i] <- x_osnd[[c(i, 2)]] - x_osnd[[c(i, 1)]]
}

# Inflexion point area ratio (For canonical waveform use c(i, 3), if using inflection point then use c(i, 4)):
ipa_ratio <- c()
for(i in 1:(length(polyWave)-1)){
  f <- c(x_osnd[[c(i, 1)]]:x_osnd[[c(i, 3)]], x_osnd[[c(i, 3)]])
  g <- predict(polyWave[[i]], f)
  j <- c(x_osnd[[c(i, 3)]]:next_o[i], next_o[i])
  k <- predict(polyWave[[i]], j)
  auc_systole <- AUC(f, g, method = "spline")
  if(sum(j < 0) > 0){
    poly_wave_diastole_subset <-CubicInterpSplineAsPiecePoly(j, k, "natural")
    zero_crossing <- solve(poly_wave_diastole_subset, b = 0) 
    zero_crossing_yval <- predict(poly_wave_diastole_subset, zero_crossing)
    first_zero_crossing <- zero_crossing[min(which(zero_crossing > (x_osnd[[c(i, 1)]] + 30)))]   # first element that crosses 0 on the waves descent
    j <- c(x_osnd[[c(i, 3)]]:first_zero_crossing, first_zero_crossing)
    k <- predict(polyWave[[i]], j)
    auc_diastole <- AUC(j, k, method = "spline")
  }else{
    auc_diastole <- AUC(j, k, method = "spline")
  }
  ipa_ratio[i] <- auc_diastole / auc_systole
  print(auc_diastole + auc_systole)
}




##Spectral Analysis

spectralanalysis <- (spectrum(baseCor))

