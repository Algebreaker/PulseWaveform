setwd("xxxxxx")

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
undetrended <- preproc(dat=data)
samplingRate <- 75
rm(data)

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
w <- find_w(d1p = deriv1Poly, deriv1 = deriv1, sp = splinePoly, sr = samplingRate)

# Find U and V:
uv <- find_u_v(dat = undetrended, wx = w$wX, wy = w$wY, d1 = deriv1, d1p = deriv1Poly, spline = splinePoly, spline_o = spline1, plot=FALSE)

# Find O in order to find the baseline
o <- find_o(wx = w$wX, inx = inflexX, iny = inflexY, d1p = deriv1Poly, sp = splinePoly)
#plot(spline1[1:500], type = "l")
#points(inflexX[o], inflexY[o])

# Preclean:
tmp <- preclean_wuv(w=w, uv=uv, o=o, samp = samplingRate, sp = spline1, q = F)
w <- tmp[[1]]
uv <- tmp[[2]]
rm(tmp)

########################################################################    

#             Step 3 : Correct Baseline and refind W, U, V and O       #

######################################################################## 

# Correct Baseline
baseCor <- baseline(inx = inflexX, iny = inflexY, o = o, dat = undetrended, sp = splinePoly, plot=F)

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
tmp <- clean_wuv(wuv = wuv, sp = splinePolyBC, inx = inflexX, o = o, samp = samplingRate, bc = baseCor, q = F)
wuv <- tmp[[1]]
ibi <- tmp[[2]]
oDiff <- tmp[[3]]
rm(tmp, w, uv)

########################################################################    

#       Step 4 : Find individual waves and the average wave            #

########################################################################

# Find the average length of a wave (and 15 since we are starting the wave from before O):
waveLen <- round(median(oDiff)+15) 

tmp <- sep_beats(odiff = oDiff, bc = baseCor, samp = samplingRate, wuv = wuv, wvlen = waveLen, ibi=ibi, o=o, inx = inflexX, scale = T, q = F) 
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
  theme(legend.position = "none") + labs( y= "PPG Output", x = "Time (Seconds)") + xlim(-0.25, max(average$x)) # + xlim(c(pulse$x[max(which(is.na(avWave)))], pulse$x[length(avWave)])) + ylim(range(avWave[!is.na(avWave)]*1.5)) 

# Create a polynomial spline for each wave:
polyWave <- list()
for(i in 2:ncol(pulse)){
  polyWave[[i-1]] <-CubicInterpSplineAsPiecePoly(pulse$x, pulse[, i], "natural")
}


########################################################################    

#                     Step 5 : Find O, S, N and D                      #

########################################################################


tmp <- diast_pk(avw = avWave, sr = samplingRate, scale = F)
dPeak <- tmp[1]
xShift <- tmp[2]
rm(tmp)

# Find OSND, then extend range as necessary (class 3 + waveforms only)
osnd <- osnd_of_average(avWave, dp = dPeak, diff = 0, sr = samplingRate)
if(dPeak == 5*samplingRate){
  dPeak <- osnd$x[4]*1.2 
}
# artefacts causing sharp inclines can create a step effect which can be falsely detected as a notch
# if N and D are too close together e.g < 1.5, reduce the diastolic peak value
#  most of these artefacts are at the end of the average wave when waves start to drop off
if((osnd$x[4]-osnd$x[3]) < 1.5 & (osnd$x[4]-osnd$x[3]) > 0){
  dPeak <- dPeak*0.95
  osnd <- osnd_of_average(avWave, dp = dPeak, diff = 0)
}

## Find OSND for each individual wave:
# Make a list of OSND for each individual wave in pulse:
scale <- 1   # Set to 1 if scaling
osnd_all <- list()
for(i in 2:ncol(pulse)){  #ncol(pulse)
  wavi <- pulse[, i][!is.na(pulse[, i])]
  if(scale == 1){
    xShift2 <- (which(abs(wavi - 0.5) == min(abs(wavi - 0.5))))  # the new 0.5 
  }else{
    xShift2 <- which.min(abs(wavi))
  }
  diff <- xShift - xShift2
  dpa <- dPeak - diff
  osnd_all[[i-1]] <- osnd_of_average(aw = wavi, dp = dpa, diff = diff, sr = samplingRate)
}

# Can plot all OSND values against the average to see if there are any obvious anomalies:
plot(avWave[!is.na(avWave)], type = "l")
for(i in 1:length(osnd_all)){
  points(osnd_all[[i]][3, 1], osnd_all[[i]][3, 2], col = "red")
}


########################################################################    

#                     Step 6 : Feature Extraction                      #

########################################################################

# First: Correct all OSNDs so that O = 0
for(i in 1:length(osnd_all)){
  osnd_all[[i]]$y <- osnd_all[[i]]$y - osnd_all[[i]]$y[1]
}

# Extract features:
features <- feature_extract(oa = osnd_all, p = pulse, pw = polyWave)
