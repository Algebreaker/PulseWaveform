# ISO data Model2 Fitting

setwd("Desktop/Scripts")

library(tidyverse)                                 
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils) 
library(spectral)
library(DescTools)

source("~/Desktop/Model2Workshop/model2.R")
source("~/Desktop/Model2Workshop/simplex.R")
source("~/Desktop/Scripts/PPG_funcs.R")

input <- "A"

config.sample_rate = 75.0
lab.time = "time (s)"
lab.ppg = "Detrended"
lab.dt = "interval (s)"

baseline_correct <- 1

# Load PPG
ppg <- read.csv("~/Desktop/Khalsa_Fletcher_collaboration copy/Physio Dial/ISO_3.0/Iso_3_PhysioData/BI506/scan_20190607/physiological_files/ECG_6_ISO_R3.1D", sep = "")   
ppg <- data.frame(
  time = (0:(nrow(ppg)-1)) / 40,
  ppg = ppg[,1]
)
names(ppg)[1] <- lab.time
names(ppg)[2] <- lab.ppg
n <- dim(ppg)[1]
vpg <- ppg[2:n,2] - ppg[1:(n-1),2]
beat <- ppg[which(vpg[1:(n-1)] < 300 & vpg[2:n] >= 300),1]
nBeat <- length(beat)
beat <- data.frame(
  beat = beat,
  dt = (1:nBeat)*0.0
)
rm(vpg)


UnDetrend <- function(ppg,factor=0,offset=1)    # these values could be changed 
{
  k <- offset * (1-factor)
  n <- nrow(ppg)
  result <- (1:n)*0
  result[1] = ppg[1,2]
  # I[i] = O[i] - O[i-1] * factor + I[i-1] - offset * (1 - factor)
  
  for (i in 2:n)
  {
    result[i] = ppg[i,2] - ppg[i-1,2] * factor - k + result[i-1]
  }
  
  return(result)
}

## Check the gradient so as to adjust downtrending parameters:
nBeats <- nrow(beat)
seg <- c(0,0,0)
i <- 1
beatTime <- beat[i,1]
nextTime <- if (i < nBeats){beat[i+1,1]}else{NA}
temp = seg[3]
seg <- model2.FindSegment(ppg,beat[i,1],nextTime)
if (temp > 0){
  seg[1] = temp + 1
}
rm(temp)
data <- model2.GetSegment(ppg,seg)
plot(data)
tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2], data[nrow(data)-2, 2], data[nrow(data)-3, 2], 
          data[nrow(data)-4, 2])    # try making tail just the last 5 points instead of last 9...
tail <- rev(tail)
xx. <- 1:length(tail)
y. <- lm(tail~xx.)
y.[[1]][2] # gradient of the tail


factor_value <- 1.01
while(y.[[1]][2][1] > -20){   #
  
  factor_value <- factor_value - 0.01
  ppg2 <- ppg
  ppg2[, 2] <- UnDetrend(ppg,factor=factor_value,offset=1)
  #plot(ppg2$`time (s)`[1:100], ppg2$Detrended[1:100], type = "l")
  nBeats <- nrow(beat)
  seg <- c(0,0,0)
  beatTime <- beat[i,1]
  nextTime <- if (i < nBeats){beat[i+1,1]}else{NA}
  temp = seg[3]
  seg <- model2.FindSegment(ppg,beat[i,1],nextTime)
  if (temp > 0){
    seg[1] = temp + 1
  }
  rm(temp)
  data <- model2.GetSegment(ppg2,seg)
  plot(data)
  if(y.[[1]][2] > 0){
    tail <- c(data[nrow(data)-5, 2], data[nrow(data)-6, 2], data[nrow(data)-7, 2], data[nrow(data)-8, 2], 
              data[nrow(data)-9, 2])
  }else{
    tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2], data[nrow(data)-2, 2], data[nrow(data)-3, 2], 
              data[nrow(data)-4, 2])
  }
  tail <- rev(tail)
  #plot(tail)
  xx. <- 1:length(tail)
  y. <- lm(tail~xx.)
  y.[[1]][2][1] # gradient of the tail
  
}

ppg3 <- data.frame(ppg[,1],UnDetrend(ppg,factor=factor_value,offset=1))

vv. <- ppg3[, 1]      
yv. <- lm(ppg3[, 2]~vv.)
yv.[[1]][2] # gradient of the tail


offset_value <- 1
while(yv.[[1]][2] > 0){
  if(yv.[[1]][2] > 5){
    offset_value <- offset_value + 1 # was 0.5
  }else{
    if(yv.[[1]][2] > 1){
      offset_value <- offset_value + 0.5
    }else{
      offset_value <- offset_value + 0.05
    }
  }
  ppg3 <- data.frame(ppg[,1],UnDetrend(ppg,factor=factor_value,offset=offset_value))
  vv. <- ppg3[, 1]      
  yv. <- lm(ppg3[, 2]~vv.)
  if(yv.[[1]][2]>0){plot(ppg[,1],UnDetrend(ppg,factor=factor_value,offset=offset_value), type = "l")}
  if(yv.[[1]][2]>0){abline(a = yv.[[1]][1], b = yv.[[1]][2], col = "red")}
  if(yv.[[1]][2]>0){print(yv.[[1]][2])} # gradient of the tail
}


# Now you have a decent decay gradient and a decent offset, set that value to the actual trace:
ppg[,2] = UnDetrend(ppg,factor=factor_value,offset=offset_value)   


# Now correct the baseline:
if(baseline_correct == 1){
  undetrended <- ppg$Detrended
  samplingRate <- 40
  # Run main script code:
  sfunction <- splinefun(1:length(undetrended), undetrended[1:length(undetrended)], method = "natural")
  deriv1 <- sfunction(seq(1, length(undetrended)), deriv = 1)
  spline1 <-  sfunction(seq(1, length(undetrended)), deriv = 0)
  splinePoly <- CubicInterpSplineAsPiecePoly(1:length(undetrended), undetrended[1:length(undetrended)], "natural")
  deriv1Poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended), deriv1, "natural") 
  inflexX <- solve(splinePoly, b = 0, deriv = 1)
  inflexY <- predict(splinePoly, inflexX)
  w <- find_w(d1p = deriv1Poly, deriv1 = deriv1, sp = splinePoly, sr = samplingRate)
  uv <- find_u_v(dat = undetrended, wx = w$wX, wy = w$wY, d1 = deriv1, d1p = deriv1Poly, spline = splinePoly, spline_o = spline1, plot=FALSE)
  o <- find_o(wx = w$wX, inx = inflexX, iny = inflexY, d1p = deriv1Poly, sp = splinePoly)
  tmp <- preclean_wuv(w=w, uv=uv, o=o, samp = samplingRate, sp = spline1, q = F)
  w <- tmp[[1]]
  uv <- tmp[[2]]
  rm(tmp)
  baseCor <- baseline(inx = inflexX, iny = inflexY, o = o, dat = undetrended, sp = splinePoly, plot=F)
  ppg[, 2] <- baseCor
}else{
  ppg[, 2] <- ppg[, 2] - min(ppg[, 2][1:100])
}


# Add output columns

ppg$Baseline = 1:nrow(ppg) * 0
ppg$Excess   = 1:nrow(ppg) * 0
ppg$Residue  = 1:nrow(ppg) * 0
beat$First      = 1:nrow(beat) * 0
beat$Last       = 1:nrow(beat) * 0
beat$Baseline   = 1:nrow(beat) * 0
beat$STime      = 1:nrow(beat) * 0
beat$SAmplitude = 1:nrow(beat) * 0
beat$SWidth     = 1:nrow(beat) * 0
beat$DTime      = 1:nrow(beat) * 0
beat$DAmplitude = 1:nrow(beat) * 0
beat$DWidth     = 1:nrow(beat) * 0
beat$NTime      = 1:nrow(beat) * 0
beat$NAmplitude = 1:nrow(beat) * 0
beat$NWidth     = 1:nrow(beat) * 0


# Estimate Parameters

nBeats <- nrow(beat)
seg <- c(0,0,0)


# Decide on a congif.rate
config.rates <- c()
for(i in 1:10){ #1:nBeats
  config.rate <- 0.95
  beatTime <- beat[i,1]
  nextTime <- if (i < nBeats){beat[i+1,1]}else{NA}
  temp = seg[3]
  seg <- model2.FindSegment(ppg,beat[i,1],nextTime)
  if (temp > 0){
    seg[1] = temp + 1                                
  }
  rm(temp)
  data <- model2.GetSegment(ppg,seg)
  tStart <- ppg[seg[1],1]
  yPrev <- ppg[max(seg[1]-1,1),2]
  amp <- max(data[, 2]) - min(data[, 2])
  constant <- 0.1092254*amp
  baseline <- min(data[,2]) - constant     
  # excess <- model2.Excess(data[,2],yPrev,baseline)
  residue <- model2.Excess(data[,2],ppg[seg[1]-1,2],-0)
  # Compare decay to post-diastolic wave:
  decay_from_s <- data.frame((which.max(data[, 2])+1):nrow(data))       
  if(nrow(decay_from_s) < 10){next}
  for(j in (which.max(data[, 2])+1):nrow(data)){
    decay_from_s[j, 1] <- data[j, 1]
    decay_from_s[j, 2] <- (baseline + config.rate^(j-which.max(data[, 2]))*config.rate*(data[j-(j-6),2]-baseline))
  }
  # last 10 points:
  a. <- data[nrow(data), 2] - decay_from_s[nrow(data), 2] # this is the last point compared to the last decay point
  b. <- data[nrow(data)-1, 2] - decay_from_s[nrow(data)-1, 2]
  c. <- data[nrow(data)-2, 2] - decay_from_s[nrow(data)-2, 2]
  d. <- data[nrow(data)-3, 2] - decay_from_s[nrow(data)-3, 2]
  e. <- data[nrow(data)-4, 2] - decay_from_s[nrow(data)-4, 2]
  f. <- data[nrow(data)-5, 2] - decay_from_s[nrow(data)-5, 2]
  g. <- data[nrow(data)-6, 2] - decay_from_s[nrow(data)-6, 2]
  h. <- data[nrow(data)-7, 2] - decay_from_s[nrow(data)-7, 2]
  i. <- data[nrow(data)-8, 2] - decay_from_s[nrow(data)-8, 2]
  x. <- c(a., b., c., d., e., f., g., h., i.)
  tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2], data[nrow(data)-2, 2], data[nrow(data)-3, 2], 
            data[nrow(data)-4, 2], data[nrow(data)-5, 2], data[nrow(data)-6, 2], data[nrow(data)-7, 2], 
            data[nrow(data)-8, 2])
  tail <- rev(tail)
  xx. <- 1:length(tail)
  y. <- lm(tail~xx.)
  tail_decay <- c(decay_from_s[nrow(data), 2], decay_from_s[nrow(data)-1, 2], decay_from_s[nrow(data)-2, 2], 
                  decay_from_s[nrow(data)-3, 2], decay_from_s[nrow(data)-4, 2], decay_from_s[nrow(data)-5, 2], 
                  decay_from_s[nrow(data)-6, 2], decay_from_s[nrow(data)-7, 2], decay_from_s[nrow(data)-8, 2])
  tail_decay <- rev(tail_decay)
  yy. <- lm(tail_decay~xx.)
  fit <- abs(sd(x.)/mean(x.)) # coefficient of variation, lower percentage = closer fit
  while(fit > 0.05 | a. < 0){
    config.rate <- config.rate - 0.01
    # Compare decay to last 10 points:
    decay_from_s <- data.frame((which.max(data[, 2])+1):nrow(data))
    for(k in (which.max(data[, 2])+1):nrow(data)){
      decay_from_s[k, 1] <- data[k, 1]
      decay_from_s[k, 2] <- (baseline + config.rate^(k-which.max(data[, 2]))*config.rate*(data[k-(k-6),2]-baseline))
    }
    # last 10 points:
    a. <- data[nrow(data), 2] - decay_from_s[nrow(data), 2] # this is the last point compared to the last decay point
    b. <- data[nrow(data)-1, 2] - decay_from_s[nrow(data)-1, 2]
    c. <- data[nrow(data)-2, 2] - decay_from_s[nrow(data)-2, 2]
    d. <- data[nrow(data)-3, 2] - decay_from_s[nrow(data)-3, 2]
    e. <- data[nrow(data)-4, 2] - decay_from_s[nrow(data)-4, 2]
    f. <- data[nrow(data)-5, 2] - decay_from_s[nrow(data)-5, 2]
    g. <- data[nrow(data)-6, 2] - decay_from_s[nrow(data)-6, 2]
    h. <- data[nrow(data)-7, 2] - decay_from_s[nrow(data)-7, 2]
    i. <- data[nrow(data)-8, 2] - decay_from_s[nrow(data)-8, 2]
    x. <- c(a., b., c., d., e., f., g., h., i.)
    subt <- data[(which.max(data[, 2])+1):nrow(data), 2] - decay_from_s[(which.max(data[, 2])+1):nrow(data), 2]
    if(abs(sd(x.)/mean(x.)) > fit & a. > 0 & all(subt > 0) | config.rate < 0.75){break}  # make sure the decay line in under the data line before breaking!
    fit <- abs(sd(x.)/mean(x.))
    print(config.rate)
    plot(data)
    lines(decay_from_s)
  }
  config.rates[i] <- config.rate
}
config.rate <- min(config.rates[!is.na(config.rates)])  # was median rather than min


####### Make the pulse dataframe: ###########
baseCor <- ppg$Detrended
sfunctionBC <- splinefun(1:length(baseCor), baseCor, method = "natural")
deriv1BC <- sfunctionBC(seq(1, length(baseCor)), deriv = 1)
spline1BC <- sfunctionBC(seq(1, length(baseCor)), deriv = 0)
splinePolyBC <- CubicInterpSplineAsPiecePoly(1:length(baseCor), baseCor, "natural")
deriv1PolyBC <- CubicInterpSplineAsPiecePoly(1:length(baseCor), deriv1BC, "natural") 
w$wY <- predict(splinePolyBC, w$wX)
uv$uY <- predict(splinePolyBC, uv$uX)
uv$vY <- predict(splinePolyBC, uv$vX)
wuv <- cbind(w, uv)
tmp <- clean_wuv(wuv = wuv, sp = splinePolyBC, inx = inflexX, o = o, samp = samplingRate, bc = baseCor, q = F)
wuv <- tmp[[1]]
ibi <- tmp[[2]]
oDiff <- tmp[[3]]
rm(tmp, w, uv)
waveLen <- round(median(oDiff)+15) 
tmp <- sep_beats(odiff = oDiff, bc = baseCor, samp = samplingRate, wuv = wuv, wvlen = waveLen, ibi=ibi, o=o, inx = inflexX, scale = T, q = F) 
pulse <- tmp[[2]]
avWave <- tmp[[1]]
rm(tmp)
# Finding OSND on the average wave:
tmp <- diast_pk(avw = avWave, sr = samplingRate, scale = T)
dPeak <- tmp[1]
xShift <- tmp[2]
rm(tmp)
osnd <- osnd_of_average(avWave, dp = dPeak, diff = 0, sr = samplingRate)
# Calculate the distance from S to D (Dtime):
D_time <- ((osnd$x[4] - osnd$x[2])/10)/samplingRate
################ finished making pulse dataframe ####################

# Find the index value that corresponds to each beat time
beat_timings_in_index_form <- c()
for(i in 1:nrow(beat)){
  beat_timings_in_index_form[i] <- which(ppg$`time (s)` == beat[i, 1])
}
beats_that_correspond_to_waves_in_pulse <- c()
for(i in 1:length(wuv$wX)){
  beats_that_correspond_to_waves_in_pulse[i] <- which.min(abs(wuv$wX[i] - beat_timings_in_index_form))
}
# Now subset 'beat' so only the beats that correspond to waves in pulse are kept:
beat <- beat[beats_that_correspond_to_waves_in_pulse, ]

########## If not incorporating the pulse dataframe, continue from here instead ##############


# Find the Excess, and the three peaks:
nBeats <- nrow(beat)
seg <- c(0,0,0)

for(i in 1:10){ #1:nBeats
  
  beatTime <- beat[i,1]
  nextTime <- if(i < nBeats){beat[i+1,1]}else{NA}
  
  temp = seg[3]
  seg <- model2.FindSegment(ppg,beat[i,1],nextTime)
  if (temp > 0){
    seg[1] = temp + 1
  }
  rm(temp)
  
  data <- model2.GetSegment(ppg,seg)
  tStart <- ppg[seg[1],1]
  yPrev <- ppg[max(seg[1]-1,1),2]
  
  amp <- max(data[, 2]) - min(data[, 2])
  constant <- 0.1092254*amp
  baseline <- min(data[,2]) - constant    
  
  residue <- model2.Excess(data[,2], ppg[seg[1]-1,2], -0)   
  
  count <- nrow(data)
  excess <- 1:count * 0.0
  excess[1] = data[1,2] - (baseline + config.rate*(yPrev-baseline))
  for (j in 2:count){
    excess[j] = data[j,2] - (baseline + config.rate*(data[j-1,2]-baseline))
  }
  rm(count)
  rm(j)
  plot(data[,1],excess, type = "l")
  
  par <- 1:10 * 0.
  par[1] = baseline
  
  residue <- excess  
  
  
  # S peak
  
  peak.w <- which(data[,1] > beat[i,1]-0.2 & data[,1] < beat[i,1]+0.2)  
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  #plot(data[,1],excess)   
  #lines(peak.t,peak.y)
  par[3] <- max(peak.y)      # height of the S peak
  par[2] <- peak.t[which(peak.y==par[3])]   # timing of the S peak
  par[4] <- 0.25    # width of the S-peak, a priori 
  #
  rm(peak.w,peak.t,peak.y)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[2:4])
  
  #plot(data[,1],excess)
  #lines(data[,1],residue)  # this plots the residual after S peak has been removed. 
  
  
  # D peak
  
  peak.w <- which(data[,1] > beat[i,1]+0.3 & data[,1] < beat[i,1]+0.6)  # this finds a range in where to find the peak of d
  peak.t <- data[peak.w,1]     # this finds the time corresponding to the peak
  peak.y <- residue[peak.w]
  #plot(data[,1],excess)   
  #lines(peak.t,peak.y)
  #
  par[6] <- max(peak.y)      # height of the D peak
  par[5] <- peak.t[which(peak.y==par[6])]    # Timing of the D peak
  par[7] <- 0.25    # Width of the D peak
  #
  rm(peak.w,peak.t,peak.y)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[5:7])
  
  #plot(data[,1],excess)
  #lines(data[,1],residue)
  
  
  # N peak
  
  t <- par[2] + c(0.25,0.75) * (par[5]-par[2])
  peak.w <- which(data[,1] > t[1] & data[,1] < t[2])
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  #plot(data[,1],excess)   
  #lines(peak.t,peak.y)
  #
  par[9] <- max(peak.y)
  par[8] <- peak.t[which(peak.y==par[9])]
  par[10] <- 0.25
  #
  rm(peak.w,peak.t,peak.y,t)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[8:10])
  
  #plot(data[,1],excess)
  #lines(data[,1],residue)
  
  
  # Store parameters
  
  w <- seg[1]:seg[3]
  ppg$Baseline[w] <- baseline
  ppg$Excess[w] <- excess
  ppg$Residue[w] <- residue
  rm(w,excess,residue,data,baseline,yPrev,nextTime,tStart)
  
  beat[i,3:4]  = c(seg[1],seg[3])
  beat[i,5:14] = par
  beat[i,9] = beat[i,9]-beat[i,6]
  beat[i,12] = beat[i,12]-beat[i,6]
  rm(par)
}
rm(seg)


# Refine Parameters (downhill simplex)      

for (i in 1:10){ #1:nBeats
  print(i)
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  
  par <- beat[i,5:14]
  #par[5] = par[5]-par[2]
  #par[8] = par[8]-par[2]
  sim <- simplex.MakeSimplex(data[,1:2],par,model2.ChiSq,0.1) 
  sim <- simplex.Run(data[,1:2],sim,model2.ChiSq)  
  par <- sim[1,]
  
  sim <- simplex.MakeSimplex(data[,1:2],sim[1,],model2.ChiSq,0.1)
  sim <- simplex.Run(data[,1:2],sim,model2.ChiSq)     
  rm(data,par)
  
  beat[i,5:14] = sim[1,]
  rm(sim)
}

# Refine the parameters again, following parameter reduction as in run_model2
fixedParams <- c(
  median(beat$SWidth[1:10]),  # S width
  median(beat$DTime[1:10]),  # D timing
  median(beat$DWidth[1:10]),  # D width
  median(beat$NTime[1:10]), # N timing
  median(beat$NWidth[1:10])  # N width
)

# Fix outputs of the simplex:
if(fixedParams[4] > 0.3){fixedParams[4] <- 0.3}
if(fixedParams[4] < 0.1){fixedParams[4] <- 0.1}
# (Optional) Make D_timing the d from the average:
fixedParams[2] <- D_time


# Re-running the downhill simplex, this time with the fixed parameters:
for(i in 1:10){ #1:nBeats
  print(i)
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  
  par <- as.double(beat[i,c(5:7,10,13)])
  sim <- simplex.MakeSimplex(data[,1:2],par,model2.ChiSqAmp,0.1,optional=fixedParams)
  sim <- simplex.Run(data[,1:2],sim,model2.ChiSqAmp,optional=fixedParams)
  par <- sim[1,]
  
  sim <- simplex.MakeSimplex(data[,1:2],sim[1,],model2.ChiSqAmp,0.1,optional=fixedParams)
  sim <- simplex.Run(data[,1:2],sim,model2.ChiSqAmp,optional=fixedParams)
  rm(data,par)
  
  beat[i,5:14] = c(sim[1,1:3],fixedParams[1:2],sim[1,4],fixedParams[3:4],sim[1,5],fixedParams[5])
  rm(sim)
}


# Plot a given wave model against its raw data form:
for(i in 1:10){
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  plot(data[,1],(data[,2]), ylim = c(-1000, 2000)) # ylim = c(-1000, 2000)
  temp<-model2.Rebuild(data,0,beat[i,5:14],TRUE)    # 5-14 are all the relevant parameters
  lines(data[,1], temp, col = "black")
  temp<-model2.Rebuild(data,0,as.double(beat[i,5:8]), TRUE)  # systolic only wave 
  lines(data[,1],temp, col = "red") 
  temp<-model2.Rebuild(data,0,as.double(c(beat[i,5],beat[i,6]+beat[i,9],beat[i,10:11])),TRUE) # diastolic
  lines(data[,1],temp, col = "blue") 
  temp<-model2.Rebuild(data,0,as.double(c(beat[i,5],beat[i,6]+beat[i,12],beat[i,13:14])),TRUE) # renal
  lines(data[,1],temp, col = "green") 
  #lines(c(1,100), c(0, 0))
  lines(c(1, 100), rep(beat[i, 5], 2))
}
