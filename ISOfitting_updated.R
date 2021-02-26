# Model2 (Updated for across time series fitting) for ISO data:

# See new functions below before starting

setwd("Desktop/Scripts")

library(tidyverse)                                 
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils) 
library(spectral)
library(DescTools)

source("~/Desktop/WorkshopCS/Workshop2/model2.R")
source("~/Desktop/WorkshopCS/Workshop2/simplex.R")
source("~/Desktop/Scripts/PPG_funcs.R")

input <- "A"

config.sample_rate = 40 # was 40 for ISO
lab.time = "time (s)"
lab.ppg = "Detrended"
lab.dt = "interval (s)"

baseline_correct <- 1



ppg <- read.csv("~/Desktop/Khalsa_Fletcher_collaboration copy/Physio Dial/ISO_3.0/Iso_3_PhysioData/AA343/scan_20190905/physiological_files/ECG_6_ISO_R1.1D", sep = "")   
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


UnDetrend <- function(ppg,factor=0,offset=1)    
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

# Look at original trace:
plot(ppg$`time (s)`[1:100], ppg$Detrended[1:100], type = "l")

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

# Now adjust factor until gradient = 0
factor_value <- 1.01
while(y.[[1]][2][1] > -10){     # default = -20  
  
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
  #abline(a = y.[[1]][1], b = y.[[1]][2])
  y.[[1]][2][1] # gradient of the tail
  
}

# Now you have a decent factor value (decay), you need a decent offset:
#plot(ppg[,1],UnDetrend(ppg,factor=factor_value,offset=1), type = "l") 
ppg3 <- data.frame(ppg[,1],UnDetrend(ppg,factor=factor_value,offset=1))

vv. <- ppg3[, 1]      
yv. <- lm(ppg3[, 2]~vv.)
#abline(a = yv.[[1]][1], b = yv.[[1]][2], col = "red")
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


# Show data
plot(ppg[1:1000,1],ppg[1:1000,2],t='l')
points(beat[,1],(1:nBeats)*0+84)



# Define how many beats you want to fit over:
beats_in <- 10  


# Find the Excess, and the three peaks:
nBeats <- nrow(beat)
seg <- c(0,0,0)

for(i in 1:100){ #1:beats_in    # 1:100    # c(1:71, 73:101)
  
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
  baseline <- min(data[,2]) - constant     # this 0.3 constant is for bioradio only, will need changing for Tulsa data
  # excess <- model2.Excess(data[,2],yPrev,baseline)
  
  residue <- model2.Excess(data[,2], ppg[seg[1]-1,2], -0)      # Sometimes this line triggers an error if seg[1] = 1
  
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


beat <- add_column(beat, Baseline2 = beat$Baseline, .after = 5)
beat <- cbind(beat, rep(0, nrow(beat)))
colnames(beat)[16] <- "config.rate"


# Start the for loop here:

beat_orig <- beat

fit_check <- list()
for(k in 1:10){         
  
beat <- beat_orig[((k*10)-9):(k*10), ]

# Extract the median NTime:
beat <- beat[1:beats_in, ]    # or whatever the number of rows of estimated parameters is 
renal_param <- median(beat$NTime)


beat <-  beat[1:beats_in, ]   
beat[, 16] <- rep(0.75, beats_in)
nBeats <- nrow(beat)
par <- as.numeric(beat[1,5:16])   

beat_start <- beat[, 3]
beat_end <- beat[, 4]
beat_vector <- list(beat_start, beat_end)
beat_vector <- c(beats_in, beat_vector)



####################################     Refine Parameters     ####################################

##################################  MAKE SIMPLEX 1  ##################################

# Make simplex for the within-beat parameters (10 currenlty):
a <- list()
for(i in 1:beats_in){                      
  
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  
  par <- as.numeric(beat[i,5:16])
  par <- model2.FixParams3(data[, 1:2], par)
  
  a[[i]] <- simplex.MakeSimplex3(data[,1:2],par,model2.ChiSq,0.1) 
}

# Make simplex for the across-beat parameters:
sim <- simplex.MakeSimplex2(data=ppg, param = par, f = model2.ChiSq3, inScale = 0.1, beat_vector = beat_vector, beat = beat, renal_param = renal_param)

# Combine results from makesimplex2 and makesimplex 3 into matrix:
sim <- make_matrix(sim, a)

######################################################################################

##################################  RUN SIMPLEX 1  ###################################

sim <- simplex.Run2(data = ppg, simplexParam = sim, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, renal_param = renal_param, run = "run 1")

# Check best result:    
#model2.ChiSq3(data = ppg, params = NULL, beats = beat_vector, beat = beat, a = sim[1, ], plot = TRUE, renal_param = renal_param)

# Extract the across and within beat parameters outputted by the simplex:
across <- sim[1, ][1:6]
within <- list()
for(i in 1:beats_in){
  temp <- rep(0, 12)
  temp[c(1:4, 7, 10)] <-  sim[1, ][((i*6)+1):((i*6)+6)]  
  within[[i]] <- temp
}

#######################################################################################

##################################  FIX PARAMS 1  #####################################

# Fix Params:
fixed <- list()
for(i in 1:beats_in){
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  fixed[[i]] <- model2.FixParams3(data, params = as.numeric(within[[i]]), across_beat_params = across)
}


# Update beat with fixed values:
new_beat <- data.frame(matrix(0, ncol = 12, nrow = beats_in))
for(i in 1:beats_in){
  new_beat[i, ] <- fixed[[i]]
}
# Across / within reassigned (having now fixed):
across <- as.numeric(new_beat[1, c(5:6, 8:9, 11:12)])
within <- list()
for(i in 1:beats_in){
  temp <- rep(0, 12)
  temp[c(1:4, 7, 10)] <- as.numeric(new_beat[i, ][c(1:4, 7, 10)])
  within[[i]] <- temp
}
# Add the first four columns from beat:
for(i in 1:4){
  new_beat <- cbind(beat[, i], new_beat)
}


#######################################################################################

##################################  MAKE SIMPLEX 2  ###################################

# Make simplex for the within-beat parameters:
a <- list()
for(i in 1:beats_in){                
  
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  
  witin <- within[[i]]                                               # You shouldn't define within from sim[1, ] here, since sim[1, ] is pre-fix!
  par <- c(witin[1:4], across[1:2], witin[7], across[3:4], witin[10], across[5:6])
  #par <- c(model2.FixParams3( data[,1:2], par ))                      # Is this line redundant since you have fixed parameters above...?
  
  a[[i]] <- simplex.MakeSimplex3(data[,1:2],par,model2.ChiSq,0.1) 
}


# Make simplex for the across-beat parameters (now new_beat can be fed in):
sim <- simplex.MakeSimplex2(data=ppg, param = par, f = model2.ChiSq3, inScale = 0.1, beat_vector = beat_vector, beat = new_beat, renal_param = renal_param)

# Make matrix again:
sim <- make_matrix(sim, a)

#######################################################################################

##################################  RUN SIMPLEX 2  ####################################

sim <- simplex.Run2(data = ppg, simplexParam = sim, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, renal_param = renal_param, run = "run 2")

# Check result:     
#model2.ChiSq3(data = ppg, params = NULL, beats = beat_vector, beat = beat, a = sim[1, ], plot = TRUE, renal_param = renal_param)

# Extract across and within beat parameters for use in next make.simplex:
across <- sim[1, 1:6]
within <- list()
for(i in 1:beats_in){
  temp <- rep(0, 12)
  temp[c(1:4, 7, 10)] <-  sim[1, ][((i*6)+1):((i*6)+6)]  
  within[[i]] <- temp
}

##################################  FIX PARAMS 2  #####################################

# Fix Params:
fixed <- list()
for(i in 1:beats_in){
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  fixed[[i]] <- model2.FixParams3(data, params = as.numeric(within[[i]]), across_beat_params = across)
}


# Update beat with fixed values:
new_beat <- data.frame(matrix(0, ncol = 12, nrow = beats_in))
for(i in 1:beats_in){
  new_beat[i, ] <- fixed[[i]]
}
# Across / within reassigned (having now fixed):
across <- as.numeric(new_beat[1, c(5:6, 8:9, 11:12)])
within <- list()
for(i in 1:beats_in){
  temp <- rep(0, 12)
  temp[c(1:4, 7, 10)] <- as.numeric(new_beat[i, ][c(1:4, 7, 10)])
  within[[i]] <- temp
}
# Add the first four columns from beat:
for(i in 1:4){
  new_beat <- cbind(beat[, i], new_beat)
}

#######################################################################################

# Check result:     
#model2.ChiSq3(data = ppg, params = as.numeric(new_beat[1, 5:16]), beats = beat_vector, beat = new_beat, a = NULL, plot = TRUE, renal_param = renal_param)

################################  Refine Parameters Again ######################################

##################################  MAKE SIMPLEX 1  ###################################

# Make simplex for the within-beat parameters:
a <- list()
for(i in 1:beats_in){                
  
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  
  witin <- within[[i]]                                            
  par <- c(witin[1:4], across[1:2], witin[7], across[3:4], witin[10], across[5:6])
  
  
  a[[i]] <- simplex.MakeSimplex3(data[,1:2],par,model2.ChiSq,0.1) 
}

# Make simplex for the across-beat parameters (now new_beat can be fed in):
sim <- simplex.MakeSimplex2(data=ppg, param = par, f = model2.ChiSq3, inScale = 0.1, beat_vector = beat_vector, beat = new_beat, renal_param = renal_param)

# Make matrix again:
sim <- make_matrix(sim, a)

#######################################################################################

##################################  RUN SIMPLEX 1  ####################################

sim <- simplex.Run2(data = ppg, simplexParam = sim, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, renal_param = renal_param, run = "run 3")

# Check result:     
#model2.ChiSq3(data = ppg, params = NULL, beats = beat_vector, beat = beat, a = sim[1, ], plot = TRUE, renal_param = renal_param)

# Extract across and within beat parameters for use in next make.simplex:
across <- sim[1, 1:6]
within <- list()
for(i in 1:beats_in){
  temp <- rep(0, 12)
  temp[c(1:4, 7, 10)] <-  sim[1, ][((i*6)+1):((i*6)+6)]  
  within[[i]] <- temp
}

##################################  MAKE SIMPLEX 2  ###################################

# Make simplex for the within-beat parameters:
a <- list()
for(i in 1:beats_in){                
  
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  
  witin <- within[[i]]                                               # You shouldn't define within from sim[1, ] here, since sim[1, ] is pre-fix!
  par <- c(witin[1:4], across[1:2], witin[7], across[3:4], witin[10], across[5:6])
  #par <- c(model2.FixParams3( data[,1:2], par ))                      # Is this line redundant since you have fixed parameters above...?
  
  a[[i]] <- simplex.MakeSimplex3(data[,1:2],par,model2.ChiSq,0.1) 
}


# Make simplex for the across-beat parameters (now new_beat can be fed in):
sim <- simplex.MakeSimplex2(data=ppg, param = par, f = model2.ChiSq3, inScale = 0.1, beat_vector = beat_vector, beat = new_beat, renal_param = renal_param)

# Make matrix again:
sim <- make_matrix(sim, a)

#######################################################################################

##################################  RUN SIMPLEX 2  ####################################

sim <- simplex.Run2(data = ppg, simplexParam = sim, f = model2.ChiSq3, optional=NULL, beat_vector = beat_vector, renal_param = renal_param, run = "run 4")

# Fix Params:
fixed <- list()
for(i in 1:beats_in){
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  fixed[[i]] <- model2.FixParams3(data, params = as.numeric(within[[i]]), across_beat_params = across)
}

# Check result:     
fit_check[[k]] <- model2.ChiSq4(data = ppg, params = NULL, beats = beat_vector, beat = beat, a = sim[1, ], plot = FALSE, renal_param = renal_param)


# Update beat with fixed values:
new_beat <- data.frame(matrix(0, ncol = 12, nrow = beats_in))
for(i in 1:beats_in){
  new_beat[i, ] <- fixed[[i]]
}
# Across / within reassigned (having now fixed):
across <- as.numeric(new_beat[1, c(5:6, 8:9, 11:12)])
within <- list()
for(i in 1:beats_in){
  temp <- rep(0, 12)
  temp[c(1:4, 7, 10)] <- as.numeric(new_beat[i, ][c(1:4, 7, 10)])
  within[[i]] <- temp
}
# Add the first four columns from beat:
for(i in 1:4){
  new_beat <- cbind(beat[, i], new_beat)
}

# Output the optimal simplex values to beat:
beat2 <- new_beat       
colnames(beat2) <- colnames(beat)
beat2 <- beat2[, -c(1:4)]

# Try plotting:
for(i in 1:beats_in){
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  yPrev <- ppg[seg[1]-1,2]
  xPrev <- ppg[seg[1]-1, 1]
  xNext <- ppg[seg[3], 1]
  rm(seg)
  temp<-model2.Rebuild2(data, yPrev, as.double(beat2[i,]),TRUE)    
  plot(data[, 1], data[, 2], ylim = c(-400, 2500), main = c("batch", k, ", wave", i))   # ylim = c(76, 86)
  lines(data[,1],temp)
  # Plot baselines:
  lines(c(xPrev, (beat2[i, 3]  + (1*beat2[i, 6]))), rep(beat2[i, 1], 2))   #0.5*beat2 for both these lines...
  lines(c((beat2[i, 3]  + (1*beat2[i, 6])), xNext), rep(beat2[i, 2], 2))
  # Plot systolic:
  # Always need 12 elements, but set amplitude and width to 0 for the peaks that aren't being used. 
  par <- as.double(beat2[i,])
  par[c(7:8, 10:11)] <- 0
  temp<-model2.Rebuild2(data,yPrev,par,TRUE)
  lines(data[,1],temp, col = "red")
  # Plot diastolic:
  par <- as.double(beat2[i,])
  par[c(4:5, 10:11)] <- 0
  temp<-model2.Rebuild2(data,yPrev,par,TRUE)      
  lines(data[,1],temp, col = "blue")
  # Plot renal:
  par <- as.double(beat2[i,])
  par[c(4:5, 7:8)] <- 0
  temp<-model2.Rebuild2(data,yPrev,par,TRUE)
  lines(data[,1],temp, col = "green")
}
}


batch <- c(fit_check[[1]][[1]], fit_check[[2]][[1]], fit_check[[3]][[1]], fit_check[[4]][[1]], fit_check[[5]][[1]], 
              fit_check[[6]][[1]], fit_check[[7]][[1]], fit_check[[8]][[1]], fit_check[[9]][[1]], fit_check[[10]][[1]])
which(batch > 500000)
mean(batch)
sd(batch)

waves <- as.numeric(c(fit_check[[1]][[2]], fit_check[[2]][[2]], fit_check[[3]][[2]], fit_check[[4]][[2]], fit_check[[5]][[2]], 
           fit_check[[6]][[2]], fit_check[[7]][[2]], fit_check[[8]][[2]], fit_check[[9]][[2]], fit_check[[10]][[2]]))
which(waves > 50000)
mean(waves)
sd(waves)

###################### New Functions: ########################


# New version of ChiSq that takes the 6-6 parameter setup:

# Description of function:

# Inputs possible:
# Data will be the ppg data
# Params - at the moment this is just being used to feed in the across_beat_params - perhaps can simplify this?
# Beats is used to feed in the beat_vector - which contains where to find the segments for each wave (and an initial value indicating how many beats to read in)
# Optional - not doing anything...? vestigeal from last ChiSq
# Beat - takes the beat dataframe or the new_beat data_frame - the info is used to retrieve input parameters for individual beats if a is null
# A - inputs a 66 parameter row for use in run.simplex
# Plot - if you want to plot the 10 beats and their fits    

# Structure:
# The parameters to be applied to all beats (across beats) will be assigned outside of the following for loop. 
# The parameters to be applied to individual beats (within beats) will be assigned inside of the for loop.
# ChiSq is calculated for each beat and the values then summated to give a value representative of all beats fitted. 

# Relation to main script:
# When parameters are first estimated they are entered into the beat dataframe. Hence the function takes beat as an input. 
# When parameters are inputted by the simplex, they are in the form of a 66*67 matrix. Hence the function can also take a as an input.
# After the simplex has run, the final parameters are reorganised into Beat2. Beat2 can be inputted through the beat argument also. 


model2.ChiSq3 <- function(data, params,debug=FALSE, beats, optional = NULL, beat, a = NULL, plot = FALSE, renal_param){  
  
  # Across-beat parameter extraction:
  if(!is.null(a)){                                            # If a 66 parameter vector has been supplied, extract the first 6 
    across_beat_params <- a[1:6]
  }else{                                                      # If not, take them from the params input
    par <- params
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]
  }
  
  # Calculation of ChiSq for all beats:
  beat_fit <- list()
  for(i in 1:beats[[1]]){                                          # The number of beats is determined by the first object of beats
    
    # Within-beat parameter extraction:
    if(!is.null(a)){                                               # If a 66 parameter vector has been supplied, take those values
      par2 <- a[((i*6)+1):((i*6)+6)]    
      par2 <- c(par2[1:4], 0, 0, par2[5], 0, 0, par2[6], 0, 0)
    }else{                                                         # If not, take the values of beat. 
      par2 <- as.numeric(beat[i, 5:16])
    }
    
    # Extract individual beat data:  
    seg <- c(beats[[2]][i],0,beats[[3]][i])
    dat <- model2.GetSegment(data,seg)
    rm(seg)
    
    # Extract systolic and diastolic parameters:
    sys <- par2[3]
    dias <- across_beat_params[2] + par2[3]
    start <- which(abs(dat[, 1]-sys) == min(abs(dat[, 1] - sys)))
    end <- which(abs(dat[, 1]-dias) == min(abs(dat[, 1] - dias)))
    
    # Find W:
    sfunction <- splinefun(1:nrow(dat), dat[, 2], method = "natural")
    deriv1 <- sfunction(seq(1, nrow(dat)), deriv = 1)
    w <- which.max(deriv1)
    
    # Fix parameters and calculate penalty:
    temp <- model2.FIX_PAR3(data = dat, params = par2, across_beat_params = across_beat_params, renal_param = renal_param)  
    penalty <- temp[1]
    fixedPar <- temp[2:length(temp)]     
    rm(temp)
    
    # Calculate fit and residue:
    fit <- model2.Rebuild2(dat,dat[1,2],params = fixedPar)    
    residue <- dat[ ,2] - fit

     # 1. Weighted region is S -> D (without slope)
     residue[start:end] <-  residue[start:end]*2
    
     # 2. Weighted region is W -> D (without slope)
     #residue[w:end] <-  residue[w:end]*2
    
     # 3. Weighted region is S -> D (with slope)
     #residue[start:end] <-  residue[start:end]*2
     #if(length(residue) > end){
     #  tail <- (end+1):length(residue)
     #  for(j in 1:length(tail)){
     #    wgt <- 2 - (0.1*j)
     #    if(wgt < 1){wgt <- 1}
     #    residue[tail[j]] <- residue[tail[j]]*wgt
     #  }
     #}
    
     # 4. Weighted region is W -> D (with slope)
     #residue[w:end] <-  residue[w:end]*2
     #if(length(residue) > end){
     #   tail <- (end+1):length(residue)
     #  for(j in 1:length(tail)){
     #     wgt <- 2 - (0.1*j)
     #   if(wgt < 1){wgt <- 1}
     #   residue[tail[j]] <- residue[tail[j]]*wgt
     #   }
     #}
    
    # Calculate Reduced Chi-Square for the beat:
    nData <- nrow(dat)    
    nPar <- length(par2)
    beat_fit[[i]] <- (sum(residue*residue) / (nData-nPar)) + as.numeric(penalty)
    
    if(plot == TRUE){
      plot(dat,  ylim = c(-150, 1600))      #ylim = c(76, 86) for bioradio data, ylim = c(-150, 1600) for ISO
      lines(dat[, 1], fit)
      #lines(dat[, 1], residue + dat[1, 2])
    }
  }
  
  # Summate individual beat ChiSq values:
  temp <- c()
  for(i in 1:length(beat_fit)){
    temp[i] <- beat_fit[[i]][1]
  }
  ts_fit <- sum(temp)
  return(ts_fit)
}

model2.ChiSq4 <- function(data, params,debug=FALSE, beats, optional = NULL, beat, a = NULL, plot = FALSE, renal_param){  
  
  # Across-beat parameter extraction:
  if(!is.null(a)){                                            # If a 66 parameter vector has been supplied, extract the first 6 
    across_beat_params <- a[1:6]
  }else{                                                      # If not, take them from the params input
    par <- params
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]
  }
  
  # Calculation of ChiSq for all beats:
  beat_fit <- list()
  for(i in 1:beats[[1]]){                                          # The number of beats is determined by the first object of beats
    
    # Within-beat parameter extraction:
    if(!is.null(a)){                                               # If a 66 parameter vector has been supplied, take those values
      par2 <- a[((i*6)+1):((i*6)+6)]    
      par2 <- c(par2[1:4], 0, 0, par2[5], 0, 0, par2[6], 0, 0)
    }else{                                                         # If not, take the values of beat. 
      par2 <- as.numeric(beat[i, 5:16])
    }
    
    # Extract individual beat data:  
    seg <- c(beats[[2]][i],0,beats[[3]][i])
    dat <- model2.GetSegment(data,seg)
    rm(seg)
    
    # Truncation of the data:
    #dat <- dat[-c(((nrow(dat)) - 10):nrow(dat)), ]
    
    # Extract systolic and diastolic parameters:
    sys <- par2[3]
    dias <- across_beat_params[2] + par2[3]
    start <- which(abs(dat[, 1]-sys) == min(abs(dat[, 1] - sys)))
    end <- which(abs(dat[, 1]-dias) == min(abs(dat[, 1] - dias)))
  
    # Fix parameters and calculate penalty:
    temp <- model2.FIX_PAR3(data = dat, params = par2, across_beat_params = across_beat_params, renal_param = renal_param)  
    penalty <- temp[1]
    fixedPar <- temp[2:length(temp)]     
    rm(temp)
    
    # Calculate fit and residue:
    fit <- model2.Rebuild2(dat,dat[1,2],params = fixedPar)    
    residue <- dat[ ,2] - fit
    # Multiply residue between systolic and diastolic peaks by 2
    residue[start:end] <-  residue[start:end]*2
    
    # Calculate Reduced Chi-Square for the beat:
    nData <- nrow(dat)    
    nPar <- length(par2)
    beat_fit[[i]] <- (sum(residue*residue) / (nData-nPar)) + as.numeric(penalty)
    
    if(plot == TRUE){
      plot(dat,  ylim = c(-150, 1600))      #ylim = c(76, 86) for bioradio data, ylim = c(-150, 1600) for ISO
      lines(dat[, 1], fit)
      #lines(dat[, 1], residue + dat[1, 2])
    }
  }
  
  # Summate individual beat ChiSq values:
  temp <- c()
  for(i in 1:length(beat_fit)){
    temp[i] <- beat_fit[[i]][1]
  }
  ts_fit <- sum(temp)
  
  fit <- list(ts_fit, beat_fit)
  
  return(fit)
}


# This version of model2.rebuild is only different in that it uses model2.Excess.Inv2 instead of model2.Excess.Inv
model2.Rebuild2 <- function(xy,offset,params,invert=TRUE){     
  result <- 1:nrow(xy) * 0.0
  # Creating the excess:
  result <- result + model2.Peak(xy[,1],params[3:5])     # Systolic parameters 
  if (length(params)>=8){
    result <- result + model2.Peak(xy[,1],params[6:8]+c(params[3],0,0))   # Diastolic parameters 
  }
  if (length(params)>=11){
    result <- result + model2.Peak(xy[,1],params[9:11]+c(params[3],0,0))  # Renal parameters
  }
  # Adding decay (config.rate + baseline parameters):
  if (invert){
    result <- model2.Excess.Inv2(xy[,1],result,offset,params[1],params[2],params[3]+1*params[6], config.rate = params[12])   # Used to be params[3]+0.5*params[6]
  }
  return(as.double(result))
}


# THIS VERSION OF MODEL2.REBUILD3 IS DESIGNED TO ONLY APPLY THE EXPONENTIAL DECAY TO THE DIASTOLIC PORTION
model2.Rebuild3 <- function(xy,offset,params,invert=TRUE){     
  result <- 1:nrow(xy) * 0.0
  # Creating the excess:
  result <- result + model2.Peak(xy[,1],params[3:5])     # Systolic parameters 
  if (length(params)>=8){
    result <- result + model2.Peak(xy[,1],params[6:8]+c(params[3],0,0))   # Diastolic parameters 
  }
  if (length(params)>=11){
    result <- result + model2.Peak(xy[,1],params[9:11]+c(params[3],0,0))  # Renal parameters
  }
  
  # Isolating individual waves:
  sys <- model2.Peak(xy[,1],params[3:5]) 
  R1 <-  model2.Peak(xy[,1],params[9:11]+c(params[3],0,0))
  R2 <- model2.Peak(xy[,1],params[6:8]+c(params[3],0,0))
  
  # Adding decay (config.rate + baseline parameters) to the systolic and diastolic peaks only:
  if (invert){
    result <- model2.Excess.Inv3(xy[,1],result,offset,params[1],params[2],params[3]+0.5*params[6], config.rate = params[12], sys, R1, R2)  
  }
  return(as.double(result))
}

# This version of model2.Excess.Inv is the one that corresponds to model2.Rebuild3
model2.Excess.Inv3 <- function(time,excess,offset,baselineStart,baselineEnd,timeBase,config.rate, sys, R1, R2){   
  nX <- length(excess)
  if (nX == 0){
    print("Help")
  }
  result <- 1:nX * 0.0
  baseline <- time * 0 + baselineStart
  baseline[which(time > timeBase)] = baselineEnd
  
  # If excess has NAs it will interrupt the reconstruction, remove them:
  temp <- which(is.nan(excess))
  if(length(temp) > 0){
    for(i in 1:length(temp)){
      excess[temp][i] <- 0 
    }
  }

  # Alternative way of adding the decay:
  # 1. Make systolic wave + decay 
  excess_1 <- sys
  #plot(excess_1)
  result_1 <- c()
  result_1[1] <- excess_1[1] + (baselineStart + config.rate*(offset-baselineStart)) 
  for (j in 2:length(excess_1)){  
    result_1[j] = excess_1[j] + (baseline[j] + config.rate*(result_1[j-1]-baseline[j]))  
  }
  #plot(result_1)
  
  # 2. Make diastolic wave + decay
  excess_2 <- R2
  #plot(excess_2)
  result_3 <- c()
  result_3[1] <- excess_2[1] + (baselineStart + config.rate*(offset-baselineStart)) 
  for (j in 2:length(excess_2)){  
    result_3[j] = excess_2[j] + (baseline[j] + config.rate*(result_3[j-1]-baseline[j]))  
  }
  #plot(result_3)
  
  # 3. Make renal wave
  excess_3 <- R1
  #plot(excess_3)
    
  # 3. Add 1 and 2. 
  result_2 <- result_1 + result_3 + excess_3
  #plot(result_2)
  
  #plot(data)
  #lines(data[, 1], result_1)
  
  return(result_2)
}


# This version of model2.Excess.Inv is only different in that is takes config.rate as an additional argument
model2.Excess.Inv2 <- function(time,excess,offset,baselineStart,baselineEnd,timeBase,config.rate){   
  nX <- length(excess)
  if (nX == 0){
    print("Help")
  }
  result <- 1:nX * 0.0
  baseline <- time * 0 + baselineStart
  baseline[which(time > timeBase)] = baselineEnd
  
  # If excess has NAs it will interrupt the reconstruction, remove them:
  temp <- which(is.nan(excess))
  if(length(temp) > 0){
    for(i in 1:length(temp)){
      excess[temp][i] <- 0 
    }
  }
  
  # Adding the decay element to the excess (one value at a time):
  result[1] = excess[1] + (baselineStart + config.rate*(offset-baselineStart)) 
  for (j in 2:nX){  
    result[j] = excess[j] + (baseline[j] + config.rate*(result[j-1]-baseline[j]))  
  }
  return(result)
}



# Different to previous FIX_PAR versions, across_beat_params is now an additional argument. 
# Conditional statements have been altered to reflect the new max number of parameters per beat (12)
# Additional constraints have been placed on waves to encourage proper fits. 

model2.FIX_PAR3 <- function(data,params,across_beat_params, debug=FALSE, renal_param){
  
  par <- params
  
  nData <- nrow(data)
  nPar <- length(par)
  
  # Transcribe parameters
  nBase <- 1
  baseline <- c( par[1], par[1] )
  if (nPar == 6 | nPar == 9 | nPar == 12){   
    baseline[2] = par[2]
    nBase <- 2
  }
  
  t <- c( par[nBase + 1], across_beat_params[2], across_beat_params[4])           # time (systolic = within, diastolic / renal = across)
  h <- c( par[nBase + 2], 0, 0 )                                                  # height (systolic = within, diastolic / renal default to 0 unless peaks supplied (see below))
  w <- c( across_beat_params[1], across_beat_params[3], across_beat_params[5])    # width (systolic / diastolic / renal = across)
  hasPeak <- c( TRUE, FALSE, FALSE )
  
  # Calculate penalty only if a peak has been supplied       
  
  if (nPar >= nBase + 7){                                    # Assign diastolic values if peak present
    hasPeak[2] = TRUE
    t[2] <- across_beat_params[2]
    h[2] <- par[nBase + 5]
    w[2] <- across_beat_params[3]
  }
  
  if (nPar >= nBase + 10){                                   # Assign renal values if peak present  
    hasPeak[3] = TRUE
    t[3] <- across_beat_params[4]
    h[3] <- par[nBase + 8]
    w[3] <- across_beat_params[5]
  }
  
  # Clamp and/or penalize parameters
  penalty <- 0
  
  # TMIN AND TMAX CAN STAY (mark beginning and end of the beat segment)
  tMin <- data[1,1]
  tMax <- data[nData,1]
  
  META_BASELINE_SHIFT <- 1.0    # penalty for how big the gap is between the two baselines
  META_MIN_PEAK_DELAY <- 0.1    # peaks cannot be following one another by less than 0.1ms
  
  p <- 1:12*0    # One penalty value for each parameter (was 1:11)
  
  # Fix height and width for each of the three waves (1:3)      
  for (i in 1:3){
    if (h[i] < 0){                                           # Heights should not be negative
      penalty <- penalty + h[i]*h[i] 
      p[3*i+1] <- h[i]*h[i]            # inside the subset was 3*i-1 - surely this would correspond to the wrong parameters...?       
      h[i] <- 0
    }
    
    if (w[i] < 0.05 | w[i] > 0.5){    
      fixed <- max( 0.05, min( w[i], 0.5 ) )                 # Width should be > 0.05 and < 0.5
      diff <- fixed - w[i]
      penalty <- penalty + diff*diff
      p[3*i+2] <- diff*diff            
      w[i] <- fixed
    }
    
    ############# My own Insertion 30/1/21 ####################
    
    # If renal peak starts before systolic peak, penalize
    #if(i==3){
    #  if(((t[1] + t[3]) - (w[3]/2)) < data[which.max(data[, 2]), 1]){   
    # can divide w[3] by 2 for a different (possibly better?) cutoff
    #    penalty <- penalty + 10000000
    #  }
    #}
    
    if(i==3){                                              # Renal width should not be greater than 0.25...
      if(w[3] > 0.25){
        diff <- 0.25 - w[3]
        penalty <- penalty + diff*diff
        w[3] <- 0.25
      }
    }
    
    if(i==3){                                              # Renal width should not be less than 0.1...
      if(w[3] < 0.1){
        diff <- 0.1 - w[3]
        penalty <- penalty + diff*diff
        w[3] <- 0.1
      }
    }
    
    if(i==2){                                             # Diastolic width should not be greater than 0.45
      if(w[2] > 0.45){
        diff <- 0.45 - w[2]
        penalty <- penalty + diff*diff
        w[2] <- 0.45
      }
    }
    
    if(i==3){                                              # Renal peak should be penalized in a graded way as its amplitude increases
      if(h[3] > 25){
        penalty <- penalty + 1000
      }
      if(h[3] > 50){
        penalty <- penalty + 5000
      }
      if(h[3] > 75){
        penalty <- penalty + 10000
      }
    }
    
    if(i==3){                                              # Renal peak timing should be similar to initial estimation
      if(t[3] > (renal_param + 0.02) ){
        penalty <- penalty + 5000
      }
    }
    
    if(i==3){                                              # Renal peak timing should be similar to initial estimation
      if(t[3] < (renal_param - 0.02) ){
        penalty <- penalty + 5000
      }
    }
    
    # SWidth:
    #if(i==1){
    #    if(w[1] > 0.27){
    #      penalty <- penalty + 150000
    #    }
    #}
    
    # S amplitude:
    if(i==1){
      max.amp <- data[which.min(abs(data[, 1] - t[1])), 2]
      if(h[1] > max.amp + 50){
        penalty <- penalty + 100000
      }
      if(h[1] > max.amp - 50){
        penalty <- penalty + 100000
      }
      if(h[1] > max.amp + 100){
        penalty <- penalty + 200000
      }
      if(h[1] > max.amp - 100){
        penalty <- penalty + 200000
      }
      
      
    }
    
    ########################################################### 
  }
  
  
  # Fix time
  
  # Systolic:
  fixed <- max( tMin, min( t[1], tMin + 1, tMax ) )      # Making sure S peak sits between min and max times of the segment
  if (debug){
    print(paste("time S: ",tMin," < ",t[1]," < min( ",tMin+1,",",tMax," )"))
  }
  if (t[1] != fixed){
    diff <- fixed - t[1]
    penalty <- penalty + diff*diff   
    p[1] <- diff*diff
    t[1] <- fixed
  }
  
  # Second S time Fix (likely to only work when systolic peak is actually the maximum of the data i.e not on class 3 waveforms:
  
  y. <- data[which.max(data[, 2]), 1]
  if(t[1] > y. + 0.04){
    penalty <- penalty + 10000
  }
  if(t[1] < y. - 0.04){
    penalty <- penalty + 10000
  }
  
  # Diastolic:
  fixed <- max( 2 * META_MIN_PEAK_DELAY, min( t[2], tMax - tMin + 0.4 * w[2] ) )   # This stops diastolic time being < 0.2
  if (debug){
    print(paste("time D: ",2 * META_MIN_PEAK_DELAY," < ",t[2]," < ",tMax - tMin + 0.4 * w[2]," )"))
  }
  if (t[2] != fixed){
    # Two peak delays between S and D
    diff <- fixed - t[2]
    if (hasPeak[2]){
      penalty <- penalty + diff*diff  
      p[4] <- diff*diff
    }
    t[2] <- fixed
  }
  
  # Renal:
  fixed <- max( META_MIN_PEAK_DELAY, min( t[3], t[2] - META_MIN_PEAK_DELAY ) )   # Stops renal peak being < 0.1 after systolic or < 0.1 before diastolic
  if (debug){
    print(paste("time R: ",META_MIN_PEAK_DELAY," < ",t[3]," < ",t[2] - META_MIN_PEAK_DELAY," )"))
  }
  if (t[3] != fixed){
    diff <- fixed - t[3]
    if (hasPeak[3]){
      penalty <- penalty + diff*diff 
      p[7] <- diff*diff
    }
    t[3] <- renal_param
  }
  
  # Don't clamp baseline shift, but penalize large shifts   (there will always be two baselines, but they will be the same if them being different results in no significantly better fit)
  diff <- abs(baseline[1] - baseline[2])    # + 100?
  #penalty <- penalty + META_BASELINE_SHIFT*diff  #*diff     
  
  ### Insertion 9/2/21 #######
  
  if(across_beat_params[6] > 0.95){
    #diff <- 0.85 - across_beat_params[6]    # diff here won't help much since you're value will be something like 0.1
    penalty <- penalty + 1000
  }
  
  ############################
  
  fixedPar <- c( baseline, t[1], h[1], w[1], t[2], h[2], w[2], t[3], h[3], w[3], across_beat_params[6])
  
  if (debug){
    print(p)
  }
  
  return( c( penalty, fixedPar ) )
}


model2.FixParams3 <- function(data,params, across_beat_params = NULL, debug=FALSE, rp = renal_param){
  
  # If across_beat_params have not been provided, extract them from params
  if(is.null(across_beat_params)){    
    across_beat_params <- params[c(5, 6, 8, 9, 11, 12)]
  }
  
  temp <- model2.FIX_PAR3(data, params, across_beat_params, debug, rp)  
  return( temp[2:length(temp)] )     # first value of temp is penalty
} 




###################################### NEW SIMPLEX FUNCTIONS #########################################

# Make Simplex 2 (for across beat parameters only):

simplex.MakeSimplex2 <- function(data,param,f,inScale,directions=NULL,inTol=-1,optional=NULL,debug=FALSE, beat_vector = beat_vector, beat = beat, renal_param = renal_param){
  if(debug){print("MakeSimplex -- debug")}                        
  nPar <- length(param)
  nScale <- length(inScale)
  if (nScale == 0)
  {
    scale <- 1:nPar * 0 + 1
  } else if (nScale == 1){
    scale <- 1:nPar * 0 + inScale
  } else if (length(inScale) == nPar){
    scale <- inScale
  } else {
    #print("Invalid scale vector length")
    return("Error: Invalid scale vector length")
  }
  if (length(inTol) == 1 & inTol > 0){
    tol <- inTol[1]
  } else {
    tol <- min(1,f(data,param, beats = beat_vector, beat = beat, renal_param = renal_param))  
  }
  
  
  chiSq <- 1:(nPar+1) * 0.0
  chiSq[1] <- f(data,param,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param)   # ChiSq[1 is the fit when no parameters are changed... 
  if (debug){ print(paste("Root chi-squared:",chiSq[1]))}
  
  result <- matrix(nrow=nPar+1,ncol=nPar)
  result[1,] <- as.double(param)
  
  useDirections = !is.null(directions)
  if (useDirections){ useDirections <- nrow(directions) == nPar & ncol(directions) == nPar }
  
  
  for (i in c(5, 6, 8, 9, 11, 12)){   # INSTEAD OF 1:nPar
    if (debug){ print(paste("Parameter",i)) }
    
    tParam <- param
    
    # Pick a direction
    delta <- 1:nPar * 0
    if (useDirections){
      delta <- scale[i] * directions[i,]
    } else {
      delta[i] <- scale[i]
    }
    
    tParam <- param - delta                                   # This part tries tweaking each parameter up or down, 
    chiSqMinus <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param)            # tParam = test parameter. 
    tParam <- param + delta                                   # The chisquare (goodness of fit) is calculated for each direction, 
    chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param)            # and presumably the direction with the smaller value is chosen. 
    
    if (debug){
      print("Select direction:")
      print(paste("chi^2(",param[i] - delta[i],") =",chiSqMinus))
      print(paste("chi^2(",param[i],") =",chiSq[1]))
      print(paste("chi^2(",param[i] + delta[i],") =",chiSq[i+1]))
      print("---")
    }
    
    if (chiSqMinus < chiSq[i+1]){      # If going down by delta is better than going up by delta, 
      delta <- -delta                  # then replace chiSq[i+1] with the lower score (ChiSqMinus)
      tParam <- param + delta
      chiSq[i+1] <- chiSqMinus
    }
    
    iKill <- 10    
    
    if (chiSq[i+1] < chiSq[1]){         # If the new fit is better than the old fit, continue to go in the direction that improved the fit
      if (debug){ print("Extending as best point") }
      while (chiSq[i+1] < chiSq[1] + tol){                 # Chisquare keeps getting iterated here (for 10 iterations)
        delta <- 2*delta    # 2* was too much here... WHY DOES IT SEEM LIKE SAMPLITUDE IS COMING DOWN???
        tParam <- param + delta
        oldScore <- chiSq[i+1]          # The current best fit gets called 'old score'
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param)   # The new fit is now designated ChiSq[i+1]
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] > oldScore){   # Check if the new fit is worse than current fit
          tParam <- param + 0.5*delta   # If so, make delta what it was one iteration previous (undoing the *2)
          chiSq[i+1] <- oldScore        # and redesignate old score to chiSq[i+1]
          break
        }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        iKill <- iKill - 1             # If the new fit is not worse than the current fit, knock of one on the interation count,
        if (iKill < 0){                # and keep iterating until either the next fit is worse (and break is called), or iKill ends (and break is also called)
          break
        }
      }
    } else if (chiSq[i+1] < chiSq[1] + tol){    # If the new fit is not better than the old fit, is it at least better than the old fit + tol?
      if (debug){ print("Extending below tolerance") }
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }          # repeated code up until here... 
        if (chiSq[i+1] - oldScore < oldScore - chiSq[1]){   # Presumably again if the fit was worse than before, reverse it
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        iKill <- iKill - 1
        if (iKill < 0){
          print("Failed to construct simplex")
          return(paste("Error: param[",i,"]",sep=""))
        }
      }
    } else {
      if (debug){ print("Shrinking above tolerance") }
      while (chiSq[i+1] > chiSq[1] + tol){              # If the new fit is much worse than the original, reduce the size of delta
        delta <- 0.5*delta
        tParam <- param + delta
        lastChiSq <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        if (iKill < 0 & (chiSq[i+1]-chiSq[1]) > 0.75 * (lastChiSq-chiSq[1])){
          if(i == 9){
            tParam[9] <- renal_param  # Ignore renal times that can't optomize
            next
          } 
          print(c("simplex constructed as per original parameter"))
          next
          #print("Failed to construct simplex")
          #return(paste("Error: param[",i,"]",sep=""))   
        }
        iKill <- iKill - 1
      }
      tParam <- param + 0.5 * delta
    }
    
    if(debug){ print(paste("Param[",i,"] =",tParam[i]))}
    result[i+1,] = as.double(tParam) 
    
  }
  
  # PARAMETER TWEAKING ENDS HERE... 
  
  if (debug){ print("/MakeSimplex") }
  return(result)
}



# Make simplex 3 (for within-beat parameters only)
simplex.MakeSimplex3 <- function(data,param,f,inScale,directions=NULL,inTol=-1,optional=NULL,debug=FALSE){
  if(debug){print("MakeSimplex -- debug")}
  nPar <- length(param)
  nScale <- length(inScale)
  if (nScale == 0)
  {
    scale <- 1:nPar * 0 + 1
  } else if (nScale == 1){
    scale <- 1:nPar * 0 + inScale
  } else if (length(inScale) == nPar){
    scale <- inScale
  } else {
    #print("Invalid scale vector length")
    return("Error: Invalid scale vector length")
  }
  if (length(inTol) == 1 & inTol > 0){
    tol <- inTol[1]
  } else {
    tol <- min(1,f(data,param))    # here is the issue, in model2.chisquare
  }
  
  chiSq <- 1:(nPar+1) * 0.0
  chiSq[1] <- f(data,param)
  if (debug){ print(paste("Root chi-squared:",chiSq[1]))}
  
  result <- matrix(nrow=nPar+1,ncol=nPar)
  result[1,] <- as.double(param)
  
  useDirections = !is.null(directions)
  if (useDirections){ useDirections <- nrow(directions) == nPar & ncol(directions) == nPar }
  
  for (i in c(1:4, 7, 10)){    # within-beat parameters only
    if (debug){ print(paste("Parameter",i)) }
    tParam <- param
    
    # Pick a direction
    delta <- 1:nPar * 0
    if (useDirections){
      delta <- scale[i] * directions[i,]
    } else {
      delta[i] <- scale[i]
    }
    
    if(i == 3){
      delta <- delta/4
    }
    
    tParam <- param - delta                                   # This part tries tweaking each parameter up or down, 
    chiSqMinus <- f(data,tParam)            # tParam = test parameter. 
    tParam <- param + delta                                   # The chisquare (goodness of fit) is calculated for each direction, 
    chiSq[i+1] <- f(data,tParam)            # and presumably the direction with the smaller value is chosen. 
    
    if (debug){
      print("Select direction:")
      print(paste("chi^2(",param[i] - delta[i],") =",chiSqMinus))
      print(paste("chi^2(",param[i],") =",chiSq[1]))
      print(paste("chi^2(",param[i] + delta[i],") =",chiSq[i+1]))
      print("---")
    }
    
    if (chiSqMinus < chiSq[i+1]){
      delta <- -delta
      tParam <- param + delta
      chiSq[i+1] <- chiSqMinus
    }
    
    iKill <- 10    
    
    if (chiSq[i+1] < chiSq[1]){
      if (debug){ print("Extending as best point") }
      while (chiSq[i+1] < chiSq[1] + tol){                 # Chisquare keeps getting iterated here (for 10 iterations)
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] > oldScore){
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        iKill <- iKill - 1
        if (iKill < 0){
          break
        }
      }
    } else if (chiSq[i+1] < chiSq[1] + tol){
      if (debug){ print("Extending below tolerance") }
      while (chiSq[i+1] < chiSq[1] + tol){
        delta <- 2*delta
        tParam <- param + delta
        oldScore <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        if (chiSq[i+1] - oldScore < oldScore - chiSq[1]){
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        iKill <- iKill - 1
        if (iKill < 0){
          #print("Failed to construct simplex")
          return(paste("Error: param[",i,"]",sep=""))
        }
      }
    } else {
      if (debug){ print("Shrinking above tolerance") }
      while (chiSq[i+1] > chiSq[1] + tol){
        delta <- 0.5*delta
        tParam <- param + delta
        lastChiSq <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }
        #print(paste(i,"-",delta,":",chiSq[i+1]))
        if (iKill < 0 & (chiSq[i+1]-chiSq[1]) > 0.75 * (lastChiSq-chiSq[1])){
          print(c("Failed to construct simplex within 10 iterations for parameter", i, "defaulting to inputted value"))
          #return(paste("Error: param[",i,"]",sep=""))
          tParam[i] <- param[i]
          next
        }
        iKill <- iKill - 1
      }
      tParam <- param + 0.5 * delta
    }
    
    if(debug){ print(paste("Param[",i,"] =",tParam[i]))}
    result[i+1,] = as.double(tParam)
  }
  
  if (debug){ print("/MakeSimplex") }
  return(result)
}





simplex.Run2 <- function(data = ppg,simplexParam = sim,f = model2.ChiSq3,optional=NULL, beat_vector = beat_vector, renal_param = renal_param, run = NULL){
  
  MAX_STEP <- 10000                                               # The number of steps to iterate through
  FTOL <- 1e-5                                  
  
  debugRtol <- 1:(MAX_STEP+1) * 0.0
  debugMin <- 1:(MAX_STEP+1) * 0.0
  debugMax <- 1:(MAX_STEP+1) * 0.0
  
  result <- simplexParam                         # Now feed in the 66*66 matrix
  nPar <- ncol(result)                           
  chiSq <- 0:nPar * 0.0
  for (i in 1:(nPar+1)){                                    # Find out the ChiSq value for each row from result
    chiSq[i] <- f(data, params = NULL, optional=NULL, a = result[i, ], beats = beat_vector, renal_param = renal_param)
  }
  
  for (iStep in 1:MAX_STEP){                             # This is a long for loop....
    extrema <- simplex.SortHighLow(chiSq)                # Finds the results which give the highest, 2nd highest and lowest ChiSq
    low <- extrema[1]
    nHigh <- extrema[2]
    high <- extrema[3]
    
    if(!is.null(run)){
      print(run)
    }
    print(iStep)
  
    chiSqMax <- chiSq[high]                        
    chiSqMin <- chiSq[low]                  
    
    print(chiSqMax)
    
    #print(paste("chi^2_min =",chiSqMin))
    #print(paste("argMax = ",high,"[",chiSqMax,"]",sep=""))
    
    rtol <- 2 * (chiSqMax - chiSqMin)/(chiSqMax + chiSqMin + 1e-10)   # Some measure of how much better high is from low...
    if (rtol < FTOL){
      bestParam <- result[low,]                     # Presumably if the difference in ChiSq (max vs min) is significant, 
      result[low,] <- result[1,]                    # the result that was changed to give the lowest ChiSq gets designated 'best Param'
      result[1,] <- bestParam                       # The changed parameter gets upgraded to first row (swapped with what is there currently)
      return(result)
    }
    debugRtol[iStep] <- rtol
    debugMin[iStep] <- chiSqMin
    debugMax[iStep] <- chiSqMax
    
    factor <- -1
    node <- simplex.HypoCentre(result,high)        # Hypocentre outputs all the parameters that are not the worst
    apex <- result[high,]                          # Apex must be the worst parameter
    test <- node - (apex - node)                   # This represents the flipping of the triangle; the whole parameter set is reversed in the direction away from the worst ChiSq point (literally subtracting one row from another here)
    score <- f(data, params = rep(0, 12),optional=optional, a = test, beats = beat_vector, renal_param = renal_param)
    
    if (score < chiSqMin){                          # If flipping improves the ChiSq, try extending further in the same direction
      test2 <- node - 2 * (apex - node)
      score2 <- f(data, params = rep(0, 12),optional=optional, a = test2, beats = beat_vector, renal_param = renal_param)
      if (score2 >= score){                       # If reflecting a further distance is better than reflecting alone, do that
        # Reflect
        #print(paste("Reflecting",high,": chi^2 ",chiSqMax,"->",score,sep=""))
        result[high,] <- test
        chiSq[high] <- score
      } else {
        # Reflect and grow
        #print(paste("Reflect-stretching",high,": chi^2 ",chiSqMax,"->",score2,sep=""))
        result[high,] <- test2
        chiSq[high] <- score2
      }
    } else if (score >= chiSq[nHigh]) {              # If reflecting is not beneficial, try shrinking instead of reflecting
      # Test for shrink with optional reflection
      factor <- 0.5
      if (score < chiSqMax)
      {
        factor <- -0.5
      }
      test2 <- node + factor * (apex - node)
      score2 <- f(data, params = rep(0, 12),optional=optional, a = test2, beats = beat_vector, renal_param = renal_param)
      if (score2 < chiSq[nHigh]){
        # Shrink (possibly reflecting)
        #print(paste("Shrinking",high,": chi^2 ",chiSqMax,"->",score2,sep=""))
        result[high,] <- test2
        chiSq[high] <- score2
      } else {
        # Shrink all
        for (i in 1:(nPar+1)){
          if (i != low){
            result[i,] <- 0.5 * (result[i,] + result[low,])
            chiSq[i] <- f(data, params = rep(0, 12),optional=optional, a = result[i, ], beats = beat_vector, renal_param = renal_param)
          }
        }
        #print(paste("General contraction: chi^2 ",chiSqMax,"->",max(chiSq),sep=""))
      }
    } else {
      # Reflect
      #print(paste("Reflecting*",high,": chi^2 ",chiSqMax,"->",score,sep=""))
      result[high,] <- test
      chiSq[high] <- score
    }
  }
  
  extrema <- simplex.SortHighLow(chiSq)
  low <- extrema[1]
  bestParam <- result[low,]
  result[low,] <- result[1,]
  result[1,] <- bestParam
  
  chiSqMax <- chiSq[extrema[3]]
  chiSqMin <- chiSq[low]
  rtol <- 2 * (chiSqMax - chiSqMin)/(chiSqMax + chiSqMin + 1e-10)
  debugRtol[MAX_STEP+1] <- rtol
  debugMin[MAX_STEP+1] <- chiSqMin
  debugMax[MAX_STEP+1] <- chiSqMax
  plot(debugMax,type='l')
  lines(debugMin)
  
  
  print(paste("Terminated downhill simplex after",MAX_STEP,"iterations."))
  print(paste("rtol =",rtol))
  return(result)
}




# Make Matrix (66*66):

make_matrix <- function(sim, a){
  # Save the top row of sim for replication:
  top_row_sim <- sim[1, c(5:6, 8:9, 11:12)]                      # GET TOP ROW OF SIM (across beat parameters)
  #Remove redundant rows and columns from sim:
  sim <- sim[c(6:7, 9:10, 12:13), c(5:6, 8:9, 11:12)]                             # GET EXPERIMENTAL ROWS OF SIM
  # You need the top_row of sim to replicate:
  top_row_sim <- matrix(data = top_row_sim, nrow = 6, ncol = 6, byrow = TRUE)     # REPLICATE TOP ROW OF SIM
  
  
  # Make the a values just the rows where within beat parameters are changed  
  # Save the top row of each matrix of a for replication...
  top_row <- list()
  for(i in 1:beats_in){                                                            # GET TOP ROWS OF A (within beat parameters)
    top_row[[i]] <- a[[i]][1, -c(5:6, 8:9, 11:12)]
    a[[i]] <- a[[i]][-c(1, 6:7, 9:10, 12:13), -c(5:6, 8:9, 11:12)]                 # AND GET EXPERIMENTAL ROWS OF A
  }
  # Top row needs to be replicated for each within beat row when the across-beat params are being changed...:
  for(i in 1:beats_in){
    top_row[[i]] <- matrix(data = top_row[[i]], ncol = 6, nrow = 6, byrow = TRUE)   # REPLICATE TOP ROWS OF A
  }
  
  
  # Assemble Matrix:
  
  # Bind replicate rows of top_row for each beat to sim:              # CREATING ACROSS BEAT ROWS (across-beat parameters)
  for(i in 1:beats_in){
    sim <- cbind(sim, top_row[[i]])
  }
  
  
  beat_rows <- list()                                                 # CREATING ROWS FOR EACH BEAT (within-beat parameters)
  for(i in 1:beats_in){                                             
    
    # ADD LEFT
    beat_rows[[i]] <- top_row_sim   
    if(i != 1){
      for(j in 1:(i-1)){
        beat_rows[[i]] <- cbind(beat_rows[[i]], top_row[[j]])  
      }
    }
    
    # ADD A
    beat_rows[[i]] <- cbind( beat_rows[[i]], a[[i]])
    
    # ADD RIGHT
    if(i == beats_in){
      break
    }else{
      for(j in (i+1):beats_in){
        beat_rows[[i]] <- cbind(beat_rows[[i]], top_row[[j]])
      }
    }
  }
  
  # Bind rows together:
  for(i in 1:beats_in){
    sim <- rbind(sim, beat_rows[[i]])
  }
  
  # Add the top row:
  final_top_row <- top_row_sim[1, ]
  for(i in 1:beats_in){
    final_top_row <- c(final_top_row, top_row[[i]][1, ])
  }
  sim <- rbind(final_top_row, sim)
  
  return(sim)
}
