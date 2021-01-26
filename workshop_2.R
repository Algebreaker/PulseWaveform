setwd("C:/Users/jojuf/Documents/GitHub/Tulsa Model 2/Workshop2")
source("model2.R")
source("simplex.R")

config.sample_rate = 75.0
lab.time = "time (s)"
lab.ppg = "Detrended"
lab.dt = "interval (s)"


# Load PPG

ppg <- read.csv("ppg.csv")
names(ppg)[1] <- lab.time
names(ppg)[2] <- lab.ppg
ppg$Baseline = 1:nrow(ppg) * 0
ppg$Excess   = 1:nrow(ppg) * 0
ppg$Residue  = 1:nrow(ppg) * 0

# Bookmark 1

# Load Intervals

beat <- read.csv("beat.csv")
nBeats <- nrow(beat)
names(beat)[1] <- lab.time
beat$First      = 1:nrow(beat) * 0
beat$Last       = 1:nrow(beat) * 0
beat$Baseline   = 1:nrow(beat) * 0
beat$Baseline2  = 1:nrow(beat) * 0
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

plot(ppg[,1],ppg[,2],t='l')
points(beat[,1],(1:nBeats)*0+84)


# Estimate Parameters

seg <- c(0,0,0)

for (i in 1:nBeats){
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
  
  baseline <- min(data[,2])-0.326

  par <- 1:11 * 0.
  par[1:2] = baseline

  excess <- model2.Excess( data[,2], yPrev, baseline )
  residue <- excess  

  
  # S peak
  
  peak.w <- which(data[,1] > beat[i,1]-0.2 & data[,1] < beat[i,1]+0.2)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  #
  par[4] <- max(peak.y)
  par[3] <- peak.t[which(peak.y==par[4])]
  par[5] <- 0.25
  #
  rm(peak.w,peak.t,peak.y)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[3:5])
  
  #plot(data[,1],excess)
  #lines(data[,1],residue)

  
  # D peak
  
  peak.w <- which(data[,1] > beat[i,1]+0.2 & data[,1] < beat[i,1]+0.6)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  w <- which(peak.y==max(peak.y))
  #
  par[7] <- max(peak.y)
  par[6] <- peak.t[which(peak.y==par[7])]
  par[8] <- 0.25
  #
  rm(peak.w,peak.t,peak.y,w)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[6:8])
  
  plot(data[,1],excess)
  lines(data[,1],residue)
  
  
  # N peak
  
  t <- par[3] + c(0.25,0.75) * (par[6]-par[3])
  peak.w <- which(data[,1] > t[1] & data[,1] < t[2])
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  #
  par[10] <- max(peak.y)
  par[9] <- peak.t[which(peak.y==par[10])]
  par[11] <- 0.25
  #
  rm(peak.w,peak.t,peak.y,t)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[9:11])
  
  plot(data[,1],excess)
  lines(data[,1],residue)
  

  # Store parameters
  
  w <- seg[1]:seg[3]
  ppg$Baseline[w] <- baseline
  ppg$Excess[w] <- excess
  ppg$Residue[w] <- residue
  rm(w,excess,residue,data,baseline,yPrev,nextTime,tStart)
  
  beat[i,3:4]  = c(seg[1],seg[3])
  beat[i,5:15] = par
  beat[i,10] = beat[i,10]-beat[i,7]
  beat[i,13] = beat[i,13]-beat[i,7]
  rm(par)
}
rm(seg)

write.csv(beat,file="beat_2.csv", row.names=FALSE)
# BOOKMARK 2
beat <- read.csv("beat_2.csv") 
nBeats <- nrow(beat)
i=1

# Refine Parameters

for (i in 1:nBeats){
  print(i)
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  yPrev <- ppg[max(seg[1]-1,1),2]
  rm(seg)  #plot(data[,1:2])
  
  #params <- as.numeric(beat[i,5:15])
  #test <- model2.Rebuild( data, yPrev, params, invert=FALSE)
  #test <- model2.Rebuild( data, yPrev, params, invert=TRUE)
  #val <- model2.ChiSq( data, params, debug=TRUE )
  #sim <- simplex.MakeSimplex(data[,1:2],params,model2.ChiSq,0.1)
  
  
  
  par <- as.numeric(beat[i,5:15])
  par <- model2.FixParams2( data[,1:2], par )
  
  sim <- simplex.MakeSimplex(data[,1:2],par,model2.ChiSq2,0.1)
  sim <- simplex.Run(data[,1:2],sim,model2.ChiSq2)
  par <- model2.FixParams2( data[,1:2], sim[1,] )
  
  sim <- simplex.MakeSimplex(data[,1:2],sim[1,],model2.ChiSq2,0.1) #run it all a second time, should be similar to values from first run-through
  sim <- simplex.Run(data[,1:2],sim,model2.ChiSq2)
  par <- model2.FixParams2( data[,1:2], sim[1,] )
  

  test <- model2.Rebuild( data, yPrev, sim[1,], invert=TRUE)
  plot(data[,1],data[,2])
  lines(data[,1],test)

  rm(data,sim)
  
  beat[i,5:15] = par
  
  rm(par)
}

write.csv(beat,file="beat_3.csv")
# BOOKMARK 3

## END OF WORKSHOP PREP


i<-3
seg <- c(beat[i,3],0,beat[i,4])
data <- model2.GetSegment(ppg,seg)
rm(seg)
plot(data[,1],data[,2])
temp<-model2.Rebuild(data,79.5,beat[i,5:14],TRUE)
lines(data[,1],temp)


# Refine Parameters Again

fixedParams <- c(
  median(beat$SWidth),  # S width
  median(beat$DTime),  # D timing
  median(beat$DWidth),  # D width
  median(beat$NTime), # N timing
  median(beat$NWidth)  # N width
)

for (i in 1:nBeats){
  print(i)
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg) #we could truncate the segment here
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

i = i
seg <- c(beat[i,3],0,beat[i,4])
yPrev <- ppg[seg[1]-1,2]
data <- model2.GetSegment(ppg,seg)
rm(seg)

plot(data[,1],data[,2])
temp<-model2.Rebuild(data,yPrev,as.double(beat[i,5:8]),TRUE)
lines(data[,1],temp)
temp<-model2.Rebuild(data,yPrev,as.double(beat[i,c(5:8,12:14)]),TRUE)
lines(data[,1],temp)
temp<-model2.Rebuild(data,yPrev,as.double(beat[i,5:14]),TRUE)
lines(data[,1],temp)

temp<-model2.Rebuild(data,yPrev,as.double(c(beat[i,5],beat[i,6]+beat[i,9],beat[i,10:11])),TRUE)
lines(data[,1],temp)
                     
temp<-model2.Rebuild(data,yPrev,as.double(c(beat[i,5],beat[i,6],0,0,beat[i,9],beat[i,10:11])),TRUE)
lines(data[,1],temp)


# ISO
ppg <- read.csv("AA343.csv")
ppg <- data.frame(
  time = (0:(nrow(ppg)-1)) / 40,
  ppg = ppg[,1]
)

n <- dim(ppg)[1]
n <- nrow(ppg)
vpg <- ppg[2:n,2] - ppg[1:(n-1),2]
beat <- ppg[which(vpg[1:(n-1)] < 300 & vpg[2:n] >= 300),1]

nBeat <- length(beat)
beat <- data.frame(
  beat = beat,
  dt = (1:nBeat)*0.0
)
rm(vpg)

ppg[,2] = UnDetrend(ppg,factor=0.9,offset=9)
