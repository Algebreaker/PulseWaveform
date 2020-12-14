source("model2.R") #assuming each of the three waveforms is motivated by constant pressure in the background, if pressure goes up, either systolic wave or reflectance wave, what we are left with accounting for this is the excess
source("simplex.R")


input <- "A"

config.sample_rate = 75.0
lab.time = "time (s)"
lab.ppg = "Detrended"
lab.dt = "interval (s)"


# Load PPG

ppg <- read.csv(paste(input,"_ppg.csv",sep=""))
plot(ppg)
names(ppg)[1] <- lab.time
names(ppg)[2] <- lab.ppg


# Load Intervals

beat <- read.csv(paste(input,"_intervals.csv",sep=""))
names(beat)[1] <- lab.time
names(beat)[2] <- lab.dt


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

for (i in 1:nBeats){
  beatTime <- beat[i,1]
  nextTime <- if (i < nBeats){beat[i+1,1]}else{NA}
  temp = seg[3]
  seg <- model2.FindSegment(ppg,beat[i,1],nextTime) #estimating which is the first and last bit to fit (preceding and following value of O)i=1
  if (temp > 0){
    seg[1] = temp + 1
  }
  rm(temp)
  

  data <- model2.GetSegment(ppg,seg)
  tStart <- ppg[seg[1],1]
  yPrev <- ppg[max(seg[1]-1,1),2] #yPrev the value of the pressure before the segment starts
  
  baseline <- min(data[,2])-0.326 #a constant, will be wildly different for datasets, when scale changes
  # excess <- model2.Excess(data[,2],yPrev,baseline)

  count <- nrow(data)
  excess <- 1:count * 0.0
  excess[1] = data[1,2] - (baseline + config.rate*(yPrev-baseline))
  for (j in 2:count){
    excess[j] = data[j,2] - (baseline + config.rate*(data[j-1,2]-baseline))
  }
  rm(count)
  rm(j) #plot the excess: plot(data[,1],excess), when it reaches 0 in the tail it means there is exponential decay

  
  par <- 1:10 * 0. 
  par[1] = baseline

  residue <- excess  

  
  # S peak
  
  peak.w <- which(data[,1] > beat[i,1]-0.2 & data[,1] < beat[i,1]+0.2)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  #
  par[3] <- max(peak.y)
  par[2] <- peak.t[which(peak.y==par[3])]
  par[4] <- 0.25 #width of the peak
  #
  rm(peak.w,peak.t,peak.y)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[2:4])
  
  #plot(data[,1],excess)
  #lines(data[,1],residue)

  
  # D peak
  
  peak.w <- which(data[,1] > beat[i,1]+0.2 & data[,1] < beat[i,1]+0.6)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  #
  par[6] <- max(peak.y)
  par[5] <- peak.t[which(peak.y==par[6])]
  par[7] <- 0.25
  #
  rm(peak.w,peak.t,peak.y)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[5:7])
  
  plot(data[,1],excess)
  lines(data[,1],residue)
  
  
  # N peak
  
  t <- par[2] + c(0.25,0.75) * (par[5]-par[2])
  peak.w <- which(data[,1] > t[1] & data[,1] < t[2])
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  #
  par[9] <- max(peak.y)
  par[8] <- peak.t[which(peak.y==par[9])]
  par[10] <- 0.25
  #
  rm(peak.w,peak.t,peak.y,t)
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,par[8:10])
  
  plot(data[,1],excess)
  lines(data[,1],residue)
  

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


# Refine Parameters

for (i in 1:nBeats){
  seg <- c(beat[i,3],0,beat[i,4])
  data <- model2.GetSegment(ppg,seg)
  rm(seg)
  
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
#Plotting all the widths and plotting them (S circles, D triangles):
#plot(beat$SWidth, ylim=c(0, 0.5))
#points(beat$DWidth, pch=17)
#points(beat$NWidth, pch=15)
