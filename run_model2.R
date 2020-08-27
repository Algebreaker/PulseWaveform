source("model2.R")
source("simplex.R")

input = "data/66102058a16"
#input = "data/Paul"

ppg <- model2.LoadPPG(input)
beat <- model2.LoadBeat(input)
nBeats = nrow(beat)

params <- matrix(nrow=nBeats,ncol=13)

count <- nrow(ppg)
fit <- 1:count * 0
residue <- 1:count * 0
baseline <- 1:count * 0
tLim <- c(1e32,-1e32)

params <- model2.Load(paste(input,"_params_a.csv",sep=""))

if (is.null(params)){
  print("Estimating fitting parameters...")
  seg <- c(0,0)
  all <- c(0,0)

  params <- matrix(nrow=nBeats,ncol=13)
  
  startClock <- Sys.time()
  reportClock <- 10
  for (i in 1:nBeats){
    if (as.double(difftime(Sys.time(),startClock,unit="secs"))>reportClock){
      print(paste("Processing beat",i))
      reportClock <- reportClock + 10
    }
    
    if (seg[1] == 0){
      seg <- model2.FindSegment(ppg,beat,i)
      all <- seg
    } else {
      tSeg <- model2.FindSegment(ppg,beat,i)
      seg[1] <- seg[2] + 1
      seg[2] <- tSeg[2]
      all[2] <- tSeg[2]
    }
    data <- model2.GetSegment(ppg,seg)
    tStart <- ppg[seg[1],1]
    yPrev <- ppg[max(seg[1]-1,1),2]
    
    tLim[1] <- min(tLim[1],tStart)
    tLim[2] <- max(tLim[2],ppg[seg[2],1])
  
    count <- seg[2] + 1 - seg[1]
    data[,1] <- ppg[seg[1],1] + (1:count - 1) / 75.0
  
    par <- model2.EstimateParams(data,yPrev,beat[i,])

    sim <- simplex.MakeSimplex(data,par,model2.ChiSq,0.1,-1)
    sim <- simplex.Run(data,sim,model2.ChiSqPPG)
    par <- sim[1,]

    params[i,1] <- tStart
    params[i,2] <- yPrev
    params[i,3:12] <- par
    params[i,13] <- model2.ChiSqPPG(data,par)

    wSeg <- seg[1]:seg[2]
    fit[wSeg] <- model2.Rebuild(data,yPrev,par)
    residue[wSeg] <- ppg[wSeg,2] - fit[wSeg]
    baseline[wSeg] <- params[i,3]
  }

  print("Completed initial parameter generation")
  model2.Save(paste(input,"_params_a.csv",sep=""),params)
  params <- model2.CleanParams(params)
} else {
  print("Using estimated fitting parameters...")
  if (ncol(params) == 12){
    temp <- matrix(nrow=nrow(params),ncol=13)
    temp[,1:12] <- params
    temp[,13] <- 0
    params <- temp
  }
}

# Limit variability in timing of D and N peaks
fixedParams = c(1.0,mean(params[,7]),sd(params[,7]),mean(params[,10]),sd(params[,10]))

seg <- c(0,0)
all <- c(0,0)

startClock <- Sys.time()
reportClock <- 10

print("Refining fitting parameters...")
for (i in 1:nBeats){
  if (as.double(difftime(Sys.time(),startClock,unit="secs"))>reportClock){
    print(paste("Processing beat",i))
    reportClock <- reportClock + 10
  }

  if (seg[1] == 0){
    seg <- model2.FindSegment(ppg,beat,i)
    all <- seg
  } else {
    tSeg <- model2.FindSegment(ppg,beat,i)
    seg[1] <- seg[2] + 1
    seg[2] <- tSeg[2]
    all[2] <- tSeg[2]
  }
  data <- model2.GetSegment(ppg,seg)
  tStart <- ppg[seg[1],1]
  yPrev <- ppg[max(seg[1]-1,1),2]

  tLim[1] <- min(tLim[1],tStart)
  tLim[2] <- max(tLim[2],ppg[seg[2],1])
  
  count <- seg[2] + 1 - seg[1]
  data[,1] <- ppg[seg[1],1] + (1:count - 1) / 75.0
  
  par <- params[i,3:12]
  # Can't make a simplex if peak amplitude is zero
  par[3] <- max(0.1,par[3])
  par[6] <- max(0.1,par[6])
  par[9] <- max(0.1,par[9])
  
  sim <- simplex.MakeSimplex(data,par[1:7],model2.ChiSqZeroN,0.1,-1,optional=fixedParams)
  sim <- simplex.Run(data,sim,model2.ChiSqZeroN,optional=fixedParams)
  best <- c(sim[1,],0,0,0)
  chiSqSD <- model2.ChiSq(data,best)

  sim <- simplex.MakeSimplex(data,par,model2.ChiSqPPGBound,0.1,-1,optional=fixedParams)
  sim <- simplex.Run(data,sim,model2.ChiSqPPGBound,optional=fixedParams)
  chiSqSND <- model2.ChiSq(data,sim[1,])

  # Accept N peak only if it improves things without disrupting the S peak
  if (chiSqSND < chiSqSD & sim[1,9] > 0 & sim[1,3] > 0.1 * sim[1,9]){
    best <- sim[1,]
  }
 
  chiSq2 <- model2.ChiSq(data,best)
  params[i,1] <- tStart
  params[i,2] <- yPrev
  params[i,3:12] <- best
  params[i,13] <- chiSq2

  if (chiSq2 > 2){
    print(paste("Beat",i,"has chi-squared of",chiSq2))
  }
  
  wSeg <- seg[1]:seg[2]
  fit[wSeg] <- model2.Rebuild(data,yPrev,best)
  residue[wSeg] <- ppg[wSeg,2] - fit[wSeg]
  baseline[wSeg] <- params[i,3]
}
print("Completed parameter generation")

model2.Save(paste(input,"_params_b.csv",sep=""),params)

# Make 'smooth' time axis
time <- ppg[1,1] + (1:nrow(ppg) - 1) / 75
w <- all[1]:all[2]
plot(time,ppg[,2],xlim=tLim,xlab="time (s)",ylab="PPG count")
lines(time[w],fit[w])
lines(time[w],baseline[w],lty="dashed")
lines(time[w],88+residue[w],lty="dotted")
