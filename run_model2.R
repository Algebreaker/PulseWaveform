source("model2.R")
source("simplex.R")

input = "data/66102058a16"
#input = "data/Toby"

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
  seg <- c(0,0,0)
  all <- c(0,0)

  if (is.null(params)){
    params <- matrix(nrow=nBeats,ncol=13)
  }

  startClock <- Sys.time()
  reportClock <- 10
  for (i in 1:nBeats){
    #print(paste("Segment",i))
    if (as.double(difftime(Sys.time(),startClock,unit="secs"))>reportClock){
      min <- as.integer(reportClock/60)
      sec <- reportClock %% 60
      time <- paste("[",min,":",sep="")
      if (sec == 0){ time <- paste(time,"00]",sep="")} else { time <- paste(time,sec,"]",sep="")}
      print(paste(time,"Processing beat",i))
      reportClock <- reportClock + 10
    }
    
    if (seg[1] == 0){
      seg <- model2.FindSegment(ppg,beat,i)
      all <- c(seg[1],seg[3])
    } else {
      tSeg <- model2.FindSegment(ppg,beat,i)
      seg[1] <- seg[3] + 1
      seg[2] <- tSeg[2]
      seg[3] <- tSeg[3]
      all[3] <- tSeg[3]
    }
    data <- model2.GetSegment(ppg,seg)
    tStart <- ppg[seg[1],1]
    yPrev <- ppg[max(seg[1]-1,1),2]
    
    tLim[1] <- min(tLim[1],tStart)
    tLim[2] <- max(tLim[2],ppg[seg[3],1])
  
    count <- seg[3] + 1 - seg[1]
    data[,1] <- ppg[seg[1],1] + (1:count - 1) / 75.0

    par <- model2.EstimateParams(data,yPrev,beat[i,])
   
    for (iRepeat in 1:2){
      sim <- simplex.MakeSimplex(data,par,model2.ChiSq,0.1)
      
      chiSqMin <- model2.ChiSq(data,sim[2,])
      argMin <- 2
      for (j in 3:nrow(sim)){
        chiSq <- model2.ChiSq(data,sim[j,])
        if (chiSq < chiSqMin){
          chiSqMin <- chiSq
          argMin <- j
        }
      }

      sim <- simplex.Run(data,sim,model2.ChiSq)
      fit <- model2.Rebuild(data,yPrev,sim[1,])
      par <- sim[1,]
    }

    if (par[8] > par[5]){
      temp <- par[5:7]
      par[5:7] <- par[8:10]
      par[8:10] <- temp
    }
 
    params[i,1] <- tStart
    params[i,2] <- yPrev
    params[i,3:12] <- par
    params[i,13] <- model2.ChiSq(data,par)

    wSeg <- seg[1]:seg[3]
    fit[wSeg] <- model2.Rebuild(data,yPrev,par)
    residue[wSeg] <- ppg[wSeg,2] - fit[wSeg]
    baseline[wSeg] <- params[i,3]
  }

  print("Completed initial parameter generation")
  model2.Save(paste(input,"_params_a.csv",sep=""),params)
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
fixedParams = c(1.0,median(params[,7]),sd(params[,7]),median(params[,10]),sd(params[,10]))

seg <- c(0,0,0)
all <- c(0,0)

startClock <- Sys.time()
reportClock <- 10

print("Refining fitting parameters...")
for (i in 1:nBeats){
  if (as.double(difftime(Sys.time(),startClock,unit="secs"))>reportClock){
    min <- as.integer(reportClock/60)
    sec <- reportClock %% 60
    time <- paste("[",min,":",sep="")
    if (sec == 0){ time <- paste(time,"00]",sep="")} else { time <- paste(time,sec,"]",sep="")}
    print(paste(time,"Processing beat",i))
    reportClock <- reportClock + 10
  }
  
  #print(paste("Refine",i))

  if (seg[1] == 0){
    seg <- model2.FindSegment(ppg,beat,i)
    all <- c(seg[1],seg[3])
  } else {
    tSeg <- model2.FindSegment(ppg,beat,i)
    seg[1] <- seg[3] + 1
    seg[2] <- tSeg[2]
    seg[3] <- tSeg[3]
    all[2] <- tSeg[3]
  }
  data <- model2.GetSegment(ppg,seg)
  
  par <- params[i,3:12]

  tStart <- ppg[seg[1],1]
  yPrev <- ppg[max(seg[1]-1,1),2]

  tLim[1] <- min(tLim[1],tStart)
  tLim[2] <- max(tLim[2],ppg[seg[3],1])
  
  count <- seg[3] + 1 - seg[1]
  data[,1] <- ppg[seg[1],1] + (1:count - 1) / 75.0
  
  par <- params[i,3:12]
  # Clamp amplitudes to be non-negative
  par[3] <- max(0.,par[3])
  par[6] <- max(0.,par[6])
  par[9] <- max(0.,par[9])
  
  # If amplitudes are zero, simplex construction can't
  # be along coordinate directions because peak position
  # will have no effect on chi-squared.  Thus, choose
  # directions which increase the amplitude while moving
  # the peak
  dir <- diag(10)
  if (par[3] < 0.1){
    dir[2:3,2:3] <- 0.707 * c(1,1,1,-1)
  }
  if (par[6] < 0.1){
    dir[5:6,5:6] <- 0.707 * c(1,1,1,-1)
  }
  if (par[3] < 0.1){
    dir[8:9,8:9] <- 0.707 * c(1,1,1,-1)
  }

  if (abs(par[5]-par[8]) < 0.1)
  {
    if (par[6] < 0){
      par[5:7] <- par[8:10]
    } else if (par[9] > 0) {
      par[5] <- (par[6]*par[5] + par[9]*par[8]) / (par[6]+par[9])
      par[7] <- (par[6]*par[7] + par[9]*par[10]) / (par[6]+par[9])
      par[6] <- par[6] + par[9]
    }
    par[8:10] <- c(-1,0,0)
  }
  
  if (par[5] < 0 | par[5] > data[nrow(data),1] - par[2] | (par[8]<0 & par[5] < 0.25)){
    if (par[8] > 0){
      par[5:7] <- par[8:10]
    }
    par[9] <- 0

    for (iRepeat in 1:2){
      sim <- simplex.MakeSimplex(data,par[1:7],model2.ChiSq,0.1,directions = dir[1:7,1:7])
      sim <- simplex.Run(data,sim,model2.ChiSq)
      par[1:7] < sim[1,]
    }

    fakeD <- c(fixedParams[2],0,0.5) # Zero amplitude, no penalties
    best <- c(sim[1,1:4],fakeD,sim[1,5:7])
  } else {
    peakD <- par[5:7]
    peakN <- par[8:10]
    fPar <- fixedParams
    if (peakD[1] < peakN[1]){
      peakD <- par[8:10]
      peakN <- par[5:7]
      par[5:7] <- peakD
      par[8:10] <- peakN
      fPar[2:3] <- fixedParams[4:5]
      fPar[4:5] <- fixedParams[2:3]
    }
    
    chiSqMin <- model2.ChiSq(data,par)
    chiSqN <- model2.ChiSq(data,c(par[1:4],peakN))
    chiSqD <- model2.ChiSq(data,c(par[1:4],peakD))

    if (chiSqN < chiSqD){
      par <- c(par[1:4],peakN,0,0,0)
    } else {
      par <- c(par[1:4],peakD,0,0,0)
    }
    
    sim <- simplex.MakeSimplex(data,par[1:7],model2.ChiSq,0.1)
    if (is.character(sim)){
      print("MakeSimplex --debug")
      iParam <- -1
      if (regexpr("^Error: param\\[[0-9]+\\]$",sim)){
        iParam <- as.integer(substr(sim,nchar("Error: param[")+1,regexpr("\\]",sim)-1))
        print(paste("Error in parameter",iParam))
        print(par[1:7])
        plot(data[,1],data[,2])
        fit <- model2.Rebuild(data,data[1,2],par[1:7])
        lines(data[,1],fit)
        tPar <- par[1:7]
        tPar[iParam] <- tPar[iParam] + 1
        fit <- model2.Rebuild(data,data[1,2],tPar)
        lines(data[,1],fit,lty="dotted")
        stop()
      }
      simplex.MakeSimplex(data,par[1:7],model2.ChiSq,0.1,optional=fixedParams,debug=TRUE)
      print(paste("Stopping: beat",i))
      stop()
    }

    sim <- simplex.Run(data,sim,model2.ChiSq,optional=fPar)
    par[1:7] <- sim[1,]

    if (sim[1,5] < 0.5 * (peakN[1] + peakD[1])){
      par <- c(sim[1,1:4],peakD,sim[1,5:7])
    } else {
      par <- c(sim[1,],peakN)
    }
    chiSqSD <- model2.ChiSq(data,par)
    best <- par

    sim <- simplex.MakeSimplex(data,par,model2.ChiSq,0.1,optional=fixedParams)
    if (is.character(sim)){
      print("MakeSimplex --debug")
      iParam <- -1
      if (regexpr("^Error: param\\[[0-9]+\\]$",sim)){
        iParam <- as.integer(substr(sim,nchar("Error: param[")+1,regexpr("\\]",sim)-1))
        print(paste("Error in parameter",iParam))
        print(par)
        plot(data[,1],data[,2])
        fit <- model2.Rebuild(data,data[1,2],par)
        lines(data[,1],fit)
        tPar <- par
        tPar[iParam] <- tPar[iParam] + 1
        fit <- model2.Rebuild(data,data[1,2],tPar)
        lines(data[,1],fit,lty="dotted")
        #stop()
      }
      simplex.MakeSimplex(data,par,model2.ChiSq,0.1,optional=fixedParams,debug=TRUE)
      print(paste("Stopping: beat",i))
      stop()
    }

    sim <- simplex.Run(data,sim,model2.ChiSq,optional=fixedParams)
    chiSqSND <- model2.ChiSq(data,sim[1,])
  
    # Accept N peak only if it improves things without disrupting the S peak
    if (chiSqSND < chiSqSD & sim[1,9] > 0 & sim[1,3] > 0.1 * sim[1,9]){
      best <- sim[1,]
    }
  }

  chiSq2 <- model2.ChiSq(data,best)
  params[i,1] <- tStart
  params[i,2] <- yPrev
  params[i,3:12] <- best
  params[i,13] <- chiSq2

  if (chiSq2 > 2){
    print(paste("Beat",i,"has chi-squared of",chiSq2))
  }
  
  wSeg <- seg[1]:seg[3]
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
lines(time[w],max(ppg[,2])+residue[w],lty="dotted")
