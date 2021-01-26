## Constants
#

lab.time = "time (s)"
lab.ppg = "Detrended"
lab.dt = "interval (s)"

config.rate = 0.95
config.a.b = 0.999
config.optimum.width = 0.3
config.penalty.width = 10

const.pi = 3.1415926535897932384626433

#
## End of constants section

model2.LoadBeat <- function(file){
  filename <- file
  if (regexpr("_intervals.csv$",file)<0){
    if (regexpr(".csv$",file)>=0){
      filename <- substring(file,0,nchar(file)-4)
    }
    filename <- paste(filename,"_intervals.csv",sep="")
  }

  result <- read.csv(filename)
  names(result)[1] <- lab.time
  names(result)[2] <- lab.dt

  return(result)
}

model2.LoadPPG <- function(file){
  filename <- file
  if (regexpr(".csv$",file)<0){
    filename <- paste(filename,".csv",sep="")
  }

  result <- read.csv(filename)
  names(result)[1] <- lab.time
  names(result)[2] <- lab.ppg

  return(result)
}

model2.FixParams <- function(data,params,debug=FALSE){
  temp <- model2.FIX_PAR( data, params, debug )
  return( temp[2:length(temp)] )#excludes first value which is the penalty
}

model2.FIX_PAR <- function(data,params,debug=FALSE){
  # Params:
  # Baseline, {LateBaseline}, t_S, h_S, w_S, t_D, h_D, w_D, {t_R, h_R, w_R}
  par <- as.numeric(params)
  
  nData <- nrow(data)
  nPar <- length(params)
  
  # Transcribe parameters
  nBase <- 1
  baseline <- c( par[1], par[1] )
  if (nPar == 8|nPar == 11){ #this means we provided 2 baselines, if we had only one it would be 7 or 10; ||nPar = 5
    baseline[2] = par[2]
    nBase <- 2
  }
  
  t <- c( par[nBase + 1], 0, 0 )
  h <- c( par[nBase + 2], 0, 0 )
  w <- c( par[nBase + 3], 0, 0 )
  hasPeak <- c( TRUE, FALSE, FALSE )
  # Calculate penalty only if a peak has been supplied
  
  if (nPar >= nBase + 6){
    hasPeak[2] = TRUE
    t[2] <- par[nBase + 4]
    h[2] <- par[nBase + 5]
    w[2] <- par[nBase + 6]
  }
  
  if (nPar >= nBase + 9){ # then this is the renal peak and we will call it peak 3
    hasPeak[3] = TRUE
    t[3] <- par[nBase + 7]
    h[3] <- par[nBase + 8]
    w[3] <- par[nBase + 9]
  }
  
  # Clamp and/or penalize parameters
  penalty <- 0
  tMin <- data[1,1]
  tMax <- data[nData,1]
  
  META_BASELINE_SHIFT <- 1.0 #might be sensible to adjust - a penalty for how big the gap is between the two baselines
  META_MIN_PEAK_DELAY <- 0.1 #the renal peak can't be less than 100ms after systolic, the diastolic can't be less than 100ms from renal...
  
  p <- 1:11*0
  
  # Fix height, width
  for (i in 1:3)
  {
    if (h[i] < 0){#if the amplitude is less than 0 we throw it away
      penalty <- penalty + h[i]*h[i]
      p[3*i-1] <- h[i]*h[i]
      
      h[i] <- 0
    }
    
    if (w[i] < 0.05 | w[i] > 0.5){#kind of metaparameters - if the width is less than 0.05 or bigger than 0.5 we throw it away
      fixed <- max( 0.05, min( w[i], 0.5 ) )
      diff <- fixed - w[i]
      
      penalty <- penalty + diff*diff
      p[3*i-0] <- diff*diff #p is for debug mode, giving a breakdown of the penalties
      
      w[i] <- fixed
    }
  }
  
  # Fix time
  #SYSTOLIC PEAK
  fixed <- max( tMin, min( t[1], tMin + 1, tMax ) ) #first peak can't be before the data or more than a second after the datat starts
  if (debug){
    print(paste("time S: ",tMin," < ",t[1]," < min( ",tMin+1,",",tMax," )"))
  }
  if (t[1] != fixed){
    diff <- fixed - t[1]
    
    penalty <- penalty + diff*diff
    p[1] <- diff*diff
    
    t[1] <- fixed
  }
  #DIASTOLIC PEAK
  fixed <- max( 2 * META_MIN_PEAK_DELAY, min( t[2], tMax - tMin + 0.4 * w[2] ) )#diastolic peak has to be at least 2 minimum peak delays after the systolic peak, i.e. numerically they have to be distinguishable peaks (regardless of physiological plausibility)
  if (debug){
    print(paste("time D: ",2 * META_MIN_PEAK_DELAY," < ",t[2]," < ",tMax - tMin + 0.4 * w[2]," )"))#tMax-tMin plus 40% of peak width (this is a delay)
  }
  if (t[2] != fixed){
    # Two peak delays between S and D
    diff <- fixed - t[2]
    
    if (hasPeak[2]){ #if we changed the vlaue of the time
      penalty <- penalty + diff*diff
      p[4] <- diff*diff
    }
    
    t[2] <- fixed
  }
  #RENAL PEAK
  fixed <- max( META_MIN_PEAK_DELAY, min( t[3], t[2] - META_MIN_PEAK_DELAY ) ) #boundaries are the min peak delay and the diastolic peak MINUS one peak delay
  if (debug){
    print(paste("time R: ",META_MIN_PEAK_DELAY," < ",t[3]," < ",t[2] - META_MIN_PEAK_DELAY," )"))
  }
  if (t[3] != fixed){
    diff <- fixed - t[3]
    
    if (hasPeak[3]){
      penalty <- penalty + diff*diff
      p[7] <- diff*diff
    }
    
    t[3] <- fixed
  }
  
  # Don't clamp baseline shift, but penalize large shifts
  diff <- abs(baseline[1] - baseline[2]) #penalty for baselines differing
  penalty <- penalty + META_BASELINE_SHIFT*diff*diff
  
  fixedPar <- c( baseline, t[1], h[1], w[1], t[2], h[2], w[2], t[3], h[3], w[3] )#assembling all params into new parameter vector, always will have 2 baselines as rebuild function expects it

  if (debug){
    print(p)
  }
  
  return( c( penalty, fixedPar ) )
}

model2.ChiSq <- function(data,params,debug=FALSE){
  temp <- model2.FIX_PAR( data, params, debug )
  penalty <- temp[1] #penalty is the first element of the array
  fixedPar <- temp[2:length(temp)] #the rest of the array are the fixed parameters
  rm(temp)

  nData <- nrow(data)
  nPar <- length(params)

  fit <- model2.Rebuild(data,data[1,2],fixedPar) #calculating what the PPG should look like based on the parameters
  residue <- data[,2] - fit #modelled data compared to original = residue

  if (debug){
    plot(data[,1],data[,2])
    lines(data[,1],fit)
    print(sum(residue*residue))
    print(p)
  }

  return(sum(residue*residue) / (nData-nPar) + as.numeric(penalty)) #if number of parameters increases, penalty is slightly weightier now
}

model2.FindPeak <- function(ppg,time){
  nppg <- nrow(ppg)
  
  w0 <- min(which(ppg[[lab.time]]>time))
  wLim <- c(max(1,w0-75),min(nppg,w0+75))
  
  a <- c(0,w0 - wLim[1],0)
  a[1] <- max(1,a[2]-20)
  a[3] <- a[2]+min(nppg-w0,10)
  
  t_time <- ppg[[lab.time]][wLim[1]:(wLim[2]-1)]
  t_ppg <- ppg[[lab.ppg]][wLim[1]:(wLim[2]-1)]
  grad <- ppg[[lab.ppg]][(wLim[1]+1):wLim[2]] - ppg[[lab.ppg]][wLim[1]:(wLim[2]-1)]
  w0 <- min(which(grad[a[1]:a[3]] == max(grad[a[1]:a[3]])) + wLim[1] + a[1] - 1)
  
  while (w0 > wLim[1] && grad[w0-wLim[1]-1] > grad[w0-wLim[1]]){
    w0 <- w0 - 1
  }
  while (w0 < wLim[2] && grad[w0-wLim[1]+1] > grad[w0-wLim[1]]){
    w0 <- w0 + 1
  }
  
  wLow <- w0-1
  while (wLow > wLim[1] && ppg[[lab.ppg]][wLow-1] < ppg[[lab.ppg]][wLow]){
    wLow <- wLow - 1
  }
  if (wLow > wLim[1]){
    wLow <- wLow - 1
  }
  
  return( c(wLow,w0) )
}

model2.FindSegment <- function(ppg,beatTime,nextBeatTime=NA){
  nppg <- nrow(ppg)

  w <- model2.FindPeak( ppg, beatTime )
  wLow <- w[1]
  w0 <- w[2]
  if (!is.na(nextBeatTime)){
    w <- model2.FindPeak( ppg, nextBeatTime )
    wHigh <- w[1]
  }else{
    wHigh <- min(nppg,w0 + 75)
    grad <- ppg[[lab.ppg]][(wLow+1):wHigh] - ppg[[lab.ppg]][wLow:(wHigh-1)]
    split <- 2*(w0-wLow)
    peakVal <- max(grad[1:split])
    nextPeakVal <- max(grad[split:length(grad)])
    nextPeak <- which(grad[split:length(grad)] == nextPeakVal) + split - 1

    if (nextPeakVal > 0.75 * peakVal){
      # Assume we're trimming off an incomplete beat
      w <- model2.FindPeak( ppg, ppg[[lab.time]][wLow+nextPeak] )
      wHigh <- w[1]
    } else {
      wHigh <- min(nppg,w0 + 75)
    }
  }

  return( c(wLow,w0,wHigh) )
}

model2.GetSegment <- function(ppg,limits){
  w <- c(limits[1]:limits[3])

  result <- matrix(nrow=length(w),ncol=2)

  result[,1] <- ppg[[lab.time]][w]
  result[,2] <- ppg[[lab.ppg]][w]

  return(result)
}

model2.EstimateParams <- function(xy,yPrev,beat){
  dt <- xy[2,1]-xy[1,1]
  
  #beat<-beats[i,] # For debugging
  #data <- NULL
  beatTime <- beat[[lab.time]]

  result <- 1:10
  result[1] <- min(xy[,2])-0.326

  residue <- model2.Excess(xy[,2],yPrev,result[1])
  #plot(xy[,1],residue)

  ## S peak
  peak.w <- which(xy[,1] > beatTime-0.2 & xy[,1] < beatTime+0.2)
  peak.t <- xy[peak.w,1]
  peak.y <- residue[peak.w]

  result[3] <- max(peak.y)
  result[2] <- peak.t[which(peak.y==result[3])]
  result[4] <- 0.25
  # fit <- model2.Rebuild(xy,yPrev,result[1:4])
  # lines(xy[,1],fit)
  
  sim <- simplex.MakeSimplex(xy,result[1:4],model2.ChiSq,0.1)
  sim <- simplex.Run(xy,sim,model2.ChiSq)
  fit <- model2.Rebuild(xy,yPrev,sim[1,])
  #plot(xy[,1],xy[,2])
  #lines(xy[,1],fit)
  #stop()
  result[1:4] = sim[1,]
  residue <- model2.SubtractExcessPeak(xy[,1],residue,result[2:4])
  
  ## D peak
  peak.w <- which(xy[,1] > beatTime+0.2 & xy[,1] < beatTime+0.6)
  peak.t <- xy[peak.w,1]
  peak.y <- residue[peak.w]
  
  #plot(xy[,1],residue)
  #lines(peak.t,peak.y)
  
  result[6] <- max(peak.y)
  result[5] <- peak.t[which(peak.y==result[6])]-result[2]
  result[7] <- 0.25
  # fit <- model2.Rebuild(xy,yPrev,result[1:7])
  # lines(xy[,1],fit)
  
  sim <- simplex.MakeSimplex(xy,c(result[1:7]),model2.ChiSq,0.1)
  sim <- simplex.Run(xy,sim,model2.ChiSq)
  result[1:7] = sim[1,]
  fit <- model2.Rebuild(xy,yPrev,sim[1,])
  #print("D")
  #plot(xy[,1],xy[,2])
  #lines(xy[,1],fit)
  #stop()
  residue <- model2.SubtractExcessPeak(xy[,1]-result[2],residue,result[5:7])
  
  ## N peak
  peak.w <- which(xy[,1] > result[2] + 0.25 * result[5] & xy[,1] < result[2] + 0.75 * result[5])
  peak.t <- xy[peak.w,1]
  peak.y <- residue[peak.w]
  
  result[8] <- max(peak.y)
  result[9] <- peak.t[which(peak.y==result[8])]-result[2]
  result[10] <- 0.25
  
  sim <- simplex.MakeSimplex(xy,c(result[1:10]),model2.ChiSq,0.1)
  sim <- simplex.Run(xy,sim,model2.ChiSq)
  result = sim[1,]
  fit <- model2.Rebuild(xy,yPrev,sim[1,])
  #print("N")
  #plot(xy[,1],xy[,2])
  #lines(xy[,1],fit)
  #stop()
  residue <- model2.SubtractExcessPeak(xy[,1]-result[2],residue,result[8:10])
  
  residue <- model2.Excess(xy[,2],yPrev,result[1])
  residue <- model2.SubtractExcessPeak(xy[,1],residue,result[2:4])
  residue <- model2.SubtractExcessPeak(xy[,1]-result[2],residue,result[5:7])
  residue <- model2.SubtractExcessPeak(xy[,1]-result[2],residue,result[8:10])
  #plot(xy[,1],xy[,2])
  #lines(xy[,1],fit)
  #stop()

  # Fix possibly broken values
  result[3] <- max( 0.05, result[3] )
  result[6] <- max( 0.05, result[6] )
  result[9] <- max( 0.05, result[9] )
  if (result[5] < 0.15){
    result[5] <- 0.30
  }
  if (result[8] < 0.05){
    result[8] <- 0.15
  }

  return(result)
}

model2.Excess <- function(y,offset,baseline){
  count <- length(y)
  if (length(offset) == 0){
    print("Help")
  }

  result <- 1:count * 0.0

  result[1] = y[1] - (baseline + config.rate*(offset-baseline))

  for (j in 2:count){
    result[j] = y[j] - (baseline + config.rate*(y[j-1]-baseline))
  }

  return(result)
}

model2.Excess.Inv <- function(time,excess,offset,baselineStart,baselineEnd,timeBase){
  nX <- length(excess)
  if (nX == 0){
    print("Help")
  }

  result <- 1:nX * 0.0
  
  baseline <- time * 0 + baselineStart
  baseline[which(time > timeBase)] = baselineEnd

  result[1] = excess[1] + (baselineStart + config.rate*(offset-baselineStart))

  for (j in 2:nX){
    result[j] = excess[j] + (baseline[j] + config.rate*(result[j-1]-baseline[j]))
  }
  return(result)
}

model2.Peak <- function(time,peakParams){
  temp <- 2*const.pi*(time - as.double(peakParams[1]))/as.double(peakParams[3])
  temp[which(temp < -const.pi)] = -const.pi
  temp[which(temp >  const.pi)] =  const.pi
  result <- as.double(peakParams[2]) * (0.5 * (1+cos(temp)))^2

  return(result)
}

model2.Plot <- function(xy,yPrev,params){
  plot(xy[,1],model2.Excess(xy[,2],yPrev,params[1]),xlab="time (s)",ylab="Excess PPG")
  for (iPeak in 1:3){
    time <- xy[,1]
    fit <- model2.Peak(xy[,1],params[(3*iPeak-1):(3*iPeak+1)])
    w <- which(fit > 0.01 * max(fit))
    lines(time[w],fit[w],lty="dotted")
  }
  #lines(xy[,1],model2.Peak(xy[,1],params[2:4]),lty="dotted")
  #lines(xy[,1],model2.Peak(xy[,1],params[5:7]),lty="dotted")
  #lines(xy[,1],model2.Peak(xy[,1],params[8:10]),lty="dotted")
  lines(xy[,1],model2.Residue(xy,params))
}

model2.Rebuild <- function(xy,offset,params,invert=TRUE){
  result <- 1:nrow(xy) * 0.0
  result <- result + model2.Peak(xy[,1],params[3:5])
  if (length(params)>=8){
    result <- result + model2.Peak(xy[,1],params[6:8]+c(params[3],0,0))
  }
  if (length(params)>=11){
    result <- result + model2.Peak(xy[,1],params[9:11]+c(params[3],0,0))
  }

  if (invert){
    result <- model2.Excess.Inv(xy[,1],result,offset,params[1],params[2],params[3]+0.5*params[6]) #params[3] is the systolic time, params[6] is the diastolic time delay, so we switch baseline halfway between the s and d peak
  }

  return(as.double(result))
}

model2.Residue <- function(xy,params){
  result <- model2.Excess(xy[,2],params[1],params[1])
  result <- model2.SubtractExcessPeak(xy[,1],result,params[2:4])
  if (length(params) >= 7){
    result <- model2.SubtractExcessPeak(xy[,1],result,params[5:7]+c(params[2],0,0))
  }
  if (length(params) >= 10){
    result <- model2.SubtractExcessPeak(xy[,1],result,params[8:10]+c(params[2],0,0))
  }
  return(result)
}

model2.Load <- function(filename){
  if (file.exists(filename)){
    params <- read.csv(filename)
    # Erase line index column
    params["X"] <- NULL

    return(data.matrix(params))
  }

  return(NULL)
}

model2.Save <- function(filename,params){
  names <- c(
    "Start Time (s)","Segment start","Start PPG",
    "Baseline",
    "S-peak time(s)","S-peak amplitude","S-peak width",
    "D-peak time(s)","D-peak amplitude","D-peak width",
    "N-peak time(s)","N-peak amplitude","N-peak width",
    "chi-squared"
  )
  output <- as.data.frame(params)
  names(output) <- names
  write.csv(output,file=filename,row.names=FALSE)
  return(0)
}

model2.SubtractExcessPeak <- function(time,residue,peakParams){
  result <- residue - model2.Peak(time,peakParams)

  return(result)
}

####################
model2.ChiSq2 <- function(data,params,debug=FALSE){
  temp <- model2.FIX_PAR2( data, params, debug )
  penalty <- temp[1] #penalty is the first element of the array
  fixedPar <- temp[2:length(temp)] #the rest of the array are the fixed parameters
  rm(temp)
  
  nData <- nrow(data)
  nPar <- length(params)
  
  fit <- model2.Rebuild(data,data[1,2],fixedPar) #calculating what the PPG should look like based on the parameters
  residue <- data[,2] - fit #modelled data compared to original = residue
  
  if (debug){
    plot(data[,1],data[,2])
    lines(data[,1],fit)
    print(sum(residue*residue))
    print(p)
  }
  
  return(sum(residue*residue) / (nData-nPar) + as.numeric(penalty)) #if number of parameters increases, penalty is slightly weightier now
}

model2.FIX_PAR2 <- function(data,params,debug=FALSE){
  # Params:
  # Baseline, {LateBaseline}, t_S, h_S, h_D, {{h_R}}
  par <- as.numeric(params)
  
  nData <- nrow(data)
  nPar <- length(params)
  
  # Transcribe parameters
  nBase <- 1
  baseline <- c( par[1], par[1] )
  if (nPar == 5|nPar == 6){ 
    baseline[2] = par[2]
    nBase <- 2
  }
  
  t <- c( par[nBase + 1], 0.285, 0.17)
  h <- c( par[nBase + 2], par[nBase + 3], 0)#passing in h_D
  w <- c( par[0.2], 0.12, 0.12 )#there is a vector to pass in params that were determined elsewhere
  hasPeak <- c( TRUE, TRUE, FALSE ) #s, d, r
  # Calculate penalty only if a peak has been supplied
  
  if (nPar >= nBase + 4){ # then this is the renal peak and we will call it peak 3
    hasPeak[3] = TRUE
    h[3] <- par[nBase + 4]
  }
  
  # Clamp and/or penalize parameters
  penalty <- 0
  tMin <- data[1,1]
  tMax <- data[nData,1]
  
  META_BASELINE_SHIFT <- 1.0 #might be sensible to adjust - a penalty for how big the gap is between the two baselines
  META_MIN_PEAK_DELAY <- 0.1 #the renal peak can't be less than 100ms after systolic, the diastolic can't be less than 100ms from renal...
  
  p <- 1:11*0
  
  # Fix height
  for (i in 1:3)
  {
    if (h[i] < 0){#make sure the amplitude is non-negative
      penalty <- penalty + h[i]*h[i]
      p[3*i-1] <- h[i]*h[i]
      
      h[i] <- 0
    }
  }
  
  # Fix time
  #SYSTOLIC PEAK
  fixed <- max( tMin, min( t[1], tMin + 1, tMax ) ) #first peak can't be before the data or more than a second after the datat starts
  if (debug){
    print(paste("time S: ",tMin," < ",t[1]," < min( ",tMin+1,",",tMax," )"))
  }
  if (t[1] != fixed){
    diff <- fixed - t[1]
    
    penalty <- penalty + diff*diff
    p[1] <- diff*diff
    
    t[1] <- fixed
  }
  # Don't clamp baseline shift, but penalize large shifts
  diff <- abs(baseline[1] - baseline[2]) #penalty for baselines differing
  penalty <- penalty + META_BASELINE_SHIFT*diff*diff
  
  fixedPar <- c( baseline, t[1], h[1], w[1], t[2], h[2], w[2], t[3], h[3], w[3] )#assembling all params into new parameter vector, always will have 2 baselines as rebuild function expects it
  
  if (debug){
    print(p)
  }
  
  return( c( penalty, fixedPar ) )
}

model2.FixParams2 <- function(data,params,debug=FALSE){
  temp <- model2.FIX_PAR2( data, params, debug )
  return( temp[2:length(temp)] )#excludes first value which is the penalty
}