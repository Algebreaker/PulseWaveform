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

model2.ChiSq <- function(data,params,optional=NULL,debug=FALSE){
  nPar <- length(params)
  fixed <- params
  penalty <- 0
  tMax <- data[nrow(data),1]
  
  p <- 1:10*0
  
  if (!is.null(optional) & length(optional) == 5){
    meanDTime = optional[2]
    sdDTime = optional[3]
    meanNTime = optional[4]
    sdNTime = optional[5]
  }

  if (params[3] < 0){
    fixed[3] <- 0
    penalty <- penalty + params[3]*params[3]
    p[1] <- params[3]*params[3]
  }
  
  if (nPar >= 7){
    if (!is.null(optional) & length(optional) >= 3){
      meanTime <- optional[2]
      sdTime <- optional[3]
     
      if (nPar == 7 & length(optional) >= 5){
        if (abs((params[5] - optional[4])/optional[5]) <  abs((params[5] - optional[2])/optional[3])){
          meanTime <- optional[4]
          sdTime <- optional[5]
        }
      }

      dTime <- (params[5] - meanTime) / sdTime
      penalty <- penalty + optional[1] * dTime*dTime
      p[2] <- optional[1] * dTime*dTime
    }

    if (params[6] < 0){
      fixed[6] <- 0
      penalty <- penalty + params[6]*params[6]
      p[3] <- params[6]*params[6]
    }
    
    fixed[7] <- max( 0.05, min( params[7], 0.5 ) )
    delta <- fixed[7] - params[7]
    penalty <- penalty + delta*delta
    p[4] <- delta*delta
    
    fixed[5] <- max( 0.1, min( params[5], (tMax - params[2]) + 0.4*fixed[7] ) )
    delta <- fixed[5] - params[5]
    penalty <- penalty + delta*delta
    p[5] <- delta*delta
  }
  
  if (nPar >= 10){
    if (!is.null(optional) & length(optional) >= 5){
      meanTime <- optional[4]
      sdTime <- optional[5]
      dTime <- (params[8] - meanTime) / sdTime
      penalty <- penalty + optional[1] * dTime*dTime
      p[6] <- optional[1] * dTime*dTime
    }
    
    if (params[9] < 0){
      fixed[9] <- 0
      penalty <- penalty + params[9]*params[9]
      p[7] <- params[9]*params[9]
    }
    
    fixed[10] <- max( 0.05, min( params[10], 0.5 ) )
    delta <- fixed[10] - params[10]
    penalty <- penalty + delta*delta
    p[8] <- delta*delta
    
    fixed[8] <- max( 0.1, min( params[8], (tMax - params[2]) + 0.4*fixed[10] ) )
    delta <- fixed[8] - params[8]
    penalty <- penalty + delta*delta
    p[9] <- delta*delta
    
    delta <- params[5] - params[8]
    if (delta < 0.1){
      delta <- (0.1 - delta) / 0.1
      penalty <- penalty + delta*delta
      p[10] <- delta*delta
    }
  }

  fit <- model2.Rebuild(data,data[1,2],fixed)
  residue <- data[,2] - fit
  
  if (debug){
    print(sum(residue*residue))
    print(p)
  }

  return(sum(residue*residue) + as.numeric(penalty))
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

model2.FindSegment <- function(ppg,beats,index){
  if (index < 1 | index > nrow(beats)){
    return(c(0,0))
  }
  
  nppg <- nrow(ppg)
  nbeat <- nrow(beats)

  w <- model2.FindPeak( ppg, beats[[index,lab.time]] )
  wLow <- w[1]
  w0 <- w[2]
  if (index < nbeat){
    w <- model2.FindPeak( ppg, beats[[index+1,lab.time]] )
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
  beatTime <- beat[[lab.time]]

  result <- 1:10
  result[1] <- min(data[,2])-0.326

  residue <- model2.Excess(xy[,2],yPrev,result[1])
  #plot(data[,1],residue)
  
  ## S peak
  peak.w <- which(data[,1] > beatTime-0.2 & data[,1] < beatTime+0.2)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]

  #print(beatTime)
  #plot(peak.t,peak.y)
  #stop()
  
  
  result[3] <- max(peak.y)
  result[2] <- peak.t[which(peak.y==result[3])]
  result[4] <- 0.25
  residue <- model2.SubtractExcessPeak(data[,1],residue,result[2:4])
  #lines(data[,1],residue)
  
  ## D peak
  peak.w <- which(data[,1] > beatTime+0.2 & data[,1] < beatTime+0.4)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]
  
  result[6] <- max(peak.y)
  result[5] <- peak.t[which(peak.y==result[6])]-result[2]
  result[7] <- 0.25
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,result[5:7])
  #lines(data[,1],residue)
  
  ## N peak
  peak.w <- which(data[,1] > beatTime+0. & data[,1] < beatTime+0.2)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]

  result[9]  <- max(peak.y)
  result[8]  <- peak.t[which(peak.y==result[9])]-result[2]
  result[10] <- 0.25
  
  residue <- model2.SubtractExcessPeak(data[,1],residue,result[8:10])
  #lines(data[,1],residue)

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

model2.Excess.Inv <- function(excess,offset,baseline){
  nX <- length(excess)
  if (nX == 0){
    print("Help")
  }
  
  result <- 1:nX * 0.0

  result[1] = excess[1] + (baseline + config.rate*(offset-baseline))
  
  for (j in 2:nX){
    result[j] = excess[j] + (baseline + config.rate*(result[j-1]-baseline))
  }
  
  return(result)
}

model2.Peak <- function(time,peakParams){
  temp <- 2*const.pi*(time - peakParams[1])/peakParams[3]
  temp[which(temp < -const.pi)] = -const.pi
  temp[which(temp >  const.pi)] =  const.pi
  result <- peakParams[2] * (0.5 * (1+cos(temp)))^2
  
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

model2.Rebuild <- function(xy,offset,params){
  excess <- 1:nrow(xy) * 0.0
  excess <- excess + model2.Peak(xy[,1],params[2:4])
  if (length(params)>=7){
    excess <- excess + model2.Peak(xy[,1],params[5:7]+c(params[2],0,0))
  }
  if (length(params)>=10){
    excess <- excess + model2.Peak(xy[,1],params[8:10]+c(params[2],0,0))
  }

  result <- model2.Excess.Inv(excess,offset,params[1])
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
    "Start Time (s)","Start PPG",
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
