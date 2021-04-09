## Constants

lab.time = "time (s)"
lab.ppg = "Detrended"
lab.dt = "interval (s)"

config.rate = 0.85 
config.a.b = 0.999
config.optimum.width = 0.3
config.penalty.width = 10

const.pi = 3.1415926535897932384626433

## End of constants section

# Functions:
# 1. model2.LoadBeat
# 2. model2.LoadPPG
# 3. model2.FindPeak
# 4. model2.GetSegment
# 5. model2.Excess
# 6. model2.Peak
# 7. model2.Plot
# 8. model2.Residue
# 9. model2.Load
# 10. model2.Save
# 11. model2.SubtractExcessPeak

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


model2.GetSegment <- function(ppg,limits){
  
  #seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(ppg$`time (s)` == nextTime))
  
  w <- c(limits[1]:limits[3])

  result <- matrix(nrow=length(w),ncol=2)

  result[,1] <- ppg[[lab.time]][w]
  result[,2] <- ppg[[lab.ppg]][w]

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
