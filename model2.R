## Constants
#

lab.time = "time (s)"
lab.ppg = "Detrended"
lab.dt = "interval (s)"

config.rate = 0.95

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

model2.ChiSq <- function(data,params,optional=NULL){
  residue <- model2.Residue(data,params)
  return(sum(residue*residue))
}

model2.ChiSqPPG <- function(data,params,optional=NULL){
  fit <- model2.Rebuild(data,data[1,2],params)
  residue <- data[,2] - fit
  return(sum(residue*residue))
}

model2.ChiSqPPGBound <- function(data,params,optional){
  penaltyScale = optional[1]
  meanDTime = optional[2]
  sdDTime = optional[3]
  meanNTime = optional[4]
  sdNTime = optional[5]
  
  dD = (params[5]-meanDTime)/sdDTime
  dN = (params[8]-meanNTime)/sdNTime
  negS = min(0,params[3])
  negD = min(0,params[6])
  negN = min(0,params[9])
  sWs = min(0.1,max(0.001,params[4]))
  sWd = min(0.1,max(0.001,params[7]))
  sWn = min(0.1,max(0.001,params[10]))
  
  par <- params
  par[3] <- max(0,params[3])
  par[6] <- max(0,params[6])
  par[9] <- max(0,params[9])
  
  penalty =
    penaltyScale * (dD*dD + dN*dN) +
    (negS*negS + negD*negD + negN*negN) +
    (1/sWs + 1/sWd + 1/sWn - 3/0.1)
  
  fit <- model2.Rebuild(data,data[1,2],par)
  residue <- data[,2] - fit
  return(penalty + sum(residue*residue))
}

model2.ChiSqZeroN <- function(data,params,optional){
  #print("0N")
  penaltyScale = optional[1]
  meanDTime = optional[2]
  sdDTime = optional[3]

  dD = (params[5]-meanDTime)/sdDTime
  negS = min(0,params[3])
  negD = min(0,params[6])
  sWs = min(0.1,max(0.001,params[4]))
  sWd = min(0.1,max(0.001,params[7]))

  ## Append zero-amplitude N peak
  par <- c(params,params[5],0,0.1)
  par[3] <- max(0,params[3])
  par[6] <- max(0,params[6])

  penalty =
    penaltyScale * (dD*dD) +
    (negS*negS + negD*negD) +
    (1/sWs + 1/sWd - 2/0.1)
  #print(paste("Penalty:",penaltyScale * (dD*dD),",",
  #              (negS*negS + negD*negD),",",
  #              (1/sWs + 1/sWd - 2/0.1)))
  
  fit <- model2.Rebuild(data,data[1,2],par)
  residue <- data[,2] - fit
  #print(paste("/0N:",penalty,sum(residue*residue)))
  return(penalty + sum(residue*residue))
}

model2.FindSegment <- function(ppg,beats,index){
  if (index < 1 | index > nrow(beats)){
    return(c(0,0))
  }
  
  tThis <- beat[[index,lab.time]]
  if (nrow(beat)>index){
    tNext <- beat[[index+1,lab.time]]
  }else{
    tNext <- tThis + 2.0
  }
  w <- which(ppg[[lab.time]]>tThis-0.125 & ppg[[lab.time]]<tNext-0.125)

  return(c(min(w),max(w)))
}

model2.GetSegment <- function(ppg,limits){
  w <- c(limits[1]:limits[2])
  
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
  peak.w <- which(data[,1] > beatTime-0.1 & data[,1] < beatTime+0.1)
  peak.t <- data[peak.w,1]
  peak.y <- residue[peak.w]

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
  excess <- excess + model2.Peak(xy[,1],params[5:7]+c(params[2],0,0))
  excess <- excess + model2.Peak(xy[,1],params[8:10]+c(params[2],0,0))
  
  result <- model2.Excess.Inv(excess,offset,params[1])
}

model2.Residue <- function(xy,params){
  result <- model2.Excess(xy[,2],params[1],params[1])
  result <- model2.SubtractExcessPeak(xy[,1],result,params[2:4])
  result <- model2.SubtractExcessPeak(xy[,1],result,params[5:7]+c(params[2],0,0))
  result <- model2.SubtractExcessPeak(xy[,1],result,params[8:10]+c(params[2],0,0))
  return(result)
}

model2.CleanParams <- function(params){
  result <- params
  
  if (mean(result[,7]-result[,4])>0){
    result[,7] <- result[,7] - result[,4]
    result[,10] <- result[,10] - result[,4]
  }
  
  result[,5] <- pmax(0,result[,5])
  result[,6] <- pmax(0.1,result[,6])
  result[,7] <- pmax(0,result[,7])
  result[,8] <- pmax(0,result[,8])
  result[,9] <- pmax(0.1,result[,9])
  result[,10] <- pmax(0,result[,10])
  result[,11] <- pmax(0,result[,11])
  result[,12] <- pmax(0.1,result[,12])
  
  return(result)
}

model2.Load <- function(filename){
  if (file.exists(filename)){
    params <- read.csv(filename)
    # Erase line index column
    params["X"] <- NULL
    
    result <- model2.CleanParmas(data.matrix(params))
    
    return(result)
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
