#Iso fitting funcs:

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

FactorAdjust <- function(ppg, beat, fs = model2.FindSegment, gs = model2.GetSegment, u = UnDetrend, factorCutoff = -20, plot = T){
  # Extract the first beat:
  beatTime <- beat[1,1]
  nextTime <- beat[2, 1] 
  seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(ppg$`time (s)` == nextTime))
  data <- gs(ppg,seg)
  if(plot == TRUE){plot(data)} 
  
  # Calculate the Gradient of the tail of the beat:
  tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2], data[nrow(data)-2, 2], data[nrow(data)-3, 2], 
            data[nrow(data)-4, 2])    
  tail <- rev(tail)
  xx. <- 1:length(tail)
  y. <- lm(tail~xx.)
  # Adjust the factor value until the gradient of the tail reaches an appropriate threshold:
  factor_value <- 1.01
  while(y.[[1]][2] > factorCutoff){    
    factor_value <- factor_value - 0.01
    ppg2 <- ppg
    ppg2[, 2] <- u(ppg,factor=factor_value,offset=1)
    beatTime <- beat[1,1]
    #a <- min(which(ppg[round(inflexX[o]), 1] > beatTime)) 
    nextTime <- beat[2, 1]    #ppg[round(inflexX[o])[a], 1]
    #nextTime <- if (i < nrow(beat)){beat[i+1,1]}else{NA}
    #seg <- fs(ppg,beat[i,1],nextTime)
    seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(ppg$`time (s)` == nextTime))
    data <- gs(ppg2,seg)
    if(plot == TRUE){plot(data)}
    if(y.[[1]][2] > 0){
      tail <- c(data[nrow(data)-5, 2], data[nrow(data)-6, 2], data[nrow(data)-7, 2], data[nrow(data)-8, 2], 
                data[nrow(data)-9, 2])
    }else{
      tail <- c(data[nrow(data), 2], data[nrow(data)-1, 2], data[nrow(data)-2, 2], data[nrow(data)-3, 2], 
                data[nrow(data)-4, 2])
    }
    tail <- rev(tail)
    xx. <- 1:length(tail)
    y. <- lm(tail~xx.)
  }
  return(factor_value)
}



OffsetAdjust <- function(ppg3, ppg, u = UnDetrend, factor_value, plot = T){
  vv. <- ppg3[, 1]      
  yv. <- lm(ppg3[, 2]~vv.)
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
    ppg3 <- data.frame(ppg[,1], u(ppg,factor=factor_value,offset=offset_value))
    vv. <- ppg3[, 1]      
    yv. <- lm(ppg3[, 2]~vv.)
    if(plot == TRUE){
      if(yv.[[1]][2]>0){plot(ppg[,1],u(ppg,factor=factor_value,offset=offset_value), type = "l")}
      if(yv.[[1]][2]>0){abline(a = yv.[[1]][1], b = yv.[[1]][2], col = "red")}
      if(yv.[[1]][2]>0){print(yv.[[1]][2])} 
    }
  }
  return(offset_value)
}

AddOutput <- function(beat){
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
  beat$config.rate= rep(0.75, nrow(beat))  
  return(beat)
}


FindStartParams <- function(batch_number, beats_in, beat, ppg, fs = model2.FindSegment, gs = model2.GetSegment, e = model2.Excess, sep = model2.SubtractExcessPeak, all_beats = FALSE, plot = FALSE){
  nBeats <- nrow(beat)
  seg <- c(0,0,0)
  if((batch_number*beats_in) > nBeats){
    print("Batch and beat values request more beats than are in time series, defaulting to max number of beats")
    maxn <- nBeats - 1
  }else{
    maxn <-(batch_number*beats_in)
    if(all_beats == TRUE){maxn <- maxn + (nrow(beat) - maxn)}
  }
  ob <- c()
  for(i in 1:(nrow(beat)-1)){
   ob[i] <- beat[i+1, 1] - beat[i, 1]    
  }
  ob_thrld <- mean(ob) + sd(ob)*4
  #segn <- 0
  for(i in 1:maxn){  
      
      beatTime <- beat[i,1]
      a <- min(which(ppg[round(inflexX[o]), 1] > beatTime))  
      nextTime <- ppg[round(inflexX[o])[a], 1]
      # If nextBeat is unreasonably late, shorten the segment:
       if(nextTime - beatTime >  ob_thrld){        
        nextTime <- beatTime + round(median(ob)) - 0.2      
       }
      seg <- c(which(ppg$`time (s)` ==  beatTime), 0, which(abs(ppg[, 1] - nextTime) == min(abs(ppg[, 1] - nextTime))))
      data <- gs(ppg,seg)
      if(plot == TRUE){plot(data)}

    if(nrow(data) < 10){next}
    tStart <- ppg[seg[1],1]
    yPrev <- ppg[max(seg[1]-1,1),2]
    
    amp <- max(data[, 2]) - min(data[, 2])
    constant <- 0.1092254*amp
    baseline <- min(data[,2]) - constant    
    residue <- e(data[,2], ppg[seg[1]-1,2], -0)     
    
    count <- nrow(data)
    excess <- 1:count * 0.0
    excess[1] = data[1,2] - (baseline + config.rate*(yPrev-baseline))
    for (j in 2:count){
      excess[j] = data[j,2] - (baseline + config.rate*(data[j-1,2]-baseline))
    }
    rm(count)
    rm(j)
    #plot(data[,1],excess, type = "l")
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
    rm(peak.w,peak.t,peak.y)
    residue <- sep(data[,1],residue,par[2:4])
    #plot(data[,1],excess)
    #lines(data[,1],residue)  
    
    # D peak
    peak.w <- which(data[,1] > beat[i,1]+0.3 & data[,1] < beat[i,1]+0.6)  # this finds a range in where to find the peak of d
    peak.t <- data[peak.w,1]     # this finds the time corresponding to the peak
    peak.y <- residue[peak.w]
    #plot(data[,1],excess)   
    #lines(peak.t,peak.y)
    par[6] <- max(peak.y)     
    par[5] <- peak.t[which(peak.y==par[6])]   
    par[7] <- 0.25  
    rm(peak.w,peak.t,peak.y)
    residue <- sep(data[,1],residue,par[5:7])
    #plot(data[,1],excess)
    #lines(data[,1],residue)
    
    # N peak
    t <- par[2] + c(0.25,0.75) * (par[5]-par[2])
    peak.w <- which(data[,1] > t[1] & data[,1] < t[2])
    peak.t <- data[peak.w,1]
    peak.y <- residue[peak.w]
    #plot(data[,1],excess)   
    #lines(peak.t,peak.y)
    par[9] <- max(peak.y)
    par[8] <- peak.t[which(peak.y==par[9])]
    par[10] <- 0.25
    rm(peak.w,peak.t,peak.y,t)
    residue <- sep(data[,1],residue,par[8:10])
    #plot(data[,1],excess)
    #lines(data[,1],residue)
    
    # Store parameters
    w <- seg[1]:seg[3]
    ppg$Baseline[w] <- baseline
    ppg$Excess[w] <- excess
    ppg$Residue[w] <- residue
    rm(w,excess,residue,data,baseline,yPrev,nextTime,tStart)
    beat[i,3:4]  = c(seg[1],seg[3])
    beat[i, 5] <- par[1]
    beat[i,6:15] = par
    beat[i,10] = beat[i,10]-beat[i,7]
    beat[i,13] = beat[i,13]-beat[i,7]
    rm(par)
  }
  rm(seg)
  temp <- list(beat, ppg)
  return(temp)
}


FindWithinParams <- function(beats_in, ppg, beat, gs = model2.GetSegment, fp = model2.FixParams3, ms = simplex.MakeSimplex3, m2 = model2.ChiSq){
  a <- list()
  for(i in 1:beats_in){         
    seg <- c(beat[i,3],0,beat[i,4])
    if(seg[1] == 0){next}
    data <- gs(ppg,seg)
    rm(seg)
    
    par <- as.numeric(beat[i,5:16])
    par <- fp(data[, 1:2], par)
    
    a[[i]] <- ms(data[,1:2],par,m2,0.1) 
  }
  return(a)
}

extractOutput <- function(beats_in, sim){
  across <- sim[1, ][1:6]
  within <- list()
  for(i in 1:beats_in){
    temp <- rep(0, 12)
    temp[c(1:4, 7, 10)] <-  sim[1, ][((i*6)+1):((i*6)+6)]  
    within[[i]] <- temp
  }
  temp <- list(across, within)
  return(temp)
}

FixOutput <- function(beats_in, beat, ppg, gs = model2.GetSegment, fp = model2.FixParams3, across = output[1], within = output[2]){
  fixed <- list()
  for(i in 1:beats_in){
    seg <- c(beat[i,3],0,beat[i,4])
    data <- model2.GetSegment(ppg,seg)
    rm(seg)
    fixed[[i]] <- model2.FixParams3(data, params = as.numeric(within[[i]]), across_beat_params = across)
  } 
  return(fixed)
}

UpdateBeat <- function(beats_in, beat, fixed){
  new_beat <- data.frame(matrix(0, ncol = 12, nrow = beats_in))
  for(i in 1:beats_in){
    new_beat[i, ] <- fixed[[i]]
  }
  new_beat <- cbind(beat[, 1:4], new_beat)
  return(new_beat)
}

PlotFits <- function(beats_in, ppg, beat2, gs = model2.GetSegment, rb = model2.Rebuild2){
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


model2.ChiSq3 <- function(data, params,debug=FALSE, beats, optional = NULL, beat, a = NULL, plot = FALSE, renal_param, dias_param){  
  
  # Across-beat parameter extraction:
  if(!is.null(a)){                                            # If a 66 parameter vector has been supplied, extract the first 6 
    across_beat_params <- a[1:6]
  }else{                                                      # If not, take them from the params input
    par <- params
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]
  }
  
  # Calculation of ChiSq for all beats:
  beat_fit <- list()
  failed_segs <- c()
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
    #dias <- across_beat_params[2] + par2[3]
    dias <- par2[3] + dias_param
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
    
    # Weighted region is W -> D (with slope)
    residue[w:end[1]] <-  residue[w:end[1]]*3
    if(length(residue) > end[1]){
      tail <- (end[1]+1):length(residue)
      for(j in 1:length(tail)){
        wgt <- 3 - (0.1*j)
        if(wgt < 1){wgt <- 1}
        residue[tail[j]] <- residue[tail[j]]*wgt
      }
    }
    
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
    if(sum(i == failed_segs) > 0){next}
    temp[i] <- beat_fit[[i]][1]
  }
  ts_fit <- sum(temp)
  return(ts_fit)
}

model2.ChiSq4 <- function(data, params,debug=FALSE, beats, optional = NULL, beat, a = NULL, plot = FALSE, renal_param, dias_param){  
  
  # Across-beat parameter extraction:
  if(!is.null(a)){                                       
    across_beat_params <- a[1:6]
  }else{                                           
    par <- params
    across_beat_params <- par[c(5, 6, 8, 9, 11, 12)]
  }
  
  # Calculation of ChiSq for all beats:
  beat_fit <- list()
  max_error <- list()
  failed_segs <- c()
  for(i in 1:beats[[1]]){                                    
    
    # Within-beat parameter extraction:
    if(!is.null(a)){                                            
      par2 <- a[((i*6)+1):((i*6)+6)]    
      par2 <- c(par2[1:4], 0, 0, par2[5], 0, 0, par2[6], 0, 0)
    }else{                                                        
      par2 <- as.numeric(beat[i, 5:16])
    }
    
    # Extract individual beat data:  
    seg <- c(beats[[2]][i],0,beats[[3]][i])
    dat <- model2.GetSegment(data,seg)
    rm(seg)
    
    # Extract systolic and diastolic parameters:
    sys <- par2[3]
    #dias <- across_beat_params[2] + par2[3]
    dias <- par2[3] + dias_param
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
    
    # Calculate fit, residue and max error:
    fit <- model2.Rebuild2(dat,dat[1,2],params = fixedPar)    
    residue <- dat[ ,2] - fit
    max_error[[i]] <- max(residue)
    
    # Weighted region is W -> D (with slope)
    residue[w:end[1]] <-  residue[w:end[1]]*3
    if(length(residue) > end[1]){
      tail <- (end[1]+1):length(residue)
      for(j in 1:length(tail)){
        wgt <- 3 - (0.1*j)
        if(wgt < 1){wgt <- 1}
        residue[tail[j]] <- residue[tail[j]]*wgt
      }
    }
    
    # Calculate Reduced Chi-Square for the beat:
    nData <- nrow(dat)    
    nPar <- length(par2)
    beat_fit[[i]] <- (sum(residue*residue) / (nData-nPar)) + as.numeric(penalty)
    
    if(plot == TRUE){
      plot(dat,  ylim = c(-150, 1600))     
      lines(dat[, 1], fit)
    }
  }
  
  # Summate individual beat ChiSq values:
  temp <- c()
  for(i in 1:length(beat_fit)){
    if(sum(i == failed_segs) > 0){next}
    temp[i] <- beat_fit[[i]][1]
  }
  ts_fit <- sum(temp)
  
  fit <- list(ts_fit, beat_fit, max_error)
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
        t[3] <- renal_param + 0.02
      }
    }
    
    if(i==3){                                              # Renal peak timing should be similar to initial estimation
      if(t[3] < (renal_param - 0.02) ){
        penalty <- penalty + 5000
        t[3] <- renal_param - 0.02
      }
    }
    
    if(i==1){                                               # S amplitude should not deviated, horizontally or vertically, from the max point of the data
      max.amp <- data[which.min(abs(data[, 1] - t[1])), 2]
      if(h[1] > max.amp + 50){
        penalty <- penalty + 100000
        h[1] <- max.amp
      }
      if(h[1] > max.amp - 50){
        penalty <- penalty + 100000
        h[1] <- max.amp
      }
      if(h[1] > max.amp + 100){
        penalty <- penalty + 200000
        h[1] <- max.amp
      }
      if(h[1] > max.amp - 100){
        penalty <- penalty + 200000
        h[1] <- max.amp
      }
    }
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
    t[1] <- y.
  }
  if(t[1] < y. - 0.04){
    penalty <- penalty + 10000
    t[1] <- y.
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
  
  # Don't clamp baseline shift (and not currently penalized)
  diff <- abs(baseline[1] - baseline[2])  
  #penalty <- penalty + META_BASELINE_SHIFT*diff  
  
  # Fix Config.rate 
  if(across_beat_params[6] > 0.95){
    penalty <- penalty + 1000
    across_beat_params[6] <- 0.95
  }

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

simplex.MakeSimplex2 <- function(data,param,f,inScale,directions=NULL,inTol=-1,optional=NULL,debug=FALSE, beat_vector = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param){
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
    tol <- min(1,f(data,param, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param))  
  }
  
  
  chiSq <- 1:(nPar+1) * 0.0
  chiSq[1] <- f(data,param,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param)   # ChiSq[1 is the fit when no parameters are changed... 
  if (debug){ print(paste("Root chi-squared:",chiSq[1]))}
  
  result <- matrix(nrow=nPar+1,ncol=nPar)
  result[1,] <- as.double(param)
  
  useDirections = !is.null(directions)
  if (useDirections){ useDirections <- nrow(directions) == nPar & ncol(directions) == nPar}
  
  
  for (i in c(5, 6, 8, 9, 11, 12)){   
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
    chiSqMinus <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param)            # tParam = test parameter. 
    tParam <- param + delta                                   # The chisquare (goodness of fit) is calculated for each direction, 
    chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param)            # and presumably the direction with the smaller value is chosen. 
    
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
    
    if (chiSq[i+1] < chiSq[1]){         # If the new fit is better than the old fit (with no parameters changed), continue to go in the direction that improved the fit
      if (debug){ print("Extending as best point") }
      while (chiSq[i+1] < chiSq[1] + tol){                 # Chisquare keeps getting iterated here (for 10 iterations)
        delta <- 2*delta    # 2* was too much here... WHY DOES IT SEEM LIKE SAMPLITUDE IS COMING DOWN???
        tParam <- param + delta
        oldScore <- chiSq[i+1]          # The current best fit gets called 'old score'
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param)   # The new fit is now designated ChiSq[i+1]
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
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param)
        if (debug){ print(paste("chi^2(",tParam[i],") =",chiSq[i+1])) }          # repeated code up until here... 
        if (chiSq[i+1] - oldScore < oldScore - chiSq[1]){   # Presumably again if the fit was worse than before, reverse it
          tParam <- param + 0.5*delta
          chiSq[i+1] <- oldScore
          break
        }
        iKill <- iKill - 1
        if (iKill < 0){
          if(i == 9){
            tParam[9] <- renal_param  # Ignore renal times that can't optomize
            break
          } 
          print(c("simplex constructed as per original parameter"))
          break
          #print("Failed to construct simplex")
          #return(paste("Error: param[",i,"]",sep=""))
        }
      }
    } else {
      if (debug){ print("Shrinking above tolerance") }
      while (chiSq[i+1] > chiSq[1] + tol){              # If the new fit is much worse than the original, reduce the size of delta
        delta <- 0.5*delta
        tParam <- param + delta
        lastChiSq <- chiSq[i+1]
        chiSq[i+1] <- f(data,tParam,optional=optional, beats = beat_vector, beat = beat, renal_param = renal_param, dias_param = dias_param)
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
          #return(paste("Error: param[",i,"]",sep=""))
          print(c("Failed to construct simplex within 10 iterations for parameter", i, "defaulting to inputted value"))
          tParam[i] <- param[i]
          next
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





simplex.Run2 <- function(data = ppg,simplexParam = sim,f = model2.ChiSq3,optional=NULL, beat_vector = beat_vector, renal_param = renal_param, dias_param = dias_param, run = NULL){
  
  MAX_STEP <- 500                                               # The number of steps to iterate through
  FTOL <- 1e-5                                  
  
  debugRtol <- 1:(MAX_STEP+1) * 0.0
  debugMin <- 1:(MAX_STEP+1) * 0.0
  debugMax <- 1:(MAX_STEP+1) * 0.0
  
  result <- simplexParam                         # Now feed in the 66*66 matrix
  nPar <- ncol(result)                           
  chiSq <- 0:nPar * 0.0
  for (i in 1:(nPar+1)){                                    # Find out the ChiSq value for each row from result
    chiSq[i] <- f(data, params = NULL, optional=NULL, a = result[i, ], beats = beat_vector, renal_param = renal_param, dias_param = dias_param)
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
    score <- f(data, params = rep(0, 12),optional=optional, a = test, beats = beat_vector, renal_param = renal_param, dias_param = dias_param)
    
    if (score < chiSqMin){                          # If flipping improves the ChiSq, try extending further in the same direction
      test2 <- node - 2 * (apex - node)
      score2 <- f(data, params = rep(0, 12),optional=optional, a = test2, beats = beat_vector, renal_param = renal_param, dias_param = dias_param)
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
      score2 <- f(data, params = rep(0, 12),optional=optional, a = test2, beats = beat_vector, renal_param = renal_param, dias_param = dias_param)
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
            chiSq[i] <- f(data, params = rep(0, 12),optional=optional, a = result[i, ], beats = beat_vector, renal_param = renal_param, dias_param = dias_param)
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


osnd_fit <- function(bf = beat_final, ppg, gs = model2.GetSegment, r = model2.Rebuild2, sf = splinefun, dp = diast_pk, oa = osnd_of_average, sr = samplingRate, plot = FALSE){
  
  osnd_diff <- list()
  for(i in 1:nrow(bf)){  
    # Find the correct data segments and corresponding model fit:
    seg <- c(bf[i,3],0,bf[i,4])  
    data <- gs(ppg,seg)
    yPrev <- ppg[seg[1]-1,2]
    xPrev <- ppg[seg[1]-1, 1]
    xNext <- ppg[seg[3], 1]
    rm(seg)
    temp <- r(data, yPrev, as.double(bf[i,-c(1:4)]),TRUE)   
    
    # Upsample fit - so that fit OSND can be calculated as the data would have been in the main script
    sfunction <- sf(1:length(temp), temp, method = "natural")
    fit <-  sfunction(seq(1, length(temp), 0.1), deriv = 0)
    # Upsample data segment:
    sfunction <- sf(1:length(data[, 2]), data[, 2], method = "natural")
    dat <-  sfunction(seq(1, length(data[, 2]), 0.1), deriv = 0)
    
    # Find OSND of fit:
    tmp <- dp(avw = fit, sr = sr, scale = T)  
    dPeak <- tmp[1]
    xShift <- tmp[2]
    rm(tmp)
    osnd_fit <- oa(fit, dp = dPeak, diff = 0, sr = sr, plot = F)
    # Find OSND of data:
    tmp <- diast_pk(avw = dat, sr = sr, scale = T)  
    dPeak <- tmp[1]
    xShift <- tmp[2]
    rm(tmp)
    osnd_dat <- oa(dat, dp = dPeak, diff = 0, sr = sr, plot = F)
    
    # Adjust x-axis for upsampling and sampling rate, then find the difference between OSND points:
    osnd_dat$x <- osnd_dat$x/(10*sr)
    osnd_fit$x <- osnd_fit$x/(10*sr)
    osnd_diff[[i]] <- osnd_dat - osnd_fit
    
    if(plot == TRUE){
      # Plot them together:
      plot((1:length(dat))/(10*samplingRate), dat, type = "l", xlab = "time", ylab = "")
      lines((1:length(dat))/(10*samplingRate), fit, col = "red")
      points(osnd_dat, pch = 19)
      points(osnd_fit, col = "red", pch = 19)
    }
  }
  
  return(osnd_diff)
}
