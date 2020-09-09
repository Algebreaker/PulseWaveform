source("model2.R")
source("simplex.R")

CalculateParams <- function(input){
  ppg <- model2.LoadPPG(input)
  beat <- model2.LoadBeat(input)
  nBeats = nrow(beat)

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
      params <- matrix(nrow=nBeats,ncol=14)
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

      beatTime <- beat[i,1]
      nextTime <- if (i < nBeats){beat[i+1,1]}else{NA}

      if (seg[1] == 0){
        seg <- model2.FindSegment(ppg,beatTime,nextTime)
        all <- c(seg[1],seg[3])
      } else {
        tSeg <- model2.FindSegment(ppg,beatTime,nextTime)
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
      params[i,2] <- seg[1]
      params[i,3] <- yPrev
      params[i,4:13] <- par
      params[i,14] <- model2.ChiSq(data,par)

      wSeg <- seg[1]:seg[3]
      fit[wSeg] <- model2.Rebuild(data,yPrev,par)
      residue[wSeg] <- ppg[wSeg,2] - fit[wSeg]
      baseline[wSeg] <- params[i,4]
    }

    print("Completed initial parameter generation")
    model2.Save(paste(input,"_params_a.csv",sep=""),params)
  } else {
    print("Using estimated fitting parameters...")
  }

  # Limit variability in timing of D and N peaks
  fixedParams = c(1.0,median(params[,8]),sd(params[,8]),median(params[,11]),sd(params[,11]))

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
    beatTime <- beat[i,1]
    nextTime <- if (i < nBeats){beat[i+1,1]}else{NA}

    if (seg[1] == 0){
      seg <- model2.FindSegment(ppg,beatTime,nextTime)
      all <- c(seg[1],seg[3])
    } else {
      tSeg <- model2.FindSegment(ppg,beatTime,nextTime)
      seg[1] <- seg[3] + 1
      seg[2] <- tSeg[2]
      seg[3] <- tSeg[3]
      all[2] <- tSeg[3]
    }
    data <- model2.GetSegment(ppg,seg)

    tStart <- ppg[seg[1],1]
    yPrev <- ppg[max(seg[1]-1,1),2]

    tLim[1] <- min(tLim[1],tStart)
    tLim[2] <- max(tLim[2],ppg[seg[3],1])

    count <- seg[3] + 1 - seg[1]
    data[,1] <- ppg[seg[1],1] + (1:count - 1) / 75.0

    par <- params[i,4:13]
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
    params[i,2] <- seg[1]
    params[i,3] <- yPrev
    params[i,4:13] <- best
    params[i,14] <- chiSq2

    if (chiSq2 > 2){
      print(paste("Beat",i,"has chi-squared of",chiSq2))
    }

    wSeg <- seg[1]:seg[3]
    fit[wSeg] <- model2.Rebuild(data,yPrev,best)
    residue[wSeg] <- ppg[wSeg,2] - fit[wSeg]
    baseline[wSeg] <- params[i,4]
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
}

CalculateFixedPeakParams <- function(input){
  ppg <- model2.LoadPPG(input)
  params <- model2.Load(paste(input,"_params_b.csv",sep=""))
  
  if (is.null(params)){
    print("No fitting parameters found")
    return()
  }

  count <- nrow(ppg)
  nBeats = nrow(params)

  all <- c(params[1,2],0)
  fit <- 1:count * 0
  residue <- 1:count * 0
  baseline <- 1:count * 0

  # Use median to avoid outliers, which might be very poorly distributed
  fixedParams <- c(
    median(params[,7]),  # S width
    median(params[,8]),  # D timing
    median(params[,10]),  # D width
    median(params[,11]), # N timing
    median(params[,13])  # N width
  )
  medianDAmplitude <- median(params[,9])

  print("Fixing parameters:")
  print(paste("  S width =",fixedParams[1]))
  print(paste("  D time  =",fixedParams[2]))
  print(paste("  D width =",fixedParams[3]))
  print(paste("  N time  =",fixedParams[4]))
  print(paste("  N width =",fixedParams[5]))

  print("Refining fitting parameters with fixed peak width...")
  startClock <- Sys.time()
  for (i in 1:nBeats){
    if (as.double(difftime(Sys.time(),startClock,unit="secs"))>reportClock){
      min <- as.integer(reportClock/60)
      sec <- reportClock %% 60
      time <- paste("[",min,":",sep="")
      if (sec == 0){ time <- paste(time,"00]",sep="")} else { time <- paste(time,sec,"]",sep="")}
      print(paste(time,"Processing beat",i))
      reportClock <- reportClock + 10
    }

    beatTime <- params[i,5]
    nextTime <- if (i < nBeats){params[i+1,1]}else{NA}

    tStart <- params[i,1]
    seg <- c(as.numeric(params[i,2]),0,0)
    yPrev <- params[i,3]

    if (i < nBeats){
      seg[3] <- params[i+1,2]-1
    } else {
      tSeg <- model2.FindSegment(ppg,beatTime,nextTime)
      seg[3] <- tSeg[3]
    }
    all[2] <- seg[3]

    count <- seg[3] + 1 - seg[1]
    data <- model2.GetSegment(ppg,seg)
    
    # Smooth times
    data[,1] <- ppg[seg[1],1] + (1:count - 1) / 75.0

    par <- params[i,c(4:6,9,12)]
    if (params[i,8] < 0.5 * (fixedParams[2] + fixedParams[4]))
    {
      par[4] <- params[i,9]
      par[3] <- 0
    }

    # Clamp initial amplitudes to be positive
    par[3:5] <- pmax(0.,par[3:5])

    sim <- simplex.MakeSimplex(data,par,model2.ChiSqAmp,0.1,optional=fixedParams)
    sim <- simplex.Run(data,sim,model2.ChiSqAmp,optional=fixedParams)
    #best <- sim[1,]
    #model2.ChiSqAmp(data,sim[1,],optional=fixedParams,debug=TRUE)
    best <- c(sim[1,1:3],fixedParams[1:2],sim[1,4],fixedParams[3:4])
    if (nPar >= 5){
      best <- c(best,sim[1,5],fixedParams[5])
    } else {
      best <- c(best,0,fixedParams[5])
    }
    #model2.ChiSq(data,best,debug=TRUE)
    chiSq2 <- model2.ChiSq(data,best)
    #print(chiSq2)

    # v Probably unchanged?
    params[i,1] <- tStart
    params[i,2] <- seg[1]
    params[i,3] <- yPrev
    # ^
    params[i,4:13] <- best
    params[i,14] <- chiSq2

    if (chiSq2 > 2){
      print(paste("Beat",i,"has chi-squared of",chiSq2))
    }

    wSeg <- seg[1]:seg[3]
    fit[wSeg] <- model2.Rebuild(data,params[i,3],best)
    residue[wSeg] <- ppg[wSeg,2] - fit[wSeg]
    baseline[wSeg] <- params[i,4]
  }
  print("Completed parameter generation")

  model2.Save(paste(input,"_params_c.csv",sep=""),params)

  # Make 'smooth' time axis
  time <- ppg[1,1] + (1:nrow(ppg) - 1) / 75
  w <- all[1]:all[2]
  plot(time[w],ppg[w,2],xlab="time (s)",ylab="PPG count",ylim=c(75,85))
  lines(time[w],fit[w])
  lines(time[w],baseline[w],lty="dashed")
  lines(time[w],max(ppg[,2])+residue[w],lty="dotted")
}


PlotFit <- function(input,ver=-1,pch=-1,xlim=NULL,xoff=0,ylim=NULL,yres=NULL){
  ppg <- model2.LoadPPG(input)
  
  filename <- ""

  vLim <- c(ver,ver)
  vLim <- if (ver < 0){c(1,5)} else {c(ver,ver)}

  for (c in unlist(strsplit("abcde",""))[vLim[1]:vLim[2]]){
    test <- paste(input,"_params_",c,".csv",sep="")
    if (file.exists(test)){
      filename <- test
    }
  }

  if (nchar(filename) == 0){
    print("No parameter files found")
    return()
  }

  params <- model2.Load(filename)
  nBeats <- nrow(params)

  count <- nrow(ppg)
  fit <- 1:count * 0
  excess <- 1:count * 0
  residue <- 1:count * 0
  baseline <- 1:count * 0

  all <- c(params[1,2],nrow(ppg))

  for (i in 1:nBeats){
    seg <- c(params[i,2],0,nrow(ppg))
    if (i < nBeats){
      seg[3] <- params[i+1,2]-1
    }

    dataseg <- model2.GetSegment(ppg,seg)

    tStart <- ppg[seg[1],1]

    yPrev <- ppg[max(seg[1],1),2]
    if (seg[1] > 1){
      yPrev <- ppg[max(seg[1]-1,1),2]
    }

    count <- seg[3] + 1 - seg[1]
    dataseg[,1] <- ppg[seg[1],1] + (1:count - 1) / 75.0

    tStart <- params[i,1]
    seg[1] <- params[i,2]
    yPrev <- params[i,3]
    par <- params[i,4:13]
    #chiSq2 <- params[i,14]
    chiSq2 <- model2.ChiSq(dataseg,par)

    wSeg <- seg[1]:seg[3]
    excess[wSeg] <- model2.Rebuild(dataseg,yPrev,par,invert=FALSE)
    fit[wSeg] <- model2.Excess.Inv(excess[wSeg],yPrev,par[1])
    residue[wSeg] <- ppg[wSeg,2] - fit[wSeg]
    baseline[wSeg] <- params[i,4]
  }

  # Make 'smooth' time axis
  time <- ppg[1,1] + (1:nrow(ppg) - 1) / 75 - xoff
  w <- all[1]:all[2]

  xlim1 <- NULL
  xlim2 <- NULL
  ylim1 <- NULL
  ylim2 <- NULL
  nHorizontal <- 1
  xTckLen <- 0.04

  if (!is.null(xlim)){
    if (length(xlim) >= 2){
      xlim1 <- xlim[1:2]
    }
    if (length(xlim) >= 4){
      xlim2 <- xlim[3:4]
      nHorizontal <- 2
      xTckLen <- 0.03
    }
  }

  if (!is.null(ylim)){
    if (length(ylim) >= 2){
      ylim1 <- ylim[1:2]
    }
    if (length(ylim) >= 4){
      ylim2 <- ylim[3:4]
    }
  }

  if (!is.null(xlim)){
    xrange <- xlim + 0.02 * c(-1,1) * (xlim[2]-xlim[1])
    w <- which(time > xrange[1] & time < xrange[2])
  }

  XMARGINS <- 1:10*0
  XMARGINS[1] <- 4.1
  XMARGINS[nHorizontal+1] <- 2.1

  for (iH in 1:nHorizontal){
    XLIM <- if(!is.null(xlim) && length(xlim) >= 2 * iH){xlim[(2*iH-1):(2*iH)]}else{NULL}
    if (!is.null(XLIM)){
      xrange <- XLIM + 0.02 * c(-1,1) * (XLIM[2]-XLIM[1])
      w <- which(time > xrange[1] & time < xrange[2])
    }

    XMAR <- XMARGINS[iH:(iH+1)]

    xBounds <- c(0,1)
    if (nHorizontal > 1){
      xBounds <- c(0,0.7,1)[iH:(iH+1)]
    }

    par(fig=c(xBounds,0.4,1),mar=c(0, XMAR[1], 4.1, XMAR[2]),new=(iH>1))

    plot(
      time[w],
      fit[w],
      axes=FALSE,
      frame=TRUE,
      t="l",
      xaxt='n',xlim=XLIM,xlab="time (s)",
      ylab=if(iH==1){"PPG count"},ylim=ylim1
      )

    axis(1,labels=FALSE,  lwd=-1,lwd.ticks=1,tck=0.02)
    axis(2,labels=(iH==1),lwd=-1,lwd.ticks=1,tck=xTckLen,las=2)
    axis(3,labels=FALSE,  lwd=-1,lwd.ticks=1,tck=0.02)
    axis(4,labels=FALSE,  lwd=-1,lwd.ticks=1,tck=xTckLen)

    if (iH == 2){
      points(time[w],ppg[w,2])
    } else if (pch != -1){
      points(time[w],ppg[w,2],pch=pch)
    }

    lines(time[w],baseline[w],lty="dashed")

    offset <- 0
    if (is.null(yres)){
      offset <- max(ppg[,2])
    } else {
      offset <- yres
    }
    lines(time[w],offset+residue[w],lty="dotted")

    par(fig=c(xBounds,0,0.4),mar=c(5.1, XMAR[1], 0, XMAR[2]), new=TRUE)

    plot(
      time[w],
      75 * excess[w],
      axes=FALSE,
      frame=TRUE,
      t="l",
      xlim=XLIM,xlab="time (s)",
      ylab="Excess gradient",ylim=ylim2
      )

    axis(1,labels=TRUE, lwd=-1,lwd.ticks=1,tck=0.04)
    axis(2,labels=(iH==1), lwd=-1,lwd.ticks=1,tck=2*xTckLen,las=2)
    axis(3,labels=FALSE,lwd=-1,lwd.ticks=1,tck=0.04)
    axis(4,labels=FALSE,lwd=-1,lwd.ticks=1,tck=2*xTckLen)
  }

  par(fig=c(0,1,0,1),mar=c(5.1, 4.1, 4.1, 2.1))
}

#for (c in unlist(strsplit("ABCDEFGH",""))){
#  CalculateParams(paste("anon/",c))
#}
#for (c in unlist(strsplit("ABCDEFGH",""))){
#  CalculateFixedPeakParams(paste("anon/",c,sep=""))
#}

ver <- 3

if (ver == 2){
  pdf("Model2_%03d.pdf",width=12,height=8,onefile=FALSE)
  PlotFit("anon/A",ver=ver,xlim=c(0,30),xoff=620,ylim=c(78.5,82.7,-3,30,-2,32),yres=78.8)
  PlotFit("anon/B",ver=ver,xlim=c(0,30),xoff=10,ylim=c(76.5,82.5,-2,42),yres=77.5)
  PlotFit("anon/C",ver=ver,xlim=c(0,30),xoff=10,ylim=c(70,88,-2,162),yres=71)
  PlotFit("anon/D",ver=ver,xlim=c(0,30),xoff=20,ylim=c(74,85,-2,102),yres=74.5)
  PlotFit("anon/E",ver=ver,xlim=c(0,30),xoff=10,ylim=c(77.7,81.5,-2,27),yres=78.2)
  PlotFit("anon/F",ver=ver,xlim=c(0,30),xoff=2580,ylim=c(83,88.3,-2,42),yres=83.3)
  PlotFit("anon/G",ver=ver,xlim=c(0,30),xoff=300,ylim=c(77.4,81,-2,27.5),yres=77.5)
  PlotFit("anon/H",ver=ver,xlim=c(0,30),xoff=10,ylim=c(77.8,82,-2,45),yres=78)
  dev.off()
}

if (ver == 3){
  dt <- c(0,1.4)
  pdf("Model2_%03d.pdf",width=12,height=8,onefile=FALSE)
  PlotFit("anon/A",ver=ver,xlim=c(0,30,5.45+dt),xoff=620,ylim=c(78.5,82.7,-3,30,-2,32),yres=78.8)
  PlotFit("anon/B",ver=ver,xlim=c(0,30,2.25+dt),xoff=10,ylim=c(76.5,82.5,-2,42),yres=77.5)
  PlotFit("anon/C",ver=ver,xlim=c(0,30,2.5+dt),xoff=10,ylim=c(70,88,-2,162),yres=71)
  PlotFit("anon/D",ver=ver,xlim=c(0,30,2.3+dt),xoff=20,ylim=c(74,85,-2,102),yres=74.5)
  PlotFit("anon/E",ver=ver,xlim=c(0,30,1.3+dt),xoff=10,ylim=c(77.7,81.5,-2,27),yres=78.2)
  PlotFit("anon/F",ver=ver,xlim=c(0,30,1.35+dt),xoff=2580,ylim=c(83,88.3,-2,42),yres=83.3)
  PlotFit("anon/G",ver=ver,xlim=c(0,30,1.1+dt),xoff=300,ylim=c(77.4,81,-2,27.5),yres=77.5)
  PlotFit("anon/H",ver=ver,xlim=c(0,30,1.8+dt),xoff=10,ylim=c(77.8,82,-2,45),yres=78)
  dev.off()
}

