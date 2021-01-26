simplex.MakeSimplex <- function(data,param,f,inScale,directions=NULL,inTol=-1,debug=FALSE){
  if (debug){print("MakeSimplex -- debug")}
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
    tol <- min(1,f(data,param))
  }

  chiSq <- 1:(nPar+1) * 0.0
  chiSq[1] <- f(data,param)
  if (debug){ print(paste("Root chi-squared:",chiSq[1]))}

  result <- matrix(nrow=nPar+1,ncol=nPar)
  result[1,] <- as.double(param)

  useDirections = !is.null(directions)
  if (useDirections){ useDirections <- nrow(directions) == nPar & ncol(directions) == nPar }

  for (i in 1:nPar){
    if (debug){ print(paste("Parameter",i)) }
    tParam <- param

    # Pick a direction
    delta <- 1:nPar * 0
    if (useDirections){
      delta <- scale[i] * directions[i,]
    } else {
      delta[i] <- scale[i]
    }

    tParam <- param - delta
    chiSqMinus <- f(data,tParam)
    tParam <- param + delta
    chiSq[i+1] <- f(data,tParam)

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
      while (chiSq[i+1] < chiSq[1] + tol){
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
          return(paste("Error: param[",i,"]",sep=""))
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
          #print("Failed to construct simplex")
          return(paste("Error: param[",i,"]",sep=""))
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

simplex.Run <- function(data,simplexParam,f){
  MAX_STEP <- 5000
  FTOL <- 1e-5

  debugRtol <- 1:(MAX_STEP+1) * 0.0
  debugMin <- 1:(MAX_STEP+1) * 0.0
  debugMax <- 1:(MAX_STEP+1) * 0.0

  result <- simplexParam
  nPar <- ncol(result)
  chiSq <- 0:nPar * 0.0
  for (i in 1:(nPar+1)){
    chiSq[i] <- f(data,result[i,])
  }

  for (iStep in 1:MAX_STEP){
    extrema <- simplex.SortHighLow(chiSq)
    low <- extrema[1]
    nHigh <- extrema[2]
    high <- extrema[3]

    chiSqMax <- chiSq[high]
    chiSqMin <- chiSq[low]

    #print(paste("chi^2_min =",chiSqMin))
    #print(paste("argMax = ",high,"[",chiSqMax,"]",sep=""))

    rtol <- 2 * (chiSqMax - chiSqMin)/(chiSqMax + chiSqMin + 1e-10)
    if (rtol < FTOL){
      bestParam <- result[low,]
      result[low,] <- result[1,]
      result[1,] <- bestParam
      return(result)
    }
    debugRtol[iStep] <- rtol
    debugMin[iStep] <- chiSqMin
    debugMax[iStep] <- chiSqMax
    #print(paste("Continuing: rtol = ",rtol,", min = ",chiSqMin,", max = ",chiSqMax,sep=""))

    factor <- -1
    node <- simplex.HypoCentre(result,high)
    apex <- result[high,]
    test <- node - (apex - node)
    score <- f(data,test)

    if (score < chiSqMin){
      test2 <- node - 2 * (apex - node)
      score2 <- f(data,test2)
      if (score2 >= score){
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
    } else if (score >= chiSq[nHigh]) {
      # Test for shrink with optional reflection
      factor <- 0.5
      if (score < chiSqMax)
      {
        factor <- -0.5
      }
      test2 <- node + factor * (apex - node)
      score2 <- f(data,test2)
      
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
            chiSq[i] <- f(data,result[i,])
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

## Internal functions

simplex.HypoCentre <- function(mat_Param,index){
  nPar <- ncol(mat_Param)

  result <- 1:nPar * 0.0
  for (i in 1:(nPar+1)){
    if (i != index){
      result <- result + mat_Param[i,]
    }
  }
  return( result / nPar )
}

simplex.SortHighLow <- function(vec_ChiSq){
  nPar <- length(vec_ChiSq)

  low <- 1
  high <- 1
  nHigh <- 2
  if (vec_ChiSq[2] > vec_ChiSq[1]){
    high <- 2
    nHigh <- 1
  }

  for (i in 2:nPar){
    if (vec_ChiSq[i] < vec_ChiSq[low]){
      low <- i
    }
    if (vec_ChiSq[i] > vec_ChiSq[high]){
      nHigh <- high
      high <- i
    }
    if (i != high & vec_ChiSq[i] > vec_ChiSq[nHigh]){
      nHigh <- i
    }
  }

  return(c(low,nHigh,high))  
}
