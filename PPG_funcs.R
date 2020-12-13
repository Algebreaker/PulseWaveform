

baseline <- function(inx, iny, o, dat, sp, plot = FALSE){
  # Making a (non-polynomial) spline to fit the baseline
  sfunction2 <- splinefun(inx[o], iny[o], method = "natural")
  splineBase <- sfunction2(seq(1, length(dat)), deriv = 0)
  
  # Plotting spline_base on spline_poly
  if(plot){
    plot(sp)
    points(inx[o], iny[o], pch = 19)
    lines(splineBase)
  }
  
  # Correcting for baseline:
  baseCor <- dat - splineBase
  if(plot){
    plot(baseCor, type = "l")
    # Plot new baseline (y = 0)
    lines(1:length(dat), seq(from = 0, to = 0, length.out = length(dat)))
  }
  
  return(baseCor)
}



clean_wuv <- function(wuv, sp, inx, o){
  
  # Remove u values that are implausibly far from baseline
  falseU <- which(wuv$uY > (median(wuv$uY) + sd(wuv$uY)) | wuv$uY < -50)
  if(length(falseU) > 1){
    wuv <- wuv[-falseU, ]
  }
  
  # Find v-u differences:
  diffVU <- wuv$vY - wuv$uY
  
  # undetected artefacts can lead to repeats in half heights (if the deriv peak is too wide), which will lead to 
  # scale factors of 0, check for these:
  if(length(which(diffVU == 0)) > 0){
    dup <- which(diffVU == 0)
    wuv <- wuv[-dup, ]
    diffVU <- diffVU[-dup]
  }
  
  # Make a vector of waves with abnormally small scale factors and remove them:
  falseScale <- which(diffVU < (median(diffVU) - 4*(sd(diffVU))))
  if(length(falseScale) > 1){
    wuv <- wuv[-falseScale, ]
    diffVU <- diffVU[-falseScale]
  }
  
  ## Find o points again:
  oY <- predict(sp, inx[o])
  #plot(splinePolyBC)
  #points(inflexX[o], o_yval, pch = 19)
  
  # Find o-w difference:
  owDiff <- c()
  for(i in 1:length(wuv$wX)){
    owDiff[i] <- wuv$wX[i] - inx[o[i]]
  }
  
  # Find distance between o_points:
  oDiff <- c()
  for(i in 1:(length(inx[o])-1)){
    oDiff[i] <- inx[o[i+1]] - inx[o[i]]
  }
  
  # Find distance between W points (Inter-beat interval):
  ibi <- c()
  for(i in 1:(length(wuv$wX))-1){
    ibi[i] <- wuv$wX[i+1] - wuv$wX[i]
  }
  
  # Remove last W:
  wuv <- wuv[-nrow(wuv), ]
  
  # Make a vector of abnormal IBIs (the waves at the end of a sequence / before an artefact): 
  #plot(ibi)
  #points(which(ibi > 1.3*median(ibi)), ibi[which(ibi > 1.3*median(ibi))], pch = 19, col = "red")
  end_waves <- which(ibi > 1.3*median(ibi))
  
  # Remove end_waves from all vectors:
  if(length(end_waves) > 1){
    wuv <- wuv[-end_waves, ]
    diffVU <- diffVU[-length(diffVU)]
    diffVU <- diffVU[-end_waves]
  }
  d <- cbind(wuv, diffVU)
  dat <- list(d, ibi, oDiff)
  return(dat)
}



diast_pk <- function(avw, sr){
  # Find the diastolic peak on the average wave to inform OSND finding (also some adjusment of x-values for removal of NA values):
  avw <- avw[!is.na(avw)]
  
  # Need to find new W position (0.5) after removing NAs
  xShift <- which(abs(avw-0.5) == min(abs(avw - 0.5)))
  avWavePoly <- CubicInterpSplineAsPiecePoly(1:length(avw), avw, "natural")
  avInflexX <- solve(avWavePoly, b = 0, deriv = 1)
  avInflexY <- predict(avWavePoly, avInflexX)
  
  # Specify limitations for where the diastolic peak can first be found i.e between 120:230 on x-axis, and below 1 on y-axis:
  peaks <- order(avInflexY[which(avInflexX < 215 & avInflexX > 120 & avInflexY < 1)], decreasing = TRUE)
  diastPk <- avInflexX[which(avInflexX < 215 & avInflexX > 120 & avInflexY < 1)][peaks[1]]
  
  # diastPk will be NA for class 3 waveforms, in which case set a default value
  if(is.na(diastPk) | diastPk < avInflexX[peaks[1]]){
    diastPk <- 10*sr
  }
  return(c(diastPk, xShift))
}



find_average <- function(p, ao){
  
  # Find the last 10 values of each wave:
  last10 <- list()
  for(i in 2:(ncol(p))){
    last10[[i-1]] <- p[, i][!is.na(p[, i])][(length(p[, i][!is.na(p[, i])])-10):(length(p[, i][!is.na(p[, i])]))]
  }
  
  # Find the average last 10 values
  mean_last10 <- c()
  for(i in 1:10){   
    row_vector <- c()
    for(j in 1:length(last10)){      
      row_vector[j] <- last10[[c(j, i)]] 
    }
    mean_last10[i] <- mean(row_vector[!is.na(row_vector)])   
  }
  
  # Find the consecutive y-axis differences (gradient essentially) of last 10 values:
  mean_diff_last10 <- c()
  for(i in 1:9){                    # 9 here since 10 values will give 9 values for differences between them
    mean_diff_last10[i] <- mean_last10[i+1] - mean_last10[i]
  }
  
  # Simply calculating the average wave by averaging each row of the dataframe doesn't work at the end of the wave, 
  # since as waves end, they no longer factor in to the calculation of the average. The average would be erroneously 
  # drawn out until the end of the last wave. Therefore, a clone dataframe of p is used to continue the trajectories
  # of waves after they have ended (by using the average gradient above) and thus make the end of the average wave
  # a more accurate approximation of the true average. 
  
  # Replace NA values with continuing downward gradient in a clone dataframe:
  p_for_finding_average <- p
  for(i in 2:(ncol(p_for_finding_average))){
    for(j in 1:length(p_for_finding_average[, i][ao[[(i-1)]][-1]])){
      p_for_finding_average[, i][ao[[(i-1)]][-1]][j] <- p_for_finding_average[, i][ao[[(i-1)]][1]] + j*mean(mean_diff_last10[1:5])
    }
  }
  
  # Calculate the average wave by row:
  average_wave <- c()
  sd_wave <- c()
  median_wave <- c()
  for(i in 1:nrow(p_for_finding_average)){
    row_vector <- c()
    for(j in 2:(ncol(p_for_finding_average))){
      row_vector[j-1] <- p_for_finding_average[i, j] 
    }
    average_wave[i] <- mean(row_vector[!is.na(row_vector)])   
    sd_wave[i] <- sd(row_vector[!is.na(row_vector)]) 
    median_wave[i] <- median(row_vector[!is.na(row_vector)])
  }
  
  # Find where the end of the average wave should end: 
  # This is done by finding the average (mode) y-value of the final value of each wave (before trajectory continuation)
  end <- c()
  for(i in 2:ncol(p)){
    end[i-1] <- p_for_finding_average[, i][ao[[(i-1)]][1]]
  }
  # define mode function
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  # Made this into a for loop because there can be waves that have no values after o, and if all waves are like that, end will be null
  if(length(end) > 1){     
    average_end <- mode(round(end[!is.na(end)], digits = 2))
    # Remove elements of average wave after where its end should be:
    below_o <- which(average_wave < average_end)
    # Find the first point in the last quarter of the wave that goes below baseline
    b. <- below_o[min(which(below_o > (length(average_wave)/(4/3))))]
    # If the wave doesn't go below baseline, don't remove any 
    if(is.na(b.) == FALSE){
      # otherwise, remove those elements from the average:
      average_wave[b.:length(average_wave)] <- NA
    }
  }
  
  # Find the median of the first x-values - make that the start of average:
  first_element <- c()
  for(i in 2:ncol(p)){
    stpqw <- p[, i][1:200]
    stq <- stpqw[is.na(stpqw)]
    first_element[i-1] <- length(stq) + 1
  }
  start <- median(first_element)-2
  average_wave[1:start] <- NA
  
  
  return(as.vector(average_wave))
}



find_u_v <- function(dat, wx, wy, d1, d1p, spline, spline_o, plot = FALSE){
  # Find half the height of w (on derivative y-axis)
  w_half_height <- predict(d1p, wx)/2   # should change this to w$w_poly_peaks_deriv1_yval
  # Find u and v:
  half_heights <- c()
  half_heights_yval <- c()
  for(i in 1:length(w_half_height)){    #length(w_half_height)
    deriv1_poly_peak_subset <- CubicInterpSplineAsPiecePoly((round(wx[i])-5):(round(wx[i])+5), d1[(round(wx[i])-5):(round(wx[i])+5)], "natural")   # Making a smaller spline for just the peak
    half_heights_precursor <- solve(deriv1_poly_peak_subset, b = w_half_height[i])   # finding the points on the small spline where the heights are half the peak
    # If only one half height detected, extend window for finding half heights:
    if(length(half_heights_precursor) < 2){
      deriv1_poly_peak_subset <- CubicInterpSplineAsPiecePoly((round(wx[i])-10):(round(wx[i])+10), d1[(round(wx[i])-10):(round(wx[i])+10)], "natural")   # Making a smaller spline for just the peak
      half_heights_precursor <- solve(deriv1_poly_peak_subset, b = w_half_height[i]) 
    }
    # If more than one half height detected:
    if(length(half_heights_precursor) > 2){
      # Find the distance between each detected half height, compare their xvals to the peak, and keep only the two that have the most similar distance to the peak... 
      a <- half_heights_precursor - wx[i]
      b <- c()
      for(j in 1:(length(a)-1)){
        b[j] <- a[j] - a[j+1]
      } 
      b[length(b) + 1] <- a[1] - a[length(a)]
      c <- which(abs(b) == min(abs(b)))
      if(c == length(b)){
        half_heights_precursor <- c(half_heights_precursor[c], half_heights_precursor[1])
      }else{
        half_heights_precursor <- c(half_heights_precursor[c], half_heights_precursor[c+1])
      }
    }
    half_heights[c((2*(i)-1), (2*(i)))] <- half_heights_precursor     # assigning half heights to a vector
    half_heights_yval[c((2*(i)-1), (2*(i)))] <- predict(deriv1_poly_peak_subset, half_heights[c((2*(i)-1), (2*(i)))])  # finding the y_values from those half heights on the smaller spline
  }
  
  # Plot u's and v's on deriv1_poly
  if(plot){
    plot(d1p)
    points(half_heights, half_heights_yval, pch = 19)
  }
  
  # Find u and v 
  u <- half_heights[seq_along(half_heights) %%2 != 0] 
  v <- half_heights[seq_along(half_heights) %%2 == 0]
  
  
  # Find u and v y-values for spline_poly
  u_yval <- c()
  for(i in 1:length(u)){
    u_yval[i] <- predict(spline, u[i])
  }
  
  v_yval <- c()
  for(i in 1:length(v)){
    v_yval[i] <- predict(spline, v[i])
  }
  
  df <- data.frame(u, u_yval, v, v_yval)
  colnames(df) <- c("uX", "uY", "vX", "vY")
  return(df)
}



find_o <- function(wx, inx, iny, d1p, sp){
  o <- c()
  for(i in 1:length(wx)){
    o[i] <- max(which(inx < wx[i]))
  }
  # Adjust for early O points:
  # First find O based on inflection point on first deriv:
  inflexD1 <- solve(d1p, b = 0, deriv = 1)
  o2 <- c()
  for(i in 1:length(wx)){
    o2[i] <- max(which(inflexD1 < wx[i]))
  }
  # Use the O derived from 1st deriv if its y-val is above 0: 
  for(i in 1:length(wx)){
    if((inx[o][i] - inflexD1[o2][i]) < 0){
      inx[o][i] <- inflexD1[o2][i]
      iny[o][i] <- predict(sp, inx[o][i]) 
    }
  }
return(o)
}



find_w <- function(d1p, deriv1, sp, sr){
  
  # sp = spline_poly
  # d1p = deriv1_poly
  # deriv1 = deriv1
  
  inflexion_points_deriv1 <- solve(d1p, b = 0, deriv = 1)                                       # Find inflexion points 1st deriv
  inflexion_points_deriv1_yval <- predict(d1p, inflexion_points_deriv1)
  w_poly_peaks <- c()                                                                           # Create vectors for w values to be stored
  w_poly_peaks_yval <- c()
  window <- list()                                                                              # Define 'window' as the timeframe you will look for peaks in
  p_p_distance <- c()
  
  # Confirm first two peaks:
  a <- 2
  while(length(w_poly_peaks) < 2){
    
    first_window <- data.frame(inflexion_points_deriv1[1]:inflexion_points_deriv1[1+a], deriv1[inflexion_points_deriv1[1]:inflexion_points_deriv1[1+a]])
    inflexion_points_window_yval <- inflexion_points_deriv1_yval[which(inflexion_points_deriv1 > first_window[1, 1] & inflexion_points_deriv1 < first_window[, 1][length(first_window[, 1])])]
    inflexion_points_window <- inflexion_points_deriv1[which(inflexion_points_deriv1 > first_window[1, 1]  & inflexion_points_deriv1 < first_window[, 1][length(first_window[, 1])])]                      # Find inflection within the window from the greater array
    threshold <- quantile(first_window[, 2], probs=c(.95))                                                                                                                                                 # Define a threshold for finding peaks, of 0.95 
    window_poly_peaks <- which(inflexion_points_window_yval > threshold)       
    
    # Very occassionely, the second peak is below the threshold of 0.95, creating a leapfrog error
    # If the number of inflexion points is significant between the first two peaks aka leapfrog error has occurred, check there are no peaks in between that are in fact above threshold
    if(length(window_poly_peaks) == 2 & window_poly_peaks[2] - window_poly_peaks[1] > 40){
      confirm_no_missed_peaks_window <- data.frame(inflexion_points_deriv1[window_poly_peaks[1]+2]:inflexion_points_deriv1[window_poly_peaks[2]-2],   deriv1[inflexion_points_deriv1[window_poly_peaks[1]+2]:inflexion_points_deriv1[window_poly_peaks[2]-2]]  ) 
      inflexion_points_confirm_window_yval <- inflexion_points_deriv1_yval[which(inflexion_points_deriv1 > confirm_no_missed_peaks_window[1, 1] & inflexion_points_deriv1 < confirm_no_missed_peaks_window[, 1][length(confirm_no_missed_peaks_window[, 1])])]
      inflexion_points_confirm_window <- inflexion_points_deriv1[which(inflexion_points_deriv1 > confirm_no_missed_peaks_window[1, 1]  & inflexion_points_deriv1 < confirm_no_missed_peaks_window[, 1][length(confirm_no_missed_peaks_window[, 1])])]                      # Find inflection within the window from the greater array
      threshold <- quantile(confirm_no_missed_peaks_window[, 2], probs=c(.95))       
      peaks_missed <- inflexion_points_confirm_window[which(inflexion_points_confirm_window_yval > threshold)]
      # now match up the first missed peak to inflexion_points in first_window
      peaks_missed <- which(inflexion_points_window == peaks_missed[1])
      # now assign the first missed peak as the second window_poly_peak
      window_poly_peaks[2] <- peaks_missed[1]
    }
    
    # Also very occassionely, the first peak is double-peaked, so make sure peaks are not too close together:
    if(length(window_poly_peaks) == 2 & window_poly_peaks[2] - window_poly_peaks[1] < 4){  
      window_poly_peaks <- window_poly_peaks[-2]   # this line removes the second element, allowing a third to become the second
    }
    
    # Find which inflection points are above threshold
    window_poly_peaks_yval <- inflexion_points_window_yval[window_poly_peaks]
    window_poly_peaks <- inflexion_points_window[window_poly_peaks]     
    #plot(first_window, type = "l")
    #points(window_poly_peaks, window_poly_peaks_yval, pch = 19)
    
    # This is a special case for when the first legitimate peak happens to have 3 very close together inflection points
    if(length(window_poly_peaks) > 2 & (window_poly_peaks[3] - window_poly_peaks[1] < 1)){
      window_poly_peaks <- window_poly_peaks[-c(1, 3)]
    }
    
    m <- mean(inflexion_points_deriv1_yval[1:(1+a)])       # Calculate mean and SD of inflection points within the window
    sd <- sd(inflexion_points_deriv1_yval[1:(1+a)])
    
    if(length(window_poly_peaks) == 2){                    # If two peaks are found, confirm that they are both significantly higher than surrounding inflection points,                                                                                  
      if(window_poly_peaks_yval[1] > (m + (1.5*sd)) & window_poly_peaks_yval[1] > (median(deriv1[1:100])+sd(deriv1[1:100]))){      # and that the second peak is greater than half the height of the first peak (unless the first peak is an artefact)       #MOST RECENT CHANGE 14/10/20: changed the | to & and changed so that only first 100 samples of deriv1 are looked at (rather than all of deriv1)
        w_poly_peaks[1] <-  window_poly_peaks[1]                                                                              
        if(window_poly_peaks_yval[2] > (m + (1.5*sd)) & window_poly_peaks_yval[2] > (window_poly_peaks_yval[1]/2) | window_poly_peaks_yval[1] > (mean(deriv1) + (5*sd(deriv1)))){  
          w_poly_peaks[2] <-  window_poly_peaks[2]                                                                                                                          
        }
      }
    }
    
    if(length(window_poly_peaks) > 3){      # If more than 3 peaks identified, assume the time series begins with an artefact and skip forward
      inflexion_points_deriv1 <- inflexion_points_deriv1[-c(1:(sr*(100/75)))]      
      inflexion_points_deriv1_yval <- inflexion_points_deriv1_yval[-c(1:(sr*(100/75)))]
      a <- -3  # reset a
    }
    
    a <- a + 5      # Widen the window looking for two peaks and continue increasing until two legitimate peaks are found
  }
  
  
  # Find remaining peaks based on p_p distance:
  artefacts <- c()
  # Calculate mean and SD of inflection points generally
  m <- mean(inflexion_points_deriv1_yval)       
  sd <- sd(inflexion_points_deriv1_yval)
  for(i in 3:length(inflexion_points_deriv1)){   # length(inflexion_points_deriv1)
    p_p_distance[i] <- w_poly_peaks[length(w_poly_peaks)] - w_poly_peaks[length(w_poly_peaks)-1]   # p_p_distance defined as distance between two previous peaks
    if(p_p_distance[i] > (mean(p_p_distance[!is.na(p_p_distance)])*2)){    # if p_p_distance is unusually, this can be because there was a large gap between an artefact and the next legitimate peak 
      p_p_distance[i] <- p_p_distance[3]      # correct any abnormal p_p_distances to avoid leapfrogging
    }
    if((w_poly_peaks[length(w_poly_peaks)] + (4*p_p_distance[i]) > length(deriv1))){                            # If the next window goes beyond the length of the data, break the loop
      break
    }                                      
    window_poly_peaks <- c()
    window_extension <- 1.35
    window_beginning <- 0.5
    printed <- NA
    # Adjust mean and standard deviation to remove any artefacts that have been identified
    if(length(artefacts) > 0){
      remove <- c()
      for(j in artefacts[which(artefacts < length(w_poly_peaks))]){
        remove[j] <- which(abs(inflexion_points_deriv1 - w_poly_peaks[j]) == min(abs(inflexion_points_deriv1 - w_poly_peaks[j])))
      }
      remove <- remove[!is.na(remove)]
      new_remove <- c()
      for(j in 1:length(remove)){
        new_remove[(length(new_remove)+1):(length(new_remove)+21)] <- (remove[j] -10): (remove[j] + 10)   # remove inflection points around artefact peaks
      }
      if(sum(new_remove < 1) > 0){  # running whats inside the for loop will result in a null new_remove vector so only run it if there are in fact negative / 0 values in new_remove
        new_remove <- new_remove[-(which(new_remove < 1))] 
      }    
      m <- mean(inflexion_points_deriv1_yval[-new_remove])       
      sd <- sd(inflexion_points_deriv1_yval[-new_remove])  
    }
    
    while(length(window_poly_peaks) < 1){
      
      if(window_extension > 10){
        window_beginning <- 2
        window_extension <- 2.5
      }
      
      window[[i]] <- data.frame((w_poly_peaks[length(w_poly_peaks)] + window_beginning*p_p_distance[i]):(w_poly_peaks[length(w_poly_peaks)] + window_extension*p_p_distance[i]),  deriv1[(w_poly_peaks[length(w_poly_peaks)] + window_beginning*p_p_distance[i]):(w_poly_peaks[length(w_poly_peaks)] + window_extension*p_p_distance[i])])    # window initially defined as (the previous peak + 0.5*p_p_distance) to (previous peak + 1.35*p_p_distance)
      inflexion_points_window_yval <- inflexion_points_deriv1_yval[which(inflexion_points_deriv1 > window[[i]][1, 1] & inflexion_points_deriv1 < window[[i]][length(window[[i]][, 1]), 1])]
      inflexion_points_window <- inflexion_points_deriv1[which(inflexion_points_deriv1 > window[[i]][1, 1] & inflexion_points_deriv1 < window[[i]][length(window[[i]][, 1]), 1])]                          
      threshold <- quantile(window[[i]][, 2], probs=c(0.95))           # If the window is somehow too long for the time series, this line will spring an error                                 
      window_poly_peaks <- which(inflexion_points_window_yval > threshold)                      
      window_poly_peaks_yval <- inflexion_points_window_yval[window_poly_peaks]
      window_poly_peaks <- inflexion_points_window[window_poly_peaks]     
      
      #plot(window[[i]], type = "l")
      #points(window_poly_peaks, window_poly_peaks_yval, pch = 19)
      
      
      if(length(window_poly_peaks) > 2){                  # If three peaks are identified, this can be due to four reasons: 1. a window does not include a genuine peak (so there are multiple secondary ones of similar height) 2. there is an artefact that causes three peaks to be in the same window 3. there are significantly large secondary peaks that they also exceed the threshold for identification 4. the peak of the 1st derivative itself has multiple inflection points e.g in participant 38 this occurs                                                                                                               # If the these criteria are met, go to the next loop
        if(max(window_poly_peaks_yval) > (m+(2*sd)) |  max(window_poly_peaks_yval) > (mean(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]) + 2*sd(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]))  ){        # If peaks are high, include the highest one as a peak. Mark it as an artefact if the other two peaks exceed m+(2*sd) aka are not secondary peaks. 
          w_poly_peaks[length(w_poly_peaks)+1] <- window_poly_peaks[which(window_poly_peaks_yval == max(window_poly_peaks_yval))]  # Mark max peak as a peak
          not_max_peaks <- window_poly_peaks_yval[order(window_poly_peaks_yval)[1:2]]  # Find the two lowest peaks
          if(not_max_peaks[1] > m+(2*sd) & not_max_peaks[2] > m+(2*sd) & (window_poly_peaks[3] - window_poly_peaks[1]) > (p_p_distance[i]/10)){    # If the two lower peaks are also high, and not so close together as to be multiple inflection points on the same peak, label the peak as an artefact. 
            cat('\n','Potential artefact',  ', plot(', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), ', deriv1[', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because the non-max peaks were high') 
            artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2) 
          }
        }else{                                           # If peaks are low, assume they are secondary and extend the window
          window_extension <- window_extension + 0.5
          window_poly_peaks <- c()
        }
      }
      
      
      if(length(window_poly_peaks) == 2){
        if(max(window_poly_peaks_yval) > (m+(2*sd)) | max(window_poly_peaks_yval) > (mean(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]) + 2*sd(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))])) ){   # at least one of the peaks high?    (occassionely there is a series of low amplitude waves that are far below average for the time series, to account for this there is also a test of height relative to the peaks local to it)
          if((window_poly_peaks[2] - window_poly_peaks[1]) < (p_p_distance[i]/3)){   # peaks close togetber?
            w_poly_peaks[length(w_poly_peaks)+1] <- window_poly_peaks[which(window_poly_peaks_yval == max(window_poly_peaks_yval))]  # If they are, take the highest peak and mark it as an artefact
            cat('\n','Potential artefact',  ', plot(', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), ', deriv1[', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because two peaks found too close together') 
          }else{      # If peaks are not too close together..
            if((window_poly_peaks_yval[1] > (m+(2*sd))) | window_poly_peaks_yval[1] > (mean(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]) + 2*sd(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))])) & window_poly_peaks_yval[1] > (window_poly_peaks_yval[2]/2)){ # Check if first peak is legit by 1. comparing to m+(2*sd)  or 2. 
              w_poly_peaks[length(w_poly_peaks)+1] <- window_poly_peaks[1]    # If so mark it as a peak
              if(window_poly_peaks_yval[2] > (m+(2*sd)) |  window_poly_peaks_yval[2] > (mean(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]) + 2*sd(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]))  & window_poly_peaks_yval[2] > (window_poly_peaks_yval[1]/2)){   # and check if the second one is legit...
                w_poly_peaks[length(w_poly_peaks)+1] <- window_poly_peaks[2]  # and if so mark it also as a a peak
              }
            }else{   # If first peak is not legit
              if(window_poly_peaks_yval[2] > (m+(2*sd)) | window_poly_peaks_yval[2] > (mean(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]) + 2*sd(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]))  & window_poly_peaks_yval[2] > (window_poly_peaks_yval[1]/2)){
                w_poly_peaks[length(w_poly_peaks)+1] <- window_poly_peaks[2] 
              }
            }
          }
        }else{     # If peaks not high, extend window
          window_extension <- window_extension + 0.5
          window_poly_peaks <- c()
        }   
      }
      
      
      if(length(window_poly_peaks) == 1){                # If one peak is identified, confirm 1. that the peak is sufficiently high so as not to be a secondary peak 
        if(window_poly_peaks_yval > (m+(2*sd)) | window_poly_peaks_yval > (mean(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))]) + 2*sd(deriv1[window_poly_peaks[1]:(window_poly_peaks[1] + (3*p_p_distance[i]))])) | window_poly_peaks_yval > (predict(d1p, w_poly_peaks[i-1])*0.9)  ){                                                                                                                           
          w_poly_peaks[length(w_poly_peaks)+1] <- window_poly_peaks 
        }else{
          if((i-1) %in% artefacts){
            w_poly_peaks[length(w_poly_peaks)+1] <- window_poly_peaks 
          }else{
            window_extension <- window_extension + 0.5
            window_poly_peaks <- c()
          }
        }
      }
      
      if(window_extension > 2 & is.na(printed)){        # If the window has been extended more than twice to find a peak, it is likely that that peak may be an artefact; label it as such
        cat('\n','Potential artefact',  ', plot(', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), ', deriv1[', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because window needed extending')  # i-1 because w_poly_peaks[i] doesn't exist yet - as such the range has also been adjusted
        printed <- 1 
        artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)
      }
      
      if(min(window[[i]]) < (m-(3*sd)) & is.na(printed)){   # If a window contains a value that drops considerably below the mean, label the peak found in that window as an artefact 
        cat('\n','Potential artefact',  ', plot(', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), ', deriv1[', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because the window contains a very low value') 
        artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)
      }
      
      window_extension <- window_extension + 0.1   # If no peaks are found in a window (or only spurious ones found) extend the window and look again
    }
  }                                                                                                                         
  w_poly_peaks_yval_deriv1 <- predict(d1p, w_poly_peaks)
  #Find m and sd for w_poly_peaks:
  m <- mean(w_poly_peaks_yval_deriv1) 
  sd <- sd(w_poly_peaks_yval_deriv1) 
  # Redefine w_poly_peaks as without artefacts
  if(length(artefacts) > 0){
    w_poly_peaks <- w_poly_peaks[-artefacts]
    w_poly_peaks_yval_deriv1 <- w_poly_peaks_yval_deriv1[-artefacts]
  }
  # Now that all peaks have been found, and artefacts removed, go through once more to look for any unusually tall peaks, and label them as artefacts
  artefacts <- c()
  for(i in 2:(length(w_poly_peaks)-1)){     
    if(w_poly_peaks_yval_deriv1[i] > 1.5*w_poly_peaks_yval_deriv1[i+1] & w_poly_peaks_yval_deriv1[i] > 1.5*w_poly_peaks_yval_deriv1[i-1] | w_poly_peaks_yval_deriv1[i] > (m+(5*sd))){  # consider too tall if greater than 1.5* the wave before and the wave after OR greater than (m+(5*sd)) 
      cat('\n','Potential artefact',  ', plot(', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), ', deriv1[', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because of a very tall peak')
      artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)
    }
    if(i > 5 & i < (length(w_poly_peaks)-5)){
      if(w_poly_peaks_yval_deriv1[i] > (1.5*mean(w_poly_peaks_yval_deriv1[c((i-5):(i+5))]))){
        cat('\n','Potential artefact',  ', plot(', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), ', deriv1[', (w_poly_peaks[i-1]-100), ':', (w_poly_peaks[i-1]+300), '], type = "l") ,', 'wave', i, '+/- 2 removed because of a very tall peak')   # Another check which can only be carried out on waves that have 5 either side of them, is to check if the ith wave is greater than 1.5 x the mean of peaks i-5 to i+5
        artefacts[length(artefacts) + c(1, 2, 3, 4, 5)] <- c(i-2, i-1, i, i+1, i+2)
      }
    }
  }       
  # Find W y_values on original trace
  w_poly_peaks_yval <- predict(sp, w_poly_peaks)  
  # Create dataframe 
  Ws <- data.frame(w_poly_peaks, w_poly_peaks_yval, w_poly_peaks_yval_deriv1)
  if(length(artefacts) > 0){
    Ws <- Ws[-artefacts, ]
  }
  Ws <- Ws[-1, ]
  Ws <-Ws[-nrow(Ws), ]
  
  colnames(Ws) <- c("wX", "wY", "wYD1")
  return(Ws)
}



osnd_of_average <- function(aw, dp, diff, sr){
  
  switch <- 0
  aw <- aw[!is.na(aw)]
  
  # Average wave doesn't need a for loop...
  average_wave_poly <- CubicInterpSplineAsPiecePoly(1:length(aw), aw, "natural")
  
  sfunction <- splinefun(1:length(aw), aw, method = "natural")
  deriv1_wave <- sfunction(1:length(aw), deriv = 1)
  deriv1_wave_poly <- CubicInterpSplineAsPiecePoly(1:length(aw), deriv1_wave, "natural") 
  
  
  # Find inflexion points on deriv1_wave_poly
  inflexion_points_deriv1_wave_poly <- solve(deriv1_wave_poly, b = 0, deriv = 1)
  inflexion_points_deriv1_wave_poly_yval <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)
  
  # Find OSND
  inflexion_points_new <- solve(average_wave_poly, b = 0, deriv = 1)
  inflexion_points_new_yval <- predict(average_wave_poly, inflexion_points_new)
  
  # You might not have to worry about correcting o with this wave since the inflection point in likely to be correct
  
  # Finding notch based on x-axis:
  # Find inflexion point closest to where the notch usually is (aka 75-80)
  notch_range <- which(inflexion_points_deriv1_wave_poly > (3.104572*sr - diff) & inflexion_points_deriv1_wave_poly < dp) # 3.104572 used to be 3.5! #  dp used to be 5*sampling rate!
  # If there is no inflexion point detected within the notch range, this could be because there is a plateu rather than a peak
  # In this case, taking the mean value of the notch range boundaries gives a reasonable approximation
  if(length(notch_range) < 1 | (length(notch_range) == 1 & inflexion_points_deriv1_wave_poly_yval[notch_range][1] < -0.02)){
    new.n <- ((3.104572*sr) + dp)/2
  }else{
    a. <- which(inflexion_points_deriv1_wave_poly_yval[notch_range] == max(inflexion_points_deriv1_wave_poly_yval[notch_range]))
    # In cases where the renal peak is higher on 1st deriv than the notch peak, make sure the notch peak is limited by x-axis
    while(inflexion_points_deriv1_wave_poly[notch_range[a.]] < 115){
      b. <- 2
      a. <- order(inflexion_points_deriv1_wave_poly_yval[notch_range], decreasing = TRUE)[b.]
      b. <- 3
    }
    # Make sure the 1st peak is not the notch:
    if(which(inflexion_points_deriv1_wave_poly_yval == max(inflexion_points_deriv1_wave_poly_yval)) == notch_range[a.]){
      notch_range <- notch_range[which(notch_range > notch_range[a.])]
    }
    new.n <- inflexion_points_deriv1_wave_poly[notch_range[a.]]
  }
  new.ny <- predict(average_wave_poly, new.n)
  #plot(deriv1_wave_poly)
  #points(new.n, predict(deriv1_wave_poly, new.n))
  
  
  # After having found the notch on all waves, you can see if there are inflexion points either side:
  # If inflexion point before is lower and inflexion point after is higher (on y axis), this must mean a second peak aka canonical wave
  # Thus if this criterion is fulfilled you can create N and D separately 
  if(length(inflexion_points_new) > 1 & inflexion_points_new_yval[max(which(inflexion_points_new < new.n))] < new.ny){
    new.n <- inflexion_points_new[max(which(inflexion_points_new < new.n))]
    new.ny <-  predict(average_wave_poly, new.n)
    d. <- inflexion_points_new[min(which(inflexion_points_new > new.n))]
    d.y <- predict(average_wave_poly, d.)
    # If no inflexion point after the notch, take instead the closest value to 0 on the 1st deriv wave after it (the next inflexion point on deriv1)
    if(is.na(d.)){
      d. <- inflexion_points_deriv1_wave_poly[notch_range[a.+1]]
      d.y <- predict(average_wave_poly, d.)
    }
    switch <- 1
  }
  
  plot(average_wave_poly)
  points(inflexion_points_new, inflexion_points_new_yval)
  points(new.n, new.ny, col = "red")
  if(switch == 1){
    points(d., d.y, col = "blue")
  }
  
  
  # Find W: the max inflection point on first deriv:
  w. <- inflexion_points_deriv1_wave_poly[which( inflexion_points_deriv1_wave_poly_yval ==  max(inflexion_points_deriv1_wave_poly_yval) & inflexion_points_deriv1_wave_poly[which(inflexion_points_deriv1_wave_poly_yval ==  max(inflexion_points_deriv1_wave_poly_yval))] < new.n)]
  w._yval <- predict(average_wave_poly, w.)
  points(w., w._yval)
  
  # Find U and V:
  # Find half the height of w (on derivative y-axis)
  hhaw <- max(inflexion_points_deriv1_wave_poly_yval)/2
  # Find u and v for derivative:
  half_heights_wave_new <- solve(deriv1_wave_poly, b = hhaw)
  #half_heights_wave_new_yval <- predict(deriv1_wave_poly, half_heights_wave_new)
  
  # If more than one half height detected:
  if(length(half_heights_wave_new) > 2){
    # Find the distance between each detected half height, compare their xvals to the peak, and keep only the two that have the most similar distance to the peak... 
    # bear in mind one has to be either side of w...
    a <- half_heights_wave_new - w.
    right_of_w <- which(a > 0)
    left_of_w <- which(a < 0)
    half_heights_wave_new <- c(  min(half_heights_wave_new[right_of_w]),  max(half_heights_wave_new[left_of_w]))
    #a <- abs(half_heights_wave_new - w.)
    #b <- c()
    #for(j in 1:(length(a)-1)){
    #  b[j] <- abs(a[j] - a[j+1])
    #} 
    #b[length(b) + 1] <- abs(a[1] - a[length(a)])
    #c <- which(abs(b) == min(abs(b)))
    #if(c == length(b)){
    #  half_heights_wave_new <- c(half_heights_wave_new[c], half_heights_wave_new[1])
    #}else{
    #  half_heights_wave_new <- c(half_heights_wave_new[c], half_heights_wave_new[c+1])
  }
  
  u <- half_heights_wave_new[1]
  v <- half_heights_wave_new[2] 
  # Find u and v y-values for original wave:
  u_v_yval_wave <- predict(average_wave_poly, half_heights_wave_new)
  u_yval <- u_v_yval_wave[1]
  v_yval <- u_v_yval_wave[2]
  points(u, u_yval)
  points(v, v_yval)
  
  # Find 0: the inflection point before W:
  if(length(which(inflexion_points_new < w.)) < 1){
    # If there are no inflexion points on deriv1 before w., or if the max of inflexion points on first deriv is above 0 on the original y-axis, use the first value as O:
    if((length(which(inflexion_points_deriv1_wave_poly < w.)) < 1) | (predict(average_wave_poly, inflexion_points_deriv1_wave_poly[max(which(inflexion_points_deriv1_wave_poly < w.))]) > 0)){
      o. <- 1
    }else{        # Otherwise, make O the maximum deriv1 inflection point before w.
      o. <- inflexion_points_deriv1_wave_poly[max(which(inflexion_points_deriv1_wave_poly < w.))]
    }
  }else{
    o. <- inflexion_points_new[max(which(inflexion_points_new < w.))]
  }
  o._yval <- predict(average_wave_poly, o.)
  points(o., o._yval, pch = 19)
  
  
  
  # Find S: the inflection point after W:   # this bit isn't working for Paul-forms - need to add in a new-s bit. 
  #s. <- inflexion_points_new[min(which(inflexion_points_new > w.))]
  #s._yval <- predict(average_wave_poly, s.)
  #points(s., s._yval, pch = 19)
  
  ## Find S:
  # Find new S:
  s.new <- w. + 2*(abs(v - w.))
  s.new_y <- predict(average_wave_poly, s.new)
  #points(s.new, s.new_y)
  # Define old s:
  s.old_y <- max(inflexion_points_new_yval)
  s.old <- inflexion_points_new[which(inflexion_points_new_yval == max(inflexion_points_new_yval))]
  # Decide which S to use...
  if((s.old - w.) < (s.new - w.)){
    s. <- s.old
    s.y <- s.old_y
  }else{
    s. <- s.new
    s.y <- s.new_y
  }
  points(s., s.y, pch = 19)
  
  if(switch == 0){
    x <- c(o., s., new.n, new.n)
  }else{
    x <- c(o., s., new.n, d.)
  }
  
  if(switch == 0){
    y <- c(o._yval, s.y, new.ny, new.ny)
  }else{
    y <- c(o._yval, s.y, new.ny, d.y)
  }
  
  osnd <- data.frame(x, y)
  return(osnd)
}



preproc <- function(data){
  dat<-data[!(data$PPG.PulseOx1=='NaN'),]
  
  #Downsample
  # The BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second, 
  # so we typically have repeated values in a pattern of 3-3-4 repeating. DownSample tries to retrieve the 
  # unique values, with an effort to be robust against variation in the repeat pattern and also against 
  # genuine repeated values.
  
  list<-rle(dat$PPG.PulseOx1)
  ID <- rep(1:length(list$values), times = list$lengths)
  data_downsampled <- c()
  
  nSrc <- nrow(dat)
  iDst <- 1
  iSrc <- 1
  iVal <- 1
  print("Deduplicating data...")
  if(TRUE){ # Mode switch in case new method is worse or somehow broken
    startClock <- Sys.time()
    reportClock <- 10
    
    data_downsampled <- cbind(dat, ID)
    
    while (iSrc <= nSrc){
      if (as.double(difftime(Sys.time(),startClock,unit="secs"))>reportClock){
        min <- as.integer(reportClock/60)
        sec <- reportClock %% 60
        time <- paste("[",min,":",sep="")
        if (sec == 0){ time <- paste(time,"00]",sep="")} else { time <- paste(time,sec,"]",sep="")}
        print(paste(time,iSrc,"of",nSrc))
        reportClock <- reportClock + 10
      }
      
      data_downsampled[iDst,] <- data_downsampled[iSrc,]
      iDst <- iDst + 1
      if (list$lengths[iVal] <= 4){
        iSrc <- iSrc + list$lengths[iVal]
        iVal <- iVal + 1
      } else {
        iSrc <- iSrc + 4
        list$lengths[iVal] <- list$lengths[iVal] - 4
      }
    }
    
    # Trim downsampled data to size
    data_downsampled <- data_downsampled[1:(iDst-1),]
    
  }else{ # Old method, using rbind which gets slow
    
    data2 <- cbind(dat, ID)
    data_downsampled <-c()
    
    for (i in 1:max(ID)){
      sub.data <- dplyr::filter(data2, ID == i)
      if(nrow(sub.data) <= 4){
        data_downsampled <- rbind(data_downsampled, sub.data[1,])
      }else if(nrow(sub.data) > 4 ){data_downsampled <- rbind(data_downsampled, sub.data[1,], sub.data[5,])}
    }
    
  } # End of method selector
  
  print("Removing DC blocker...")
  #Undetrend
  # Analysis of device output indicates that the PPG signal is detrended by application of the following
  # formula: OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1]), where the constant 0.96875 is
  # an approximation fitted to the data.
  # Individual pulse events are more comprehensible if the detrending is not used, so this function 
  #removes it by inverting the above function. 
  undetrended <-replicate(length(data_downsampled$PPG.PulseOx1)-1, 0) 
  undetrended<-c(data_downsampled$PPG.PulseOx1[1],undetrended) #add first detrended value to vector
  for (i in 2:length(data_downsampled$PPG.PulseOx1)){
    undetrended[i]<-((data_downsampled$PPG.PulseOx1[i]-80) - ((data_downsampled$PPG.PulseOx1[i-1]-80) * 0.96875) + (undetrended[i-1]))
  }
  print("Done")
  return(undetrended)
}



sep_beats <- function(odiff, bc, dat, samp, wuv, wvlen){
  # Redefine baseline corrected data:
  sourcedata <- baseCor[1:length(undetrended)]
  
  # Define a dataframe to contain individual waves (first column is the x-axis (in seconds) - currently set for bioradio data):
  pulse <- data.frame(seq((-141/(samplingRate*10)), ((waveLen*10 -9)-142)/(samplingRate*10), by = 1/(samplingRate*10)))   
  
  
  afterO <- list()
  beforeO <- list()
  extra_long_wave <- c()
  for(i in 1:(length(wuv$wX))){  
    
    # Make a polynomial spline of rounded u - 15 : rounded u + waveLen - 10:   # Now row 141 in xxxx = 0, therefore u = 0 
    splPolySub <- CubicInterpSplineAsPiecePoly((round(wuv$uX[i])-15):(round(wuv$uX[i]) + (waveLen-10)), sourcedata[(round(wuv$uX[i])-15):(round(wuv$uX[i]) + (waveLen-10))], "natural")
    
    # Turn into discrete form
    splSub <- predict(splPolySub, c(seq((wuv$uX[i]-14), (wuv$uX[i]+(waveLen-15)), 0.1)))  
    
    # Make into dataframe:
    splSub <-  as.data.frame(splSub)
    splSub <- cbind(splSub, c(seq((wuv$uX[i]-14), (wuv$uX[i]+(waveLen-15)), 0.1)))
    colnames(splSub) <- c('y', 'x') 
    # Scale so that v-u = 1
    splSub$y <- splSub$y/(wuv$diffVU[i])     
    # Adjust such that u = 0, v = 1 on y-axis
    yDiff <- splSub$y[141]  #??????          
    splSub$y <- splSub$y - yDiff
    
    # Find the x-value for each wave that corresponds to when it = 0.5 in height (this requires making a spline):
    splPolySub2 <- CubicInterpSplineAsPiecePoly(splSub$x, splSub$y, "natural")
    halfCross <- solve(splPolySub2, b = 0.5, deriv = 0)
    halfCross <- halfCross[which(abs(halfCross - wuv$wX[i]) == min(abs(halfCross - wuv$wX[i])))]    
    
    # Convert to discrete form again: (need to redefine splSub)
    splSub2 <- predict(splPolySub, c(seq((halfCross-14), (halfCross+(waveLen-15)), 0.1)))  
    splSub2 <-  as.data.frame(splSub2)
    splSub2 <- cbind(splSub2, c(seq((halfCross-14), (halfCross+(waveLen-15)), 0.1)))  
    colnames(splSub2) <- c('y', 'x') 
    
    # Scale again
    splSub2$y <- splSub2$y/(wuv$diffVU[i]) 
    # Adjust y-axis such that u = 0, v = 1
    yDiff <- wuv$uY[i] / wuv$diffVU[i]
    splSub2$y <- splSub2$y - yDiff
    
    # Find next_o
    afterO[[i]] <- which(splSub2$x > inflexX[o][min(which(inflexX[o] > wuv$wX[i]))])
    
    # Occassionely can get some unusually long waves that have been baseline corrected i.e two systolic peaks merged, these will return integer(o) for the above line - the below checks if o-o difference is unusally large for a wave
    if( (inflexX[o][min(which(inflexX[o] > wuv$wX[i]))]) -  (inflexX[o][max(which(inflexX[o] < wuv$wX[i]))]) > (median(ibi)*1.3)  ){
      extra_long_wave[length(extra_long_wave) + 1] <- i
    }
    
    # Find values before the o of the wave itself 
    beforeO[[i]] <- which(splSub2$x < inflexX[o][max(which(inflexX[o] < wuv$wX[i]))])
    
    # Correct such that x column and wave column are correctly aligned
    splSub3 <- c()
    for(i in 1:nrow(splSub2)){
      splSub3[i+1] <- splSub2$y[i]
    }
    
    # Following lines probably inefficient way of getting everything aligned
    
    # If splSub3 and nrow(pulse) are the same length, you need only adjust afterO
    if(length(splSub3) == nrow(pulse)){
      if(length(afterO[[i]]) > 0){
        diff2 <- length(splSub3) - max(afterO[[i]])
        for(j in 1:diff2){
          afterO[[i]] <- c(afterO[[i]], (max(afterO[[i]]) + 1))
        }
      }
    }
    
    # Or
    # Adjust such that splSub3 is the same length as pulse
    if(length(splSub3) > nrow(pulse)){
      diff <- length(splSub3) - nrow(pulse)
      len <- length(splSub3)
      splSub3 <- splSub3[-((len - (diff-1)):len)]
      if(diff > 1){     # must correct the afterO values so that they also do not contain values beyond the length of splSub3 (include case where length of afterO[[i]] is one so the code works...)
        if(length(afterO[[i]]) > 1 ){
          afterO[[i]] <- afterO[[i]][-(which(afterO[[i]] > length(splSub3)))]  #afterO[[i]][1:(which(afterO[[i]] == (len - (diff-1))) - 1) 
        }else{
          afterO[[i]] <- afterO[[i]][-(which(afterO[[i]] > length(splSub3)))]
        }
      }
    }
    
    if(length(splSub3) < nrow(pulse)){
      diff <- nrow(pulse) - length(splSub3)
      splSub3 <- c(splSub3, rep(NA, diff))
      if(length(afterO[[i]]) > 0){
        diff2 <- length(splSub3) - max(afterO[[i]])
        for(j in 1:diff2){
          afterO[[i]] <- c(afterO[[i]], (max(afterO[[i]]) + 1))
        }
      }
    }
    
    
    # Add column to dataframe
    pulse <- cbind(pulse, splSub3)
  }
  
  
  for(i in 1:(ncol(pulse) -1)){ 
    colnames(pulse)[i+1] <- paste("wave", i, sep = "_")       
  }
  colnames(pulse)[1] <- "x"
  
  # Remove any values after O for each wave:
  for(i in 2:(ncol(pulse))){
    pulse[, i][afterO[[(i-1)]][-1]] <- NA  
  }
  
  # Remove values before O before each wave:
  for(i in 2:(ncol(pulse))){
    pulse[, i][beforeO[[(i-1)]][-1]] <- NA  
  }
  
  # Remove any extra long waves (i.e where a distance of 2 Os has been counted as one wave):
  if(length(extra_long_wave) > 0){
    pulse <- pulse[, -(extra_long_wave + 1)]
  }
  
  # Remove extra tall waves:
  tall_waves <- c()
  for(i in 2:ncol(pulse)){
    if(max(abs(pulse[, i][!is.na(pulse[, i])])) > 1.5){
      tall_waves[i] <- i
    }
  }
  tall_waves <- tall_waves[!is.na(tall_waves)]
  if(length(tall_waves) > 0){
    pulse <- pulse[, -c(tall_waves)]
  }
  average_wave <- find_average(p = pulse, ao = afterO)
  
  dat <- list(average_wave, pulse)
  return(dat)
}



spectrum <- function(baseline_corrected){
  #filtered <- filter.fft(spline_poly,fc=0.0925,BW=0.0525,n=50)
  #plot.fft(filtered)
  #spectral<-spec.fft(spline_poly)
  #plot(data_downsampled$PPG.PulseOx1)
  #plot(baseline_corrected)
  
  #powerspectrum<-spectrum(data_downsampled$PPG.PulseOx1)
  #powerspectrum<-spectrum(baseline_corrected)
  
  #LF<-ffilter(data_downsampled$PPG.PulseOx1, f=75, from = 0.04, to = 0.145, bandpass = TRUE) #low bandpass, sampling frequency of 75 Hz (75 times per second)
  #powerspectrum<-spectrum(LF)
  #HF<-ffilter(data_downsampled$PPG.PulseOx1, f=75, from = 0.145, to = 0.45, bandpass = TRUE) #high bandpass
  #powerspectrum<-spectrum(HF)
  
  #spectralratio<-sum(LF)/sum(HF) #ratio of LF/HF
  
  #Y1 <- fft(baseline_corrected[1:1000])
  
  #plot(abs(Y1), type="h")
  
  #bands<-c(0.04, 0.145, 0.45)
  #frebands<-Freq(baseline_corrected, breaks = bands)
  #Y1 <- fft(frebands[1])
  
  library('signal')
  low<-0.04
  high<-0.145
  bf <- butter(2, c(low, high), type = "pass")
  signal.filtered <- filtfilt(bf, baseline_corrected)
  fourier <- fft(signal.filtered)
  plot(abs(fourier), type="h")
  
  low<-0.145
  high<-0.45
  bf <- butter(2, c(low, high), type = "pass")
  signal.filtered <- filtfilt(bf, baseline_corrected)
  fourier2 <- fft(signal.filtered)
  plot(abs(fourier2), type="h")
  
  spectralratio<-sum(fourier)/sum(fourier2) #ratio of LF/HF
  spectralratio
  return(spectralratio)
  
}
