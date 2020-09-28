# Sorting out the ISO3 data so that the dosage levels are ordered and plotted:

library(pracma)
library(SplinesUtils) 
library(splines2)
library(reshape)
library(ggplot)
library(dplyr)
library(tidyverse)

source("/Users/simonwilliamson/Desktop/Scripts/Find_W_Revised.R")
source("/Users/simonwilliamson/Desktop/Scripts/Rscript_as_function.R")
source("/Users/simonwilliamson/Desktop/Scripts/preproc.R")
source("/Users/simonwilliamson/Desktop/Scripts/find_osnd.R")
source("/Users/simonwilliamson/Desktop/Scripts/spectrum.R")
source("/Users/simonwilliamson/Desktop/Scripts/find_w_u_v.R")
source("/Users/simonwilliamson/Desktop/Scripts/fit_sines.R")
source("/Users/simonwilliamson/Desktop/Scripts/refit_peaks.R")
source("/Users/simonwilliamson/Desktop/Scripts/Find_W_Revised.R")
source("/Users/simonwilliamson/Desktop/Scripts/find_average.R")
source("/Users/simonwilliamson/Desktop/Scripts/osnd_of_average.R")

# Read in the run_order txt file
dose <-read.table("/Users/simonwilliamson/Desktop/ISO_3_run_dose_order.txt", header=TRUE)

# Run ecg attempt script to get the nam and test lists

setwd("/Users/simonwilliamson/Desktop/Khalsa_Fletcher_collaboration copy/Physio Dial/ISO_3.0/Iso_3_PhysioData/")

#List directories in current wd
direc <- list.dirs( full.names = TRUE)
#Find those with ECG files in them
direc <- direc[grep("physiological", direc)]  
dose[,3] <- as.character(dose[,3])
test_nam <- c()
for(i in 1:length(direc)){test_nam[i] <- substr(direc[i], start= 3, stop = 7)}
#Create a huge list with participant codes and each ppg reading 
test <- list()
for(i in 1:length(direc)){
  p <- list.files(direc[i], pattern="ECG_", full.names = TRUE)
  p <- subset(p, !(grepl("REST", p)), drop = TRUE)
  if (length(p) == 6){
    varnam = substr(direc[i], start= 3, stop = 7)
    order <- unlist(strsplit(dose[dose[,1]==varnam,3], ","))
    s <- list()
    for(j in 1:6){
      label <- order[j]
      num <- sum(order[1:j]==label)
      
      t <- read.table(p[j])
      t <- unlist(t, use.names = FALSE)
      s[[sprintf("%s_%i",label,num)]] <- t
    }
    test[[varnam]] <- s
  }
}
test <- Filter(function(x) length(x) == 6, test) ##We're only interested in 6 measures - the iso measures
nam = names(test)

ppg_spl <- list()
dv1 <- list()
deriv1_dis <- list()
#Create splines of the ppg trace, the first derivative and discrete vals of the first deriv
for(i in 1:length(test)){
  len = length(test[[i]])
  spl <- list()
  d1 <- list()
  d1d <- list()
  for(j in 1:len){
    sfunction <- splinefun(1:length(test[[i]][[j]]), test[[i]][[j]], method = "natural")
    spl[[j]] <- CubicInterpSplineAsPiecePoly(1:length(test[[i]][[j]]), test[[i]][[j]], "natural")
    d1_temp <- sfunction(1:length(test[[i]][[j]]), deriv = 1)
    d1d[[j]] <- d1_temp
    d1[[j]] <- CubicInterpSplineAsPiecePoly(1:length(test[[i]][[j]]), d1_temp, 'natural')
  }
  varnam = substr(direc[i], start= 3, stop = 7)
  ppg_spl[[varnam]] <- spl
  dv1[[varnam]] <- d1
  deriv1_dis[[varnam]] <- d1d
}
excl <- c()
excl2 <- c()
#Go through and look for the traces where there are lots of occurences of floor (-2048)
# and ceiling (2046) values - if there are more than 100, add participant to the excl/excl2 arrays
for(i in 1:length(test)){
  len = length(test[[i]])
  for(j in 1:len){
    n <- floor(length(test[[i]][[j]])/1000)
    nGood <- 0
    for(k in 1:n)
    {
      temp <- test[[i]][[j]][(1000*k-999):(1000*k)]
      maxvals <- sum(temp == 2046 | temp == -2048)
      if (maxvals == 0){
        nGood <- nGood + 1
      }
    }
    if (nGood < 2){
      excl <- append(excl, (paste(nam[i], j)))
      excl2 <- append(excl2, (paste(nam[i])))
    }
  }
}
#These participants have all been visually identified as needing to be removed
remall <- c("AL826", "AO139", "AP493", "AU602", "AZ985", "AZ883")
#Remove them from all of the different lists
for(i in 1:length(remall)){
  test[[remall[i]]] <- NULL
  ppg_spl[[remall[i]]] <- NULL
  dv1[[remall[i]]] <- NULL
  deriv1_dis[[remall[i]]] <- NULL
}

#Change the nam array to reflect exclusions 
nam <- names(test)


# Create a data_frame based on nam and test
run_orders <- data.frame(nam)

nums <- list()
for(i in 1:length(test)){
  nums[[i]] <- rownames(summary(test[[i]]))
}

dose_unordered <- c()
for(i in 1:length(nums)){
  matches <- regmatches(nums[[i]], gregexpr("[[:digit:]]+", nums[[i]]))
  doses <-  as.numeric(unlist(matches))[c(1, 3, 5, 7, 9, 11)]
  doses[which(doses == 5)] <- 0.5
  dose_unordered[i] <- toString(doses)
}

run_orders <- cbind(run_orders, dose_unordered)

# Now run that dataframe through these decoding lines to get the correct order for each participant:
dose_orders <- list()
for(i in 1:nrow(run_orders)){
  x <- run_orders[i, 2]
  y <- 1:6 
  matches <- regmatches(x, gregexpr("[[:digit:]]+", x))
  x <- as.numeric(unlist(matches))
  x <- x[-(which(x == 5)-1)]
  x[which(x == 5)] <- 0.5
  x <- order(x)
  y <- y[x]
  dose_orders[[i]] <- y
}

# Now you have dose_orders, a list of orders which corresponds to nam. 

# Run through the main R script to get averaged waves, then arrange and plot the averaged waves in order:
# Black = 0mg, blue = 0.5mg, red = 2mg

for(k in 99:99){   # 1:length(test)       
  len = length(test[[k]])   
  # Make a list of the average waves 
  dose_ordered_average_waves <- list()
  for(l in 1:len){         
    data <- data.frame(test[[c(k, l)]])
    data <- cbind(data, (1:length(data)))
    undetrended_data <- data
    colnames(undetrended_data)[1] <- "undetrended"
    sampling_rate <- 40
    sfunction <- splinefun(1:length(undetrended_data$undetrended), undetrended_data$undetrended[1:length(undetrended_data$undetrended)], method = "natural")
    deriv1 <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 1)
    spline_1 <-  sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 0)
    
    # Create polynomial splines: 
    spline_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), undetrended_data$undetrended[1:length(undetrended_data$undetrended)], "natural")
    deriv1_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), deriv1, "natural") 
    
    # Finding inflexion points on spline_poly (points on deriv1 where x = 0):
    inflexion_points <- solve(spline_poly, b = 0, deriv = 1)
    inflexion_points_yval <- predict(spline_poly, inflexion_points)
    
    # Find W:
    w <- find_w(d1p = deriv1_poly, deriv1 = deriv1, sp = spline_poly)
    #plot(deriv1_poly)
    #points(w$w_poly_peaks, w$w_poly_peaks_yval_deriv1)
    
    # Need to address the warning message here... BG754 of ear data is a good example
    
    # Find U and V:
    u_v <- find_u_v(dat = undetrended_data$undetrended, wx = w$w_poly_peaks, wy = w$w_poly_peaks_yval, d1 = deriv1, d1p = deriv1_poly, spline = spline_poly, spline_o = spline_1, plot=FALSE)
    #plot(6400:6600, spline_1[6400:6600], type = "l")
    #points(u_v$u, u_v$u_yval)
    #points(u_v$v, u_v$v_yval)
    #points(w$w_poly_peaks, w$w_poly_peaks_yval)
    
    # Find O in order to find the baseline:
    o <- c()
    for(i in 1:length(w$w_poly_peaks)){
      o[i] <- max(which(inflexion_points < w$w_poly_peaks[i]))
    }
    #plot(spline_poly)
    #points(inflexion_points[o], inflexion_points_yval[o], pch = 19)
    # Find O based on inflexion point on first deriv:
    #plot(1:500, deriv1[1:500], type = "l")
    inflexion_points_d1 <- solve(deriv1_poly, b = 0, deriv = 1)
    #points(inflexion_points_d1, predict(deriv1_poly, inflexion_points_d1))
    o2 <- c()
    for(i in 1:length(w$w_poly_peaks)){
      o2[i] <- max(which(inflexion_points_d1 < w$w_poly_peaks[i]))
    }
    #points(inflexion_points_d1[o2], predict(deriv1_poly, inflexion_points_d1[o2]), col = "red")
    #plot(spline_1[1:500], type = "l")
    #points(inflexion_points_d1[o2], predict(spline_poly, inflexion_points_d1[o2]), col = "red")
    #  only use the O based on 1st deriv if its y-val is above 0 
    for(i in 1:length(w$w_poly_peaks)){
      if((inflexion_points[o][i] - inflexion_points_d1[o2][i]) < 0){
        inflexion_points[o][i] <- inflexion_points_d1[o2][i]
        inflexion_points_yval[o][i] <- predict(spline_poly, inflexion_points[o][i]) 
      }
    }
    #plot(spline_1[1:500], type = "l")
    #points(inflexion_points[o], inflexion_points_yval[o])
    
    
    # Adjust for early O points: #Simon will figure this bit out 
    #for(i in 1:length(w$w_poly_peaks)){
    #  o_decider <- w$w_poly_peaks[i] - 2*(w$w_poly_peaks[i] - u_v$u[i])
    #  if(abs(o_decider - inflexion_points[o[i]]) > 1.5){
    #    inflexion_points[o[i]] <- o_decider
    #    inflexion_points_yval[o[i]] <- predict(spline_poly, o_decider)
    #  }
    #}
    
    # Find where W lies from U to V on the x-axis, for each wave:
    percentage_distance_to_v <- c()
    for(i in 1:length(w$w_poly_peaks)){
      percentage_distance_to_v[i] <- (w$w_poly_peaks[i] - u_v$u[i]) / (u_v$v[i] - u_v$u[i])*100
    }
    #plot(w$w_poly_peaks, percentage_distance_to_v, type = "l")
    sdpdtv <- sd(percentage_distance_to_v)
    #plot(percentage_distance_to_v)
    #points(which(percentage_distance_to_v > (median(percentage_distance_to_v) + sdpdtv)), percentage_distance_to_v[which(percentage_distance_to_v > (median(percentage_distance_to_v) + sdpdtv))], pch = 19, col = "red")
    #points(which(percentage_distance_to_v < (median(percentage_distance_to_v) - sdpdtv)), percentage_distance_to_v[which(percentage_distance_to_v < (median(percentage_distance_to_v) - sdpdtv))], pch = 19, col = "red")
    
    # If u and v have been incorrectly identified because of artefactual inteference in the systolic upstroke, 
    # the position of w relative to them is likely to be different to normal, so such waves can be identified and removed. 
    
    # Make a vector of abnormal pdtv (allowing any values between 35 and 65): 
    pdtv_waves <- c(which(percentage_distance_to_v > (sdpdtv + median(percentage_distance_to_v)) & percentage_distance_to_v > 65), which(percentage_distance_to_v < (median(percentage_distance_to_v) - sdpdtv) & percentage_distance_to_v < 35))
    
    # Remove pdtv_waves from all vectors that include their info...:
    if(length(pdtv_waves) > 0){
      w <- w[-pdtv_waves, ]
      u_v <- u_v[-pdtv_waves, ]
    }
    
    
    ########################################################################    
    
    #             Step 3 : Correct Baseline and refind W, U, V and O       #
    
    ######################################################################## 
    
    # Correct Baseline
    baseline_corrected <- baseline(plot=FALSE)
    
    # Redefine discrete splines:
    sfunction_bc <- splinefun(1:length(baseline_corrected), baseline_corrected, method = "natural")
    deriv1_bc <- sfunction_bc(seq(1, length(baseline_corrected)), deriv = 1)
    deriv2_bc <- sfunction_bc(seq(1, length(baseline_corrected)), deriv = 2)
    spline_1_bc <- sfunction_bc(seq(1, length(baseline_corrected)), deriv = 0)
    
    # Redefine polynomial splines: 
    spline_poly_bc <- CubicInterpSplineAsPiecePoly(1:length(baseline_corrected), baseline_corrected, "natural")
    deriv1_poly_bc <- CubicInterpSplineAsPiecePoly(1:length(baseline_corrected), deriv1_bc, "natural") 
    
    # Refind W y values:
    w$w_poly_peaks_yval <- predict(spline_poly_bc, w$w_poly_peaks)
    
    # Refind U and V y-values:
    u_v$u_yval <- predict(spline_poly_bc, u_v$u)
    u_v$v_yval <- predict(spline_poly_bc, u_v$v)
    
    
    # Check how far u is from baseline...:
    #plot(u_v$u_yval)
    median(u_v$u_yval) - sd(u_v$u_yval)
    # Remove u values that are implausibly far from baseline
    false_u_vals <- which(u_v$u_yval > (median(u_v$u_yval) + sd(u_v$u_yval)) | u_v$u_yval < -50)
    if(length(false_u_vals) > 1){
      w <- w[-false_u_vals, ]
      u_v <- u_v[-false_u_vals, ]
      # v_minus_u <- v_minus_u[-false_u_vals]
    }
    
    # Find v-u differences:
    v_minus_u <- u_v$v_yval - u_v$u_yval
    
    # undetected artefacts can lead to repeats in half heights (if the deriv peak is too wide), which will lead to 
    # scale factors of 0, check for these:
    if(length(which(v_minus_u == 0)) > 0){
      dup <- which(v_minus_u == 0)
      w <- w[-dup, ]
      u_v <- u_v[-dup, ]
      v_minus_u <- v_minus_u[-dup]
    }
    
    # Undetected artefactual waves can cause unusally small scale factors - these should be removed:
    #plot(v_minus_u)
    #points(which(v_minus_u < (median(v_minus_u) - 3*(sd(v_minus_u)))), v_minus_u[which(v_minus_u < (median(v_minus_u) - 3*(sd(v_minus_u))))], pch = 19, col = "red")
    # Make a vector of waves with such scale factors and remove them:
    false_scale_waves <- which(v_minus_u < (median(v_minus_u) - 3*(sd(v_minus_u))))
    if(length(false_scale_waves) > 1){
      w <- w[-false_scale_waves, ]
      u_v <- u_v[-false_scale_waves, ]
      v_minus_u <- v_minus_u[-false_scale_waves]
    }
    
    # Find o points again:     
    o_yval <- predict(spline_poly_bc, inflexion_points[o])
    #plot(spline_poly_bc)
    #points(inflexion_points[o], o_yval, pch = 19)
    
    # Find o-w difference:
    o_w_difference <- c()
    for(i in 1:length(w$w_poly_peaks)){
      o_w_difference[i] <- w$w_poly_peaks[i] - inflexion_points[o[i]]
    }
    
    # Find distance between o_points:
    o_difference <- c()
    for(i in 1:(length(inflexion_points[o])-1)){
      o_difference[i] <- inflexion_points[o[i+1]] - inflexion_points[o[i]]
    }
    
    # Find distance between W points (Inter-beat interval):
    ibi <- c()
    for(i in 1:(length(w$w_poly_peaks))-1){
      ibi[i] <- w$w_poly_peaks[i+1] - w$w_poly_peaks[i]
    }
    # Also remove the last wave:
    w <- w[-nrow(w), ]
    #plot(ibi)
    #points(which(ibi > 1.3*median(ibi)), ibi[which(ibi > 1.3*median(ibi))], pch = 19, col = "red")
    
    # Make a vector of abnormal IBIs (the waves at the end of a sequence / before an artefact): 
    end_waves <- which(ibi > 1.3*median(ibi))    
    
    # Remove end_waves from all vectors that include their info...:
    if(length(end_waves) > 1){
      w <- w[-end_waves, ]
      u_v <- u_v[-nrow(u_v), ]
      u_v <- u_v[-end_waves, ]
      v_minus_u <- v_minus_u[-length(v_minus_u)]
      v_minus_u <- v_minus_u[-end_waves]
    }
    
    ########################################################################    
    
    #       Step 4 : Find individual waves and the average wave            #
    
    ######################################################################## 
    
    # Find the average length of a wave (and 15 since we are starting the wave from before O):
    source_data_column_length <- round(median(o_difference)+15) 
    
    # Redefine baseline corrected data:
    sourcedata <- baseline_corrected[1:length(undetrended_data$undetrended)]
    
    # Define a dataframe to contain individual waves (first column is the x-axis (in seconds)):
    pulse <- data.frame(seq((-141/(sampling_rate*10)), ((source_data_column_length*10 -9)-142)/(sampling_rate*10), by = 1/(sampling_rate*10)))   
    
    # Add waves to the dataframe:
    after_o <- list()
    before_o <- list()
    extra_long_wave <- c()
    for(i in 1:(length(w$w_poly_peaks))  ){   #(length(w$w_poly_peaks))  
      
      # Make a polynomial spline of rounded u - 15 : rounded u + source_data_column_length - 10:   # Now row 141 in xxxx = 0, therefore u = 0 
      spline_poly_wave_subset <- CubicInterpSplineAsPiecePoly((round(u_v$u[i])-15):(round(u_v$u[i]) + (source_data_column_length-10)), sourcedata[(round(u_v$u[i])-15):(round(u_v$u[i]) + (source_data_column_length-10))], "natural")
      
      # Turn into discrete form
      xxxx <- predict(spline_poly_wave_subset, c(seq((u_v$u[i]-14), (u_v$u[i]+(source_data_column_length-15)), 0.1)))  
      
      # Make into dataframe
      xxxx <-  as.data.frame(xxxx)
      xxxx <- cbind(xxxx, c(seq((u_v$u[i]-14), (u_v$u[i]+(source_data_column_length-15)), 0.1)))
      colnames(xxxx) <- c('y', 'x') 
      # Scale so that v-u = 1
      xxxx$y <- xxxx$y/(v_minus_u[i])     
      # Adjust such that u = 0, v = 1 on y-axis
      y_axis_difference <- xxxx$y[141]           
      xxxx$y <- xxxx$y - y_axis_difference
      
      # At this point, the 141st element is u. This means you know the exact x value for u. 
      # You should also have the exact x value for all the o's
      # The second o, i.e the first next_o, is 279.2701
      # You can find out which element of xxxx is the first to exceed this value, and mark it
      # once the final data frame is formed, you can then make NAs all elements after the marked one
      
      # Find the x-value for each wave that corresponds to when it = 0.5 in height (this requires making a spline)
      spline_poly_wave_subset_2 <- CubicInterpSplineAsPiecePoly(xxxx$x, xxxx$y, "natural")
      half_crossing <- solve(spline_poly_wave_subset_2, b = 0.5, deriv = 0)
      half_crossing <- half_crossing[which(abs(half_crossing - w$w_poly_peaks[i]) == min(abs(half_crossing - w$w_poly_peaks[i])))]    # this should be the one that is closest to w (sections from the previous and following waves can also cross 0.5)
      
      
      # Convert to discrete form again: (need to redefine xxxx):
      xxxx2 <- predict(spline_poly_wave_subset, c(seq((half_crossing-14), round((half_crossing+(source_data_column_length-15))), 0.1)))    
      xxxx2 <-  as.data.frame(xxxx2)
      xxxx2 <- cbind(xxxx2, c(seq((half_crossing-14), round((half_crossing+(source_data_column_length-15))), 0.1))) 
      colnames(xxxx2) <- c('y', 'x') 
      
      # Scale again
      xxxx2$y <- xxxx2$y/(v_minus_u[i]) 
      # Adjust y-axis such that u = 0, v = 1
      y_axis_difference <- u_v$u_yval[i] / v_minus_u[i]
      xxxx2$y <- xxxx2$y - y_axis_difference
      
      # Find next_o
      # On some waves there is next identifiable o is after the wave, so it doesn't get shortened e.g 84
      after_o[[i]] <- which(xxxx2$x > inflexion_points[o][min(which(inflexion_points[o] > w$w_poly_peaks[i]))])
      
      # Occassionely can get some unusually long waves that have been baseline corrected i.e two systolic peaks merged, these will return integer(o) for the above line - the below checks if o-o difference is unusally large for a wave
      if( (inflexion_points[o][min(which(inflexion_points[o] > w$w_poly_peaks[i]))]) -  (inflexion_points[o][max(which(inflexion_points[o] < w$w_poly_peaks[i]))]) > (median(ibi)*1.3)  ){
        extra_long_wave[length(extra_long_wave) + 1] <- i
      }
      
      # Find values before the o of the wave itself 
      before_o[[i]] <- which(xxxx2$x < inflexion_points[o][max(which(inflexion_points[o] < w$w_poly_peaks[i]))])
      
      # Correct such that x column and wave column are correctly aligned
      xxxx3 <- c()
      for(j in 1:nrow(xxxx2)){
        xxxx3[j+1] <- xxxx2$y[j]
      }
      #xxxx3 <- xxxx3[-length(xxxx3)]
      
      # If xxxx3 and nrow(pulse) are the same length, you need only adjust after_o
      if(length(xxxx3) == nrow(pulse)){
        if(length(after_o[[i]]) > 0){
          diff2 <- length(xxxx3) - max(after_o[[i]])
          for(j in 1:diff2){
            after_o[[i]] <- c(after_o[[i]], (max(after_o[[i]]) + 1))
          }
        }
      }
      
      # Or
      # Adjust such that xxxx3 is the same length as pulse
      if(length(xxxx3) > nrow(pulse)){
        diff <- length(xxxx3) - nrow(pulse)
        len1 <- length(xxxx3)
        xxxx3 <- xxxx3[-((len1 - (diff-1)):len1)]
        if(diff > 1){     # must correct the after_o values so that they also do not contain values beyond the length of xxxx3 (include case where length of after_o[[i]] is one so the code works...)
          if(length(after_o[[i]]) > 1 ){
            after_o[[i]] <- after_o[[i]][-(which(after_o[[i]] > length(xxxx3)))]  #after_o[[i]][1:(which(after_o[[i]] == (len - (diff-1))) - 1) 
          }else{
            after_o[[i]] <- after_o[[i]][-(which(after_o[[i]] > length(xxxx3)))]
          }
        }
      }
      
      if(length(xxxx3) < nrow(pulse)){
        diff <- nrow(pulse) - length(xxxx3)
        xxxx3 <- c(xxxx3, rep(NA, diff))
        # You do need to adjust after_o as well, since you are increasing length of xxxx2 to fit the nrow of pulse I guess
        if(length(after_o[[i]]) > 0){
          diff2 <- length(xxxx3) - max(after_o[[i]])
          for(j in 1:diff2){
            after_o[[i]] <- c(after_o[[i]], (max(after_o[[i]]) + 1))
          }
        }
      }
      
      # Add column to dataframe
      pulse <- cbind(pulse, xxxx3)
    }
    for(i in 1:(ncol(pulse) -1) ){   #(ncol(pulse) -1)       
      colnames(pulse)[i+1] <- paste("wave", i, sep = "_")       
    }
    colnames(pulse)[1] <- "x"
    
    # Remove any values after O for each wave:
    for(i in 2:(ncol(pulse))){
      pulse[, i][after_o[[(i-1)]][-1]] <- NA  
    }
    # (xxxx3 is just xxxx2 moved down a row - hence the [-1] above)
    # Remove values before O before each wave:
    for(i in 2:(ncol(pulse))){
      pulse[, i][before_o[[(i-1)]][-1]] <- NA  
    }
    
    # Remove any extra long waves (i.e where a distance of 2 Os has been counted as one wave):
    if(length(extra_long_wave) > 0){
      pulse <- pulse[, -(extra_long_wave + 1)]
    }
    
    # Some waves are still abnormally large (one reason is the rise in baseline being recognized as a peak?)
    # Remove tall waves:
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
    
    # Find the average wave:
    average_wave <- find_average(p = pulse, ao = after_o)
    #plot(average_wave)
    dose_ordered_average_waves[[l]] <- average_wave    
    
  }   
  # Now correct the order of the waves:
  dose_ordered_average_waves <- dose_ordered_average_waves[dose_orders[[k]]]
  plot(dose_ordered_average_waves[[1]], type = "l", main = nam[k], ylab = "averaged values", xlab = "")
  lines(dose_ordered_average_waves[[2]])
  lines(dose_ordered_average_waves[[3]], col = "blue")
  lines(dose_ordered_average_waves[[4]], col = "blue")
  lines(dose_ordered_average_waves[[5]], col = "red")
  lines(dose_ordered_average_waves[[6]], col = "red")
}
