# Calculating pulse transit time: 

# Note need packages from main R script loaded

# Import:
data <- read.csv(file.choose())
data <- data[1:10000, ]

# Visualise
par(mfrow = c(3, 1))
plot(data$Elapsed.Time, data$ECG, type = "l")
plot(data$Elapsed.Time, data$PPG.PulseOx1, type = "l")
plot(data$Elapsed.Time, data$RIP, type = "l")

# Should make sure the ECG events can be related to the PPG time series even after the latter has been pre-processed

ecg_and_time <- data[, c(1, 5)]
undetrended_data <- data.frame(preproc(dat=data))

par(mfrow = c(2, 1))
plot(undetrended_data$Elapsed.Time, undetrended_data$PPG.PulseOx1, type = "l")
plot(ecg_and_time$Elapsed.Time, ecg_and_time$ECG, type = "l")

# Find peaks on ECG:

sfunction <- splinefun(1:length(ecg_and_time$ECG), ecg_and_time$ECG[1:length(ecg_and_time$ECG)], method = "natural")
deriv1 <- sfunction(seq(1, length(ecg_and_time$ECG)), deriv = 1)

# The derivative is better for peak detection, but this method will probably need revision
par(mfrow = c(1, 1))
plot(deriv1, type = "l")
x. <- findpeaks(deriv1, minpeakheight = 0.0001)

# Use the peak of the derivatives to help locate the actual peaks
ecg_peaks <- c()
for(i in 1:nrow(x.)){
  y. <- which(ecg_and_time$ECG[x.[, 2][i]:(x.[, 2][i] + 20)] == max(ecg_and_time$ECG[x.[, 2][i]:(x.[, 2][i] + 20)])) 
  ecg_peaks[i] <-  x.[, 2][i] + (y.-1)
}
plot(ecg_and_time$ECG, type = "l")
points(ecg_peaks, ecg_and_time$ECG[ecg_peaks])

R_peak_timings <- ecg_and_time$Elapsed.Time[ecg_peaks]

# Calculate interbeat interval
R_R_interval <- c()
for(i in 1:length(R_peak_timings)){
  R_R_interval[i] <- R_peak_timings[i+1] - R_peak_timings[i]
}

# At this point you can remove artefacts (interbeat intervals that are unusually small):
artifacts <- which(R_R_interval < (mean(R_R_interval[!is.na(R_R_interval)]) - sd(R_R_interval[!is.na(R_R_interval)])))
R_R_interval[artifacts] <- NA 
R_peak_timings[artifacts] <- NA

# Combine R-R_intervals with the ecg_and_time dataframe, lining them up with the peaks
rri <- 1:nrow(ecg_and_time)
rri[1:nrow(ecg_and_time)] <- NA
rri[ecg_peaks] <- R_R_interval

ecg_and_time <- cbind(ecg_and_time, rri)
colnames(ecg_and_time)[3] <- "interbeat interval"

# Plot the interbeat column:
par(mfrow = c(2, 1))
plot(ecg_and_time$Elapsed.Time[(ecg_peaks)], ecg_and_time$`interbeat interval`[(ecg_peaks)], type = "l", xlim = c(ecg_and_time$Elapsed.Time[1], ecg_and_time$Elapsed.Time[nrow(ecg_and_time)]), ylim = c())   
plot(ecg_and_time$Elapsed.Time, ecg_and_time$ECG, type = "l")


# Calculate PTT (currently calculating this as time from R peak on ECG to W on PPG):

# Change the name of the PPG column so it can be fed into main R script:
colnames(undetrended_data)[7] <- "undetrended"
# Create discrete splines:
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


# Converting W points back to elapsed_time values:
# Note: W is in index values with decimal places, whereas elapsed time 
# is discrete. R will normally round w to the nearest elapsed time value in the dataframe.
actual_w_time <- c()
for(i in 1:length(w$w_poly_peaks)){
  dif <- w$w_poly_peaks[i] - floor(w$w_poly_peaks[i])
  n. <- undetrended_data$Elapsed.Time[floor(w$w_poly_peaks[i])]
  n.. <- undetrended_data$Elapsed.Time[floor(w$w_poly_peaks[i])+1]
  actual_w_time[i] <- time_dif*dif + n.
}
par(mfrow = c(1, 1))
plot(undetrended_data$Elapsed.Time, undetrended_data$undetrended, type = "l")
points(actual_w_time, w$w_poly_peaks_yval)


# Plot ECG and PPG waves:
par(mfrow = c(1, 1))
plot(ecg_and_time$Elapsed.Time, (ecg_and_time$ECG - mean(ecg_and_time$ECG)), type = "l")
points(R_peak_timings, (ecg_and_time$ECG[ecg_peaks] -mean(ecg_and_time$ECG)), pch = 19)
lines(undetrended_data$Elapsed.Time, (undetrended_data$undetrended - mean(undetrended_data$undetrended))/15000, col = "red")
points(actual_w_time, (w$w_poly_peaks_yval - mean(undetrended_data$undetrended))/15000, pch = 19, col = "red")


# Removing the waves:
plot(R_peak_timings, rep(0.2, length(R_peak_timings)), pch = 19, ylim = c(0.15, 0.35))
points(actual_w_time, rep(0.2, length(actual_w_time)), pch = 19, col = "red")


# Calculate pulse transit times:
pulse_transit_time <- c()
for(i in 1:length(R_peak_timings)){
  if(is.na(R_peak_timings[i])){
    next
  }else{
    higher_vals <- actual_w_time - R_peak_timings[i]
    higher_vals <- higher_vals[which(higher_vals > 0)]
    if(length(higher_vals) < 1){
      break
    }
    next_w <- higher_vals[which(higher_vals == min(higher_vals))]
    next_w <- which(actual_w_time - R_peak_timings[i] == next_w)
    if(is.na(R_peak_timings[i+1])){
      pulse_transit_time[i] <- actual_w_time[next_w] - R_peak_timings[i] 
      next
    }
    if(R_peak_timings[i+1] < actual_w_time[next_w]){
      next
    }else{
      pulse_transit_time[i] <- actual_w_time[next_w] - R_peak_timings[i] 
    }
  }
}

#Plot:
lines(R_peak_timings, pulse_transit_time[1:length(R_peak_timings)])

# Mean
mn <- mean(pulse_transit_time[!is.na(pulse_transit_time)])

# Relation to respiration:
par(mfrow = c(2, 1))
plot(R_peak_timings, pulse_transit_time[1:length(R_peak_timings)], type = "l")
plot(data$Elapsed.Time, data$RIP, type = "l")

par(mfrow = c(1, 1))
plot(R_peak_timings, (pulse_transit_time[1:length(R_peak_timings)] - mn), type = "l")
lines(data$Elapsed.Time, (data$RIP - mean(data$RIP)), col = "red")

# Relevant literature on this relationship: https://link.springer.com/article/10.1007/s11517-006-0064-y
