setwd("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/PulseAnalysis/Data/Craig")

library(tidyverse)                                 
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils) #SplinesUtils is best downloaded directly from Github
library(spectral)
library(seewave)

source("preproc.R")
source("find_osnd.R")
source("spectrum.R")
source("find_w_u_v.R")
source("fit_sines.R")
source("refit_peaks.R")
source("Find_W_Revised.R")


########################################################################    

#             Step 1 : Reading in Data and Preprocessing               #

######################################################################## 


#### Bioradio (finger) ####
data <- read.csv("~/Desktop/PulseAnalysis/Data/Craig/Source.csv", header = T)
#downsampling data and undetrending 
undetrended_data <- data.frame(preproc(dat=data))


#### ISO-3 (finger) #### 
data <- read.delim("~/Desktop/tests/test.txt", header = T).  
data <- cbind(data, (1:nrow(data)))
undetrended_data <- data
colnames(undetrended_data)[1] <- "undetrended"


#### Intero-battery (ear) ####
data <- read.csv("~/Desktop/Khalsa_Fletcher_collaboration copy/Intero_Battery/AN Float/data-organized/AT053/T0/behavioral_session/AT053-T0-__BH-R1-PHYS.csv")
data <- data[, c(1, 4)]
data <- data[seq(1, nrow(data), 100), ]
undetrended_data <- data
colnames(undetrended_data)[2] <- "undetrended"
# Most ear PPG traces have large periods of nothing followed by large artefact followed by waveforms, remove the first two:
undetrended_data <- undetrended_data[-c(1:(max(which(undetrended_data$undetrended == max(undetrended_data$undetrended)))+1000)), ]  


########################################################################    

#             Step 2 : Create splines, find W, U, V and O              #

######################################################################## 


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

# Find U and V:
u_v <- find_u_v(dat = undetrended_data$undetrended, wx = w$w_poly_peaks, wy = w$w_poly_peaks_yval, d1 = deriv1, d1p = deriv1_poly, spline = spline_poly, spline_o = spline_1, plot=FALSE)

## Find o in order to find the baseline
o <- c()
for(i in 1:length(w$w_poly_peaks)){
  o[i] <- max(which(inflexion_points < w$w_poly_peaks[i]))
}
#plot(spline_poly)
#points(inflexion_points[o], inflexion_points_yval[o], pch = 19)

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

# Refind W:
w_bc <- find_w(d1p = deriv1_poly_bc, deriv1 = deriv1_bc, sp = spline_poly_bc)

# Refind U and V:
u_v_bc <- find_u_v(dat = baseline_corrected, wx = w_bc$w_poly_peaks, wy = w_bc$w_poly_peaks_yval, d1 = deriv1_bc, d1p = deriv1_poly_bc, spline = spline_poly_bc, spline_o = spline_1_bc, plot = FALSE)

# Find v-u differences:
v_minus_u <- u_v_bc$v_yval - u_v_bc$u_yval

## Find o points again:
o_yval <- predict(spline_poly_bc, inflexion_points[o])
#plot(spline_poly_bc)
#points(inflexion_points[o], o_yval, pch = 19)

# Find o-w difference:
o_w_difference <- c()
for(i in 1:length(w_bc$w_poly_peaks)){
  o_w_difference[i] <- w_bc$w_poly_peaks[i] - inflexion_points[o[i]]
}

# Find distance between o_points:
o_difference <- c()
for(i in 1:(length(inflexion_points[o])-1)){
  o_difference[i] <- inflexion_points[o[i+1]] - inflexion_points[o[i]]
}


########################################################################    

#       Step 4 : Find individual waves and the average wave            #

######################################################################## 

# Find the average length of a wave (and 15 since we are starting the wave from before O):
source_data_column_length <- round(mean(o_difference)+15) 

# Redefine baseline corrected data:
sourcedata <- baseline_corrected[1:length(undetrended_data$undetrended)]

# Define a dataframe to contain individual waves (first column is the x-axis (in seconds)):
pulse <- data.frame(seq((-141/750), ((source_data_column_length*10 -9)-142)/750, by = 1/750)) 

# Create a dataframe of waveforms
after_o <- list()
for(i in 1:(length(w_bc$w_poly_peaks)-1)){  
  
  # Make a polynomial spline of rounded u - 15 : rounded u + source_data_column_length - 10:   # Now row 141 in xxxx = 0, therefore u = 0 
  spline_poly_wave_subset <- CubicInterpSplineAsPiecePoly((round(u_v_bc$u[i])-15):(round(u_v_bc$u[i]) + (source_data_column_length-10)), sourcedata[(round(u_v_bc$u[i])-15):(round(u_v_bc$u[i]) + (source_data_column_length-10))], "natural")
  
  # Turn into discrete form
  xxxx <- predict(spline_poly_wave_subset, c(seq((u_v_bc$u[i]-14), (u_v_bc$u[i]+(source_data_column_length-15)), 0.1)))  
  
  # Make into dataframe:
  xxxx <-  as.data.frame(xxxx)
  xxxx <- cbind(xxxx, c(seq((u_v_bc$u[i]-14), (u_v_bc$u[i]+(source_data_column_length-15)), 0.1)))
  colnames(xxxx) <- c('y', 'x') 
  # Scale so that v-u = 1
  xxxx$y <- xxxx$y/(v_minus_u[i])     
  # Adjust such that u = 0, v = 1 on y-axis
  y_axis_difference <- xxxx$y[141]           
  xxxx$y <- xxxx$y - y_axis_difference
 
  # Find the x-value for each wave that corresponds to when it = 0.5 in height (this requires making a spline):
  spline_poly_wave_subset_2 <- CubicInterpSplineAsPiecePoly(xxxx$x, xxxx$y, "natural")
  half_crossing <- solve(spline_poly_wave_subset_2, b = 0.5, deriv = 0)
  half_crossing <- half_crossing[which(abs(half_crossing - w_bc$w_poly_peaks[i]) == min(abs(half_crossing - w_bc$w_poly_peaks[i])))]    
  
  # Convert to discrete form again: (need to redefine xxxx)
  xxxx2 <- predict(spline_poly_wave_subset, c(seq((half_crossing-14), (half_crossing+(source_data_column_length-15)), 0.1)))  
  xxxx2 <-  as.data.frame(xxxx2)
  xxxx2 <- cbind(xxxx2, c(seq((half_crossing-14), (half_crossing+(source_data_column_length-15)), 0.1)))  
  colnames(xxxx2) <- c('y', 'x') 
  
  # Scale again
  xxxx2$y <- xxxx2$y/(v_minus_u[i]) 
  # Adjust y-axis such that u = 0, v = 1
  y_axis_difference <- u_v_bc$u_yval[i] / v_minus_u[i]
  xxxx2$y <- xxxx2$y - y_axis_difference
  
  # Find next_o
  after_o[[i]] <- which(xxxx2$x > inflexion_points[o][i+1])
  
  # Correct such that x column and wave column are correctly aligned
  xxxx3 <- c()
  for(i in 1:nrow(xxxx2)){
    xxxx3[i+1] <- xxxx2$y[i]
  }
  xxxx3 <- xxxx3[-length(xxxx3)]
  
  # Add to column to dataframe
  pulse <- cbind(pulse, xxxx3)
}
for(i in 1:(length(w_bc$w_poly_peaks)-1)){ 
  colnames(pulse)[i+1] <- paste("wave", i, sep = "_")       
}
colnames(pulse)[1] <- "x"

# Remove any values after O for each wave:
for(i in 2:(ncol(pulse))){
  pulse[, i][after_o[[(i-1)]][-1]] <- NA  
}

# Find the average wave:
# First replace NA values with last value that was not NA i.e o, which will require a clone dataframe:
pulse_for_finding_average <- pulse
for(i in 2:(ncol(pulse_for_finding_average))){
  pulse_for_finding_average[, i][after_o[[(i-1)]][-1]] <- pulse_for_finding_average[, i][after_o[[(i-1)]][1]]  
}
average_wave <- c()
sd_wave <- c()
median_wave <- c()
for(i in 1:nrow(pulse_for_finding_average)){
  row_vector <- c()
  for(j in 2:(ncol(pulse_for_finding_average))){
    row_vector[j-1] <- pulse_for_finding_average[i, j] 
  }
  average_wave[i] <- mean(row_vector[!is.na(row_vector)])   
  sd_wave[i] <- sd(row_vector[!is.na(row_vector)]) 
  median_wave[i] <- median(row_vector[!is.na(row_vector)])
}

# Can now stack and plot the mean +/- median wave 

pulse_stacked <- gather(pulse, key = "wave_ID", value = "values", -c("x"))

average <- data.frame(seq((-141/750), ((source_data_column_length*10 -9)-142)/750, by = 1/750))  
average <- cbind(average, average_wave)
colnames(average)[1] <- "x"

ggplot(data = pulse_stacked, aes(x, values, col = wave_ID), col = "black") +
  scale_color_manual(values = rep("black", ncol(pulse))) +  
  geom_line(size = 1.5, alpha = ((1/length(w_bc$w_poly_peaks)*10)-(1/length(w_bc$w_poly_peaks))))  + ylim(range(average_wave[-1]*1.5)) + geom_line(data = average, aes(x, average_wave), size = 1.125, color = "red") +  # ylim will vary based on source data
  theme(legend.position = "none") + labs( y= "PPG Output", x = "Time (Seconds)") + xlim(range(pulse$x))

# Create a polynomial spline for each wave:
poly_wave <- list()
for(i in 2:ncol(pulse)){
  poly_wave[[i-1]] <-CubicInterpSplineAsPiecePoly(pulse$x, pulse[, i], "natural")
}


########################################################################    

#                     Step 5 : Find O, S, N and D                      #

########################################################################


#Use new polynomial splines to find w/u/v/notch values for each waveform 
wuv <- find_wuv(p=pulse, col_len = source_data_column_length, p_w = poly_wave)

## Find O, S, N, D on the new polynomial splines:
osnd_xy <- find_osnd(p = pulse, p_w = poly_wave, col_len = source_data_column_length, wuvn = wuv)
osnd_y <- osnd_xy[1:(length(osnd_xy)/2)]
osnd_x <- osnd_xy[(length(osnd_xy)/2+1):length(osnd_xy)]


#Fit S and D sines using the OSND points

sd_sines <- find_sd_sine(p = pulse, wuvn = wuv, osndx = osnd_x, osndy = osnd_y, pw = poly_wave, plot=FALSE)
s_sines <- sd_sines[1:(length(sd_sines)/2)]
d_sines <- sd_sines[(length(sd_sines)/2+1):length(sd_sines)]

#Fit N sines using the S and D sines

n_sines <- fit_n_sine(p=pulse, ss = s_sines, ds = d_sines, osndx = osnd_x, osndy = osnd_y, wuvn = wuv, plot=TRUE)


##Refitting the SND peaks 

avg_period <- 3*(mean(wuv$v_x)-mean(wuv$u_x))
refitted_snd <- refit_peaks(p=pulse, ss= s_sines, ds= d_sines, ns= n_sines, period= avg_period)

##Plotting the refitted SND peaks back on the waveforms 
#for(i in 2:(ncol(pulse2))){
#  plot(1:nrow(pulse2), pulse2[,i], type='l')
#  points(refitted_snd$s_x[i-1], refitted_snd$s_y[i-1], pch = 19, col ='red')
#  points(refitted_snd$n_x[i-1], refitted_snd$n_y[i-1], pch = 19, col = 'blue')
#  points(refitted_snd$d_x[i-1], refitted_snd$d_y[i-1], pch = 19, col = 'yellow')
#}

#Creating a new osnd_x/y list incorporating the reffited SND peaks 
refit_osnd_x <- list()
refit_osnd_y <- list()

for(i in 1:length(osnd_x)){
  refit_osnd_x[[i]] <- c(NA, NA, NA, NA)
  refit_osnd_x[[c(i,1)]] <- osnd_x[[c(i,1)]]
  refit_osnd_x[[c(i,2)]] <- refitted_snd$s_x[i]
  refit_osnd_x[[c(i,3)]] <- refitted_snd$n_x[i]
  refit_osnd_x[[c(i,4)]] <- refitted_snd$d_x[i]
  refit_osnd_y[[i]] <- c(NA, NA, NA, NA)
  refit_osnd_y[[c(i,1)]] <- osnd_y[[c(i,1)]]
  refit_osnd_y[[c(i,2)]] <- refitted_snd$s_y[i]
  refit_osnd_y[[c(i,3)]] <- refitted_snd$n_y[i]
  refit_osnd_y[[c(i,4)]] <- refitted_snd$d_y[i]
}



#Refit sines using the refitted SND values
refit_sd_sines <- find_sd_sine(p = pulse, wuvn = wuv, osndx = refit_osnd_x, osndy = refit_osnd_y, pw = poly_wave, plot=FALSE)
refit_s_sines <- refit_sd_sines[1:(length(sd_sines)/2)]
refit_d_sines <- refit_sd_sines[(length(sd_sines)/2+1):length(sd_sines)]

refit_n_sines <- fit_n_sine(p=pulse, ss = s_sines, ds = d_sines, osndx = refit_osnd_x, osndy = refit_osnd_y, wuvn = wuv, plot=FALSE)

#Finding the residual after subtracting all the sines from the original pulse 
sine_resid <- list()
const <- c()
for(i in 2:ncol(pulse)){
  const <- rep(refit_s_sines[[i-1]][[length(refit_s_sines[[i-1]])]],length(refit_s_sines[[i-1]]))
  sine_resid[[i-1]] <- pulse[,i]
  sine_resid[[i-1]] <- (((((sine_resid[[i-1]] - refit_s_sines[[i-1]]) + const) - refit_d_sines[[i-1]]) + const) - refit_n_sines[[i-1]]) + const
}

#Plotting all the sines on the pulse trace, as well as the residual
#for(i in 2:ncol(pulse2)){
#  plot(1:(nrow(pulse2)), pulse2[,i], type = 'l')
#  lines(1:(nrow(pulse2)), refit_s_sines[[i-1]], col='red')
#  lines(1:(nrow(pulse2)), refit_d_sines[[i-1]], col='purple')
#  lines(1:(nrow(pulse2)), refit_n_sines[[i-1]], col='green')
#  lines(1:(nrow(pulse2)), sine_resid[[i-1]], col='yellow')
#}


#Correct all OSNDs so that O = 0
for(i in 1:length(osnd_y)){
  osnd_correction <- osnd_y[[i]]
  osnd_correction <- osnd_correction - osnd_correction[1]
  osnd_y[[i]] <- osnd_correction
  print(osnd_y[[i]])
}




## Features (canonical waveform first):

# Peak to notch time (waveforms need to be normalized on o-o interval first (or can divide by the pulse duration)):
pn_time <- c()
for(i in 1:length(osnd)){
  pn_time[i] <- (x_osnd[[c(i, 3)]] - x_osnd[[c(i, 2)]])/(next_o[i] - x_osnd[[c(i, 1)]])
}

# Peak to peak time (PPT):
PPT <- c()
for(i in 1:length(osnd)){
  PPT[i] <- x_osnd[[c(i, 4)]] - x_osnd[[c(i, 2)]]
}

# Stiffness Index = height / PPT
SI <- c()
for(i in 1:length(osnd)){
  SI[i] <- # Insert height here # / PPT[i]
}

# Notch to peak ratio  = height of notch / height of primary peak (this was referred to as RI in multivariate paper)
np_ratio <- c()
for(i in 1:length(osnd)){
  np_ratio[i] <-  osnd[[c(i, 3)]] / osnd[[c(i, 2)]]
}

# Notch-time ratio = time interval from notch to end of pulse / time interval from notch to beginning of pulse
nt_ratio <- c()
for(i in 1:length(osnd)){
  nt_ratio[i] <- (next_o[i] - x_osnd[[c(i, 3)]]) /  (x_osnd[[c(i, 3)]] - x_osnd[[c(i, 1)]])
}

# Systolic amplitude (aka amplitude):
sa <- c()
for(i in 1:length(osnd)){
  sa[i] <-  osnd[[c(i, 2)]]
}

# Reflectance peak to forward peak ratio (Augmentation index as per Takazawa et al, Reflection index as per Padilla et al) )
AI <- c()
for(i in 1:length(osnd)){
  AI[i] <-  osnd[[c(i, 4)]] / osnd[[c(i, 2)]]
}

# Alternate augmentation index (Rubins et al) ( (x-y)/x ):
aAI <- c()
for(i in 1:length(osnd)){
  aAI[i] <- (osnd[[c(i, 2)]] - osnd[[c(i, 4)]]) / osnd[[c(i, 2)]]
}

# Crest time (time from foot of waveform (o) to peak (S)):
CT <- c()
for(i in 1:length(osnd)){
  CT[i] <- x_osnd[[c(i, 2)]] - x_osnd[[c(i, 1)]]
}

# Inflexion point area ratio (For canonical waveform use c(i, 3), if using inflection point then use c(i, 4)):
ipa_ratio <- c()
for(i in 1:(length(poly_wave)-1)){
  f <- c(x_osnd[[c(i, 1)]]:x_osnd[[c(i, 3)]], x_osnd[[c(i, 3)]])
  g <- predict(poly_wave[[i]], f)
  j <- c(x_osnd[[c(i, 3)]]:next_o[i], next_o[i])
  k <- predict(poly_wave[[i]], j)
  auc_systole <- AUC(f, g, method = "spline")
  if(sum(j < 0) > 0){
    poly_wave_diastole_subset <-CubicInterpSplineAsPiecePoly(j, k, "natural")
    zero_crossing <- solve(poly_wave_diastole_subset, b = 0) 
    zero_crossing_yval <- predict(poly_wave_diastole_subset, zero_crossing)
    first_zero_crossing <- zero_crossing[min(which(zero_crossing > (x_osnd[[c(i, 1)]] + 30)))]   # first element that crosses 0 on the waves descent
    j <- c(x_osnd[[c(i, 3)]]:first_zero_crossing, first_zero_crossing)
    k <- predict(poly_wave[[i]], j)
    auc_diastole <- AUC(j, k, method = "spline")
  }else{
    auc_diastole <- AUC(j, k, method = "spline")
  }
  ipa_ratio[i] <- auc_diastole / auc_systole
  print(auc_diastole + auc_systole)
}




##Spectral Analysis

spectralanalysis <- (spectrum(baseline_corrected))


########
