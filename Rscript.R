setwd("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/PulseAnalysis/Data/Craig")
data <- read.table("Source.csv", header=T, sep=",")
library(tidyverse)                                 
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils) #SplinesUtils is best downloaded directly from Github
library(spectral)
library(seewave)

source("preproc.R")
source("osnd.R")
source("spectrum.R")
source("wuv.R")
source("sines.R")
source("refit_peaks.R")
source("Find_W_Revised.R")


#Preprocessing which involves downsampling data and undetrending 
undetrended_data <- data.frame(preproc(dat=data))

sfunction <- splinefun(1:length(undetrended_data$undetrended), undetrended_data$undetrended[1:length(undetrended_data$undetrended)], method = "natural")
deriv1 <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 1)

#Creating polynomial splines from the data 
spline_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), undetrended_data$undetrended[1:length(undetrended_data$undetrended)], "natural")
deriv1_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), deriv1, "natural") 
## Finding inflexion points on spline_poly (points on deriv1 where x = 0):
inflexion_points <- solve(spline_poly, b = 0, deriv = 1)
inflexion_points_yval <- predict(spline_poly, inflexion_points)

w <- find_w_r(dat=undetrended_data$undetrended, d1p = deriv1_poly)

u_v <- find_u_v(dat = undetrended_data$undetrended, wx = w$w_poly_peaks, wy = w$w_poly_peaks_yval, d1 = deriv1, d1p = deriv1_poly, plot=FALSE)

## Find o in order to find the baseline
o <- c()
for(i in 1:length(w$w_poly_peaks)){
  o[i] <- max(which(inflexion_points < w$w_poly_peaks[i]))
}
#plot(spline_poly)
#points(inflexion_points[o], inflexion_points_yval[o], pch = 19)


# Adjust for early O points: #Simon will figure this bit out 
for(i in 1:length(w$w_poly_peaks)){
  o_decider <- w$w_poly_peaks[i] - 2*(w$w_poly_peaks[i] - u_v$u[i])
  if(abs(o_decider - inflexion_points[o[i]]) > 1.5){
    inflexion_points[o[i]] <- o_decider
    inflexion_points_yval[o[i]] <- predict(spline_poly, o_decider)
  }
}
#plot(spline_poly)
#points(inflexion_points[o], inflexion_points_yval[o], pch = 19)

baseline_corrected <- baseline(plot=FALSE)

# Redefine splines now that baseline corrected:
sfunction_bc <- splinefun(1:length(baseline_corrected), baseline_corrected, method = "natural")
deriv1_bc <- sfunction_bc(seq(1, length(baseline_corrected)), deriv = 1)
deriv2_bc <- sfunction_bc(seq(1, length(baseline_corrected)), deriv = 2)

# Turning the baseline_corrected data into a piece-wise polynomial spline (non-discrete): 
spline_poly_bc <- CubicInterpSplineAsPiecePoly(1:length(baseline_corrected), baseline_corrected, "natural")
deriv1_poly_bc <- CubicInterpSplineAsPiecePoly(1:length(baseline_corrected), deriv1_bc, "natural") 


w_bc <- find_w(dat=baseline_corrected, d1 = deriv1_bc, d1p = deriv1_poly_bc, sp = spline_poly_bc, plot = FALSE)

u_v_bc <- find_u_v(dat = baseline_corrected, wx = w_bc$w_poly_peaks, wy = w_bc$w_poly_peaks_yval, d1 = deriv1_bc, d1p = deriv1_poly_bc, plot = FALSE)


# Find v-u differences
v_minus_u <- u_v_bc$v_yval - u_v_bc$u_yval

## Find o points again:
o_yval <- predict(spline_poly_bc, inflexion_points[o])
plot(spline_poly_bc)
points(inflexion_points[o], o_yval, pch = 19)

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

source_data_column_length <- round(max(o_difference)+15)   # Add 15 since we are starting the wave from before O


## Chopping up the original data_undetrended (now baseline_corrected) into individual waves:

## Remove incomplete waves from beginning / end 
if(w_bc$w_poly_peaks[1] < 15){
  w_bc  <- w_bc[-1,]      # First w peak should be greater than 10
  o_w_difference <- o_w_difference[-1]
  o_difference <- o_difference[-1]
  v_minus_u <- v_minus_u[-1]
}
#if the last element is less than length of baseline corrected
if(w_bc$w_poly_peaks[length(w_bc$w_poly_peaks)] > (length(baseline_corrected) - 77)){
  w_bc <- w_bc[-length(w_bc$w_poly_peaks), ]
  o_w_difference <- o_w_difference[-length(o_w_difference)]
  o_difference <- o_difference[-length(o_difference)]
  v_minus_u <- v_minus_u[-length(v_minus_u)]
}

## Create a dataframe of waveforms

sourcedata <- baseline_corrected[1:length(undetrended_data$undetrended)]

pulse <- data.frame(seq((-141/750), ((source_data_column_length*10 -9)-142)/750, by = 1/750))  # First column of pulse is the x-axis (in seconds)

# change this from 1:length(w_poly_peaks) to Q1 when averaging quartiles
for(i in 1:(length(w_bc$w_poly_peaks))){              
  
  # Make a polynomial spline of rounded w - 15 : rounded w + source_data_column_length - 10
  spline_poly_wave_subset <- CubicInterpSplineAsPiecePoly((round(w_bc$w_poly_peaks[i])-15):(round(w_bc$w_poly_peaks[i]) + (source_data_column_length-10)), sourcedata[(round(w_bc$w_poly_peaks[i])-15):(round(w_bc$w_poly_peaks[i]) + (source_data_column_length-10))], "natural")
  
  # Turn into discrete form
  xxxx <- predict(spline_poly_wave_subset, c(seq((w_bc$w_poly_peaks[i]-14), (w_bc$w_poly_peaks[i]+(source_data_column_length-15)), 0.1)))  
  
  xxxx <-  as.data.frame(xxxx)
  xxxx <- cbind(xxxx, c(seq(c(seq((w_bc$w_poly_peaks[i]-14), (w_bc$w_poly_peaks[i]+(source_data_column_length-15)), 0.1)))))
  colnames(xxxx) <- c('y', 'x') 
  xxxx[, 2] <- xxxx[, 2] - (xxxx$x[1]-1) 
  xxxx$wave <- i 
  # Need to scale so that v-u = 1
  xxxx$y <- xxxx$y/(v_minus_u[i])     
  # Need to calculate difference between this w point and the first w point
  y_axis_difference <- w_bc$w_poly_peaks_yval[1] - xxxx$y[142]
  xxxx$y <- xxxx$y + y_axis_difference
  # Need to adjust x values so that all w's line up on x-axis
  x_axis_difference <- w_bc$w_poly_peaks[1] - xxxx$x[142]
  xxxx$x <- xxxx$x + x_axis_difference
  # Adjust such that w = 0 on x-axis
  xxxx$x <- xxxx$x - xxxx$x[142]  
  # Adjust such that w = 0.5 on y-axis
  y_axis_difference <- xxxx$y[142] - 0.5
  xxxx$y <- xxxx$y - y_axis_difference
  
  pulse <- cbind(pulse, xxxx$y)
}
for(i in 1:(length(w_bc$w_poly_peaks))){ 
  colnames(pulse)[i+1] <- paste("wave", i, sep = "_")       
}
colnames(pulse)[1] <- "x"
                       

## Now find average waveform by averaging each row in pulse + additional variance parameters:
average_wave <- c()
sd_wave <- c()
median_wave <- c()
for(i in 1:nrow(pulse)){
  row_vector <- c()
  for(j in 2:(ncol(pulse))){
    row_vector[j-1] <- pulse[i, j] 
  }
  average_wave[i] <- mean(row_vector)
  sd_wave[i] <- sd(row_vector) 
  median_wave[i] <- median(row_vector)
}


## Can now stack and plot the mean +/- median wave 

pulse_stacked <- gather(pulse, key = "wave_ID", value = "values", -c("x"))

average <- data.frame(seq((-141/750), ((source_data_column_length*10 -9)-142)/750, by = 1/750))  
average <- cbind(average, average_wave)
colnames(average)[1] <- "x"

ggplot(data = pulse_stacked, aes(x, values, col = wave_ID), col = "black") +
  scale_color_manual(values = rep("black", ncol(pulse))) +  
  geom_line(size = 1.5, alpha = ((1/length(w_bc$w_poly_peaks)*10)-(1/length(w_bc$w_poly_peaks)))) + geom_line(data = average, aes(x, average_wave), size = 1.125, color = "red") + ylim(-0.5, 1.75) +   # ylim will vary based on source data
  theme(legend.position = "none") + labs( y= "PPG Output", x = "Time (Seconds)")


## Now data is scaled and y-axis normalized*, create a new polynomial spline for each wave

poly_wave <- list()
for(i in 2:ncol(pulse)){
  poly_wave[[i-1]] <-CubicInterpSplineAsPiecePoly(pulse$x, pulse[, i], "natural")
}


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
