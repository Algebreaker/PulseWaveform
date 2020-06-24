setwd("/home/johanna/Documents/Ninja theory/PulseAnalysis/Data/Craig")
data <- read.table("Source2.csv", header=T, sep=",") #first line of the csv file needs to be deleted


source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/preproc.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/find_w.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/find_u_v.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/baseline.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/find_osnd.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/spectrum.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/find_wuv.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/fit_sd_sine.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/fit_n_sine.R")
source("/home/johanna/Documents/Ninja theory/PulseAnalysis/git/PulseWaveform/refit_peaks.R")


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

w <- find_w(dat=undetrended_data$undetrended, d1 = deriv1, d1p = deriv1_poly, sp = spline_poly, plot=FALSE)

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
o_yval <- predict(spline_poly, inflexion_points[o])
plot(spline_poly)
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

source_data_column_length_precursor <- c(71, 81, 91, 101, 111, 121)
new_vector <- which(abs(source_data_column_length_precursor - (mean(o_difference)+15)) == min(abs(source_data_column_length_precursor - (mean(o_difference)+15))))
if((mean(o_difference)+15) > source_data_column_length_precursor[new_vector]){
  source_data_column_length <- source_data_column_length_precursor[new_vector+1]
}else{
    source_data_column_length <- source_data_column_length_precursor[new_vector]
}


## Chopping up the original data_undetrended (now baseline_corrected) into individual waves:

## Remove incomplete waves from beginning / end 
if(w_bc$w_poly_peaks[1] < 15){
  w_bc  <- w_bc[-1,]      # First w peak should be greater than 10
  o_w_difference <- o_w_difference[-1]
  o_difference <- o_difference[-1]
}
#if the last element is less than length of baseline corrected
if(w_bc$w_poly_peaks[length(w_bc$w_poly_peaks)] > (length(baseline_corrected) - 77)){
  w_bc <- w_bc[-length(w_bc$w_poly_peaks), ]
  o_w_difference <- o_w_difference[-length(o_w_difference)]
  o_difference <- o_difference[-length(o_difference)]
}

## Create a stacked waveform first so that waves can be matched exactly according to w

sourcedata <- baseline_corrected[1:length(undetrended_data$undetrended)]

pulse <- data.frame()

# change this from 1:length(w_poly_peaks) to Q1 when averaging quartiles
for(i in 1:(length(w_bc$w_poly_peaks))){              
#### Something isn't working in this (see intermediate_poly_Wave[[21]]) source 2 data....
  spline_poly_wave_subset <- CubicInterpSplineAsPiecePoly((round(w_bc$w_poly_peaks[i])-15):(round(w_bc$w_poly_peaks[i]) +  (source_data_column_length-10)) , sourcedata[(round(w_bc$w_poly_peaks[i])-15):(round(w_bc$w_poly_peaks[i]) + (source_data_column_length-10))], "natural")
  
  # Now get the y-values to fill the dataframe using the predict function
  
  xxxx <- predict(spline_poly_wave_subset, c(seq((round(w_bc$w_poly_peaks[i])-15), 
                                                 (round(w_bc$w_poly_peaks[i])-1), 0.1),   # this is now 141 in length
                                             w_bc$w_poly_peaks[i],                        # this is now number 142
                                             seq((round(w_bc$w_poly_peaks[i])+1), 
                                                 (round(w_bc$w_poly_peaks[i])+(source_data_column_length-10)), 0.1)))  
  xxxx <-  as.data.frame(xxxx)
  xxxx <- cbind(xxxx, c(seq((round(w_bc$w_poly_peaks[i])-15), (round(w_bc$w_poly_peaks[i])-1), 0.1), w_bc$w_poly_peaks[i],  seq(   (round(w_bc$w_poly_peaks[i])+1), (round(w_bc$w_poly_peaks[i])+(source_data_column_length-10)), 0.1)))
  colnames(xxxx) <- c('y', 'x') 
  xxxx[, 2] <- xxxx[, 2] - (xxxx$x[1]-1) 
  xxxx$wave <- i 
  # Need to scale so that v-u = 1
  xxxx$y <- xxxx$y/(v_minus_u[i])      ## scaling is an issue - did you calculate u and v?
  # Need to calculate difference between this w point and the first w point
  y_axis_difference <- w_bc$w_poly_peaks_yval[1] - xxxx$y[142]
  xxxx$y <- xxxx$y + y_axis_difference
  # Need to adjust x values so that all w's line up on x-axis
  x_axis_difference <- w_bc$w_poly_peaks[1] - xxxx$x[142]
  xxxx$x <- xxxx$x + x_axis_difference
  # Adjust such that w = 0 on x-axis
  xxxx$x <- xxxx$x - xxxx$x[142]   # 92nd element of each wave is w
  
  pulse <- rbind(pulse, xxxx)
}


## Create intermediate splines for converting to non-stacked data-frame:     
intermediate_poly_wave <- list()
for(i in 1:(length(w_bc$w_poly_peaks))){             
  intermediate_poly_wave[[i]] <- CubicInterpSplineAsPiecePoly(pulse$x[(((i-1)*(
    ((source_data_column_length*10)+33)))+1):((((i-1)*(
      ((source_data_column_length*10)+33)))+1)+(
        ((source_data_column_length*10)+32)))], pulse$y[(((i-1)*(
          ((source_data_column_length*10)+33)))+1):((((i-1)*(
            ((source_data_column_length*10)+33)))+1)+
              ((source_data_column_length*10)+32))], "natural")  
}

pulse2 <- data.frame(-14:(source_data_column_length-11))      ## -9 isn't far back enough - o's are being missed... is - 14?
for(i in 1:(length(w_bc$w_poly_peaks))){
  
  yval <- predict(intermediate_poly_wave[[i]], c(-14:(source_data_column_length-11)))   ## Since w = 0 on x-axis, you can specify relative to 0 
  pulse2 <- cbind(pulse2, yval)
  colnames(pulse2)[i+1] <- paste("wave", i, sep = "_") 
  
}
colnames(pulse2)[1] <- 'x'

# adjust such that w = 0.5 on y-axis
pulse2[, -1] <- pulse2[, -1] - pulse2$wave_1[15] + 0.5


## Now find average waveform by averaging each row in pulse2 + additional variance parameters:
average_wave <- c()
sd_wave <- c()
median_wave <- c()
for(i in 1:nrow(pulse2)){
  row_vector <- c()
  for(j in 2:(ncol(pulse2))){
    row_vector[j-1] <- pulse2[i, j] 
  }
  average_wave[i] <- mean(row_vector)
  sd_wave[i] <- sd(row_vector) 
  median_wave[i] <- median(row_vector)
}


## Can now stack and plot the mean +/- median wave 

pulse_stacked <- gather(pulse2, key = "wave_ID", value = "values", -c("x"))

average <- data.frame(-14:(source_data_column_length-11))    ## should differ depending on length of waveform
average <- cbind(average, average_wave)
colnames(average)[1] <- "x"

median <- data.frame(-14:(source_data_column_length-11))
median <- cbind(median, median_wave)
colnames(median)[1] <- "x"

ggplot(data = pulse_stacked, aes(x, values, col = wave_ID), col = "black") +
  scale_color_manual(values = rep("black", ncol(pulse2))) +  
  geom_line(size = 1.5, alpha = ((1/length(w_bc$w_poly_peaks)*10)-(1/length(w_bc$w_poly_peaks)))) + geom_line(data = average, aes(x, average_wave), size = 1.125, color = "red") + ylim(-0.5, 1.75) +   # ylim will vary based on source data
  theme(legend.position = "none")


## Now data is scaled and y-axis normalized*, create a new polynomial spline for each wave

poly_wave <- list()
for(i in 2:ncol(pulse2)){
  poly_wave[[i-1]] <-CubicInterpSplineAsPiecePoly(1:length(pulse2[, i]), pulse2[, i], "natural")
}


#Use new polynomial splines to find w/u/v/notch values for each waveform 
wuv <- find_wuv(p=pulse2, col_len = source_data_column_length, p_w = poly_wave)

## Find O, S, N, D on the new polynomial splines:

osnd_xy <- find_osnd(p = pulse2, p_w = poly_wave, col_len = source_data_column_length, wuvn = wuv)
osnd_y <- osnd_xy[1:(length(osnd_xy)/2)]
osnd_x <- osnd_xy[(length(osnd_xy)/2+1):length(osnd_xy)]


#Fit S and D sines using the OSND points

sd_sines <- find_sd_sine(p = pulse2, wuvn = wuv, osndx = osnd_x, osndy = osnd_y, pw = poly_wave, plot=FALSE)
s_sines <- sd_sines[1:(length(sd_sines)/2)]
d_sines <- sd_sines[(length(sd_sines)/2+1):length(sd_sines)]

#Fit N sines using the S and D sines

n_sines <- fit_n_sine(p=pulse2, ss = s_sines, ds = d_sines, osndx = osnd_x, osndy = osnd_y, wuvn = wuv, plot=FALSE)


##Refitting the SND peaks 

avg_period <- 3*(mean(wuv$v_x)-mean(wuv$u_x))
refitted_snd <- refit_peaks(p=pulse2, ss= s_sines, ds= d_sines, ns= n_sines, period= avg_period)

##Plotting the refitted SND peaks back on the waveforms 
for(i in 2:(ncol(pulse2))){
  plot(1:nrow(pulse2), pulse2[,i], type='l')
  points(refitted_snd$s_x[i-1], refitted_snd$s_y[i-1], pch = 19, col ='red')
  points(refitted_snd$n_x[i-1], refitted_snd$n_y[i-1], pch = 19, col = 'blue')
  points(refitted_snd$d_x[i-1], refitted_snd$d_y[i-1], pch = 19, col = 'yellow')
}
  


#Print all osnd's
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

