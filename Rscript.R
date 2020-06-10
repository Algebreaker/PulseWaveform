setwd("/home/johanna/Documents/Ninja theory/PulseAnalysis/Data/Craig")
data <- read.table("Source2.csv", header=T, sep=",") #first line of the csv file needs to be deleted

library(tidyverse)                                  #Have tidyverse packages installed and call tidyverse in library()
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils)
library(spectral)
library(seewave)

#remove NaNs
data<-data[!(data$PPG.PulseOx1=='NaN'),]


#Downsample
# The BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second, 
# so we typically have repeated values in a pattern of 3-3-4 repeating. DownSample tries to retrieve the 
# unique values, with an effort to be robust against variation in the repeat pattern and also against 
# genuine repeated values.

list<-rle(data$PPG.PulseOx1)
ID <- rep(1:length(list$values), times = list$lengths)
data2 <- cbind(data, ID)
data_downsampled <-c()

for (i in 1:max(ID)){
  sub.data <- filter(data2, ID == i)
  if(nrow(sub.data) <= 4){
    data_downsampled <- rbind(data_downsampled, sub.data[1,])
  }else if(nrow(sub.data) > 4 ){data_downsampled <- rbind(data_downsampled, sub.data[1,], sub.data[5,])}
}

#Undetrend
# Analysis of device output indicates that the PPG signal is detrended by application of the following
# formula: OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1]), where the constance 0.96875 is
# an approximation fitted to the data.
# Individual pulse events are more comprehensible if the detrending is not used, so this function 
#removes it by inverting the above function. 
   # can simply reverse the equatino as in the c++ script. 
undetrended <-replicate(length(data_downsampled$PPG.PulseOx1)-1, 0) 
undetrended<-c(data_downsampled$PPG.PulseOx1[1],undetrended) #add first detrended value to vector
for (i in 2:length(data_downsampled$PPG.PulseOx1))   
{
  undetrended[i]<-((data_downsampled$PPG.PulseOx1[i]-80) - ((data_downsampled$PPG.PulseOx1[i-1]-80) * 0.96875) + (undetrended[i-1]))
}
undetrended_data<-cbind(data_downsampled,undetrended)


## create spline + derivatives 
#to get y value from x value just use sfunction(x)
sfunction <- splinefun(1:length(undetrended_data$undetrended), undetrended_data$undetrended[1:length(undetrended_data$undetrended)], method = "natural")
spline <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 0)
deriv1 <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 1)
deriv2 <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 2)

########## 

# Turning the undetrended data into a piece-wise polynomial spline (non-discrete): 
spline_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), undetrended_data$undetrended[1:length(undetrended_data$undetrended)], "natural")


## Finding inflexion points on spline_poly
# Find all x values for inflexion points (points on deriv1 that equal 0):
inflexion_points <- solve(spline_poly, b = 0, deriv = 1)
# Find the y values for inflexion points:
inflexion_points_yval <- predict(spline_poly, inflexion_points)
# Plot the y values:
plot(spline_poly)
points(inflexion_points, inflexion_points_yval, pch = 19)

## Finding W  (note that the quantile threshold 0.95 needs adjusting for some datasets)

# Finding inflexion points on deriv1 requires redefining 1st deriv as a piece-wise spline
deriv1_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), deriv1, "natural") 
# Find inflexion points on deriv1_poly
inflexion_points_deriv1 <- solve(deriv1_poly, b = 0, deriv = 1)
inflexion_points_deriv1_yval <- predict(deriv1_poly, inflexion_points_deriv1)
plot(deriv1_poly)
points(inflexion_points_deriv1, inflexion_points_deriv1_yval, pch = 19)
# Find correct threshold using histogram
# hdat<-hist(deriv1)
quantiles<-quantile(deriv1,probs=c(.025,.95))      ## this still needs adjusting based on source data              
threshold<-quantiles[2]
# Identifying peaks of deriv1_poly:
w_poly_peaks <- predict(deriv1_poly, inflexion_points_deriv1)
w_poly_peaks <- which(w_poly_peaks > threshold)     
w_poly_peaks <- inflexion_points_deriv1[w_poly_peaks]
# Plot on deriv1_poly
plot(deriv1_poly)
w_poly_peaks_yval <- predict(deriv1_poly, w_poly_peaks)
points(w_poly_peaks, w_poly_peaks_yval, pch = 19)
# Plot back on spline_poly:
plot(spline_poly)
w_poly_peaks_yval <- predict(spline_poly, w_poly_peaks)
points(w_poly_peaks, w_poly_peaks_yval, pch = 19)

## Find o in order to find the baseline
o <- c()
for(i in 1:length(w_poly_peaks)){
  o[i] <- max(which(inflexion_points < w_poly_peaks[i]))
}
plot(spline_poly)
points(inflexion_points[o], inflexion_points_yval[o], pch = 19)


## Finding U and V

# Find half the height of w (on derivative y-axis)
w_half_height <- predict(deriv1_poly, w_poly_peaks)/2
# Find u and v:
half_heights <- c()
half_heights_yval <- c()
for(i in 1:length(w_half_height)){
  deriv1_poly_peak_subset <- CubicInterpSplineAsPiecePoly((w_poly_peaks[i]-10):(w_poly_peaks[i]+10), deriv1[(w_poly_peaks[i]-10):(w_poly_peaks[i]+10)], "natural") 
  half_heights_precursor <- solve(deriv1_poly_peak_subset, b = w_half_height[i])
  half_heights[c((2*(i)-1), (2*(i)))] <- half_heights_precursor
  half_heights_yval[c((2*(i)-1), (2*(i)))] <- predict(deriv1_poly_peak_subset, half_heights[c((2*(i)-1), (2*(i)))])
}
# Plot u's and v's on deriv1_poly
plot(deriv1_poly)
points(half_heights, half_heights_yval, pch = 19)
# Find u and v 
u <- half_heights[seq_along(half_heights) %%2 != 0] 
v <- half_heights[seq_along(half_heights) %%2 == 0]   
# Find u and v y-values for spline_poly
u_v_yval <- c()
for(i in 1:length(w_poly_peaks)){
  spline_poly_peak_subset <- CubicInterpSplineAsPiecePoly((w_poly_peaks[i]-10):(w_poly_peaks[i]+10), spline[(w_poly_peaks[i]-10):(w_poly_peaks[i]+10)], "natural") 
  u_v_yval[c((2*(i)-1), (2*(i)))] <- predict(spline_poly_peak_subset, half_heights[c((2*(i)-1), (2*(i)))])
}
# Plot u's and v's on spline_poly
plot(spline_poly)
points(half_heights, u_v_yval, pch = 19)


# Adjust for early O points:
for(i in 1:length(w_poly_peaks)){
  o_decider <- w_poly_peaks[i] - 2*(w_poly_peaks[i] - u[i])
  if(abs(o_decider - inflexion_points[o[i]]) > 1.5){
    inflexion_points[o[i]] <- o_decider
    inflexion_points_yval[o[i]] <- predict(spline_poly, o_decider)
  }
}
plot(spline_poly)
points(inflexion_points[o], inflexion_points_yval[o], pch = 19)


# Making a (non-polynomial) spline to fit the baseline
sfunction3 <- splinefun(inflexion_points[o], inflexion_points_yval[o], method = "natural")
spline_base <- sfunction3(seq(1, length(undetrended_data$undetrended)), deriv = 0)

# Plotting spline_base on spline_poly
plot(spline_poly)
points(inflexion_points[o], inflexion_points_yval[o], pch = 19)
lines(spline_base)

# Correcting for baseline:
baseline_corrected <- undetrended_data$undetrended[1:length(undetrended_data$undetrended)] - spline_base
plot(baseline_corrected, type = "l")
# Plot new baseline (y = 0)
lines(1:length(undetrended_data$undetrended), seq(from = 0, to = 0, length.out = length(undetrended_data$undetrended)))

# Redefine splines now that baseline corrected:
sfunction <- splinefun(1:length(undetrended_data$undetrended), baseline_corrected, method = "natural")
spline <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 0)
deriv1 <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 1)
deriv2 <- sfunction(seq(1, length(undetrended_data$undetrended)), deriv = 2)

# Turning the baseline_corrected data into a piece-wise polynomial spline (non-discrete): 
spline_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), baseline_corrected[1:length(undetrended_data$undetrended)], "natural")

# then re-find W for corrected baseline
# Finding inflexion points on deriv1 requires redefining 1st deriv as a piece-meal spline
deriv1_poly <- CubicInterpSplineAsPiecePoly(1:length(undetrended_data$undetrended), deriv1, "natural") 
# Find inflexion points on deriv1_poly
inflexion_points_deriv1 <- solve(deriv1_poly, b = 0, deriv = 1)
inflexion_points_deriv1_yval <- predict(deriv1_poly, inflexion_points_deriv1)
plot(deriv1_poly)
points(inflexion_points_deriv1, inflexion_points_deriv1_yval, pch = 19)
# Find correct threshold using histogram
# hdat<-hist(deriv1)
quantiles<-quantile(deriv1,probs=c(.025,.95))
threshold<-quantiles[2]
# Identifying peaks of deriv1_poly:
w_poly_peaks <- predict(deriv1_poly, inflexion_points_deriv1)
w_poly_peaks <- which(w_poly_peaks > threshold)     
w_poly_peaks <- inflexion_points_deriv1[w_poly_peaks]
# Plot on deriv1_poly
plot(deriv1_poly)
w_poly_peaks_yval <- predict(deriv1_poly, w_poly_peaks)
points(w_poly_peaks, w_poly_peaks_yval, pch = 19)
# Plot back on spline_poly:  ## this only exists in the old form
plot(spline_poly)
w_poly_peaks_yval <- predict(spline_poly, w_poly_peaks)
points(w_poly_peaks, w_poly_peaks_yval, pch = 19)


##Find U and V again

# Find half the height of w (on derivative y-axis)
w_half_height <- predict(deriv1_poly, w_poly_peaks)/2
# Find u and v:
half_heights <- c()
half_heights_yval <- c()
for(i in 1:length(w_half_height)){
  deriv1_poly_peak_subset <- CubicInterpSplineAsPiecePoly((w_poly_peaks[i]-10):(w_poly_peaks[i]+10), deriv1[(w_poly_peaks[i]-10):(w_poly_peaks[i]+10)], "natural") 
  half_heights_precursor <- solve(deriv1_poly_peak_subset, b = w_half_height[i])
  half_heights[c((2*(i)-1), (2*(i)))] <- half_heights_precursor
  half_heights_yval[c((2*(i)-1), (2*(i)))] <- predict(deriv1_poly_peak_subset, half_heights[c((2*(i)-1), (2*(i)))])
}
# Plot u's and v's on deriv1_poly
plot(deriv1_poly)
points(half_heights, half_heights_yval, pch = 19)
# Find u and v 
u <- half_heights[seq_along(half_heights) %%2 != 0] 
v <- half_heights[seq_along(half_heights) %%2 == 0]   
# Find u and v y-values for spline_poly
u_v_yval <- c()
for(i in 1:length(w_poly_peaks)){
  spline_poly_peak_subset <- CubicInterpSplineAsPiecePoly((w_poly_peaks[i]-10):(w_poly_peaks[i]+10), spline[(w_poly_peaks[i]-10):(w_poly_peaks[i]+10)], "natural") 
  u_v_yval[c((2*(i)-1), (2*(i)))] <- predict(spline_poly_peak_subset, half_heights[c((2*(i)-1), (2*(i)))])
}
# Plot u's and v's on spline_poly
plot(spline_poly)
points(half_heights, u_v_yval, pch = 19)

## Finding a scalar for each wave:
# Find individual u and v y-values 
u_yval <- u_v_yval[seq_along(u_v_yval) %%2 != 0] 
v_yval <- u_v_yval[seq_along(u_v_yval) %%2 == 0]   
# Find v-u differences
v_minus_u <- v_yval - u_yval


## Find o points again:
o_yval <- predict(spline_poly, inflexion_points[o])
plot(spline_poly)
points(inflexion_points[o], o_yval, pch = 19)

# Find o-w difference:
o_w_difference <- c()
for(i in 1:length(w_poly_peaks)){
  o_w_difference[i] <- w_poly_peaks[i] - inflexion_points[o[i]]
}


# Find mean distance between o_points:
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
if(w_poly_peaks[1] < 15){
  w_poly_peaks <- w_poly_peaks[-1]      # First w peak should be greater than 10
  o_w_difference <- o_w_difference[-1]
  o_difference <- o_difference[-1]
}
if(w_poly_peaks[length(w_poly_peaks)] > (length(baseline_corrected) - 77)){
  w_poly_peaks <- w_poly_peaks[-length(w_poly_peaks)]
  o_w_difference <- o_w_difference[-length(o_w_difference)]
  o_difference <- o_difference[-length(o_difference)]
}

## Create a stacked waveform first so that waves can be matched exactly according to w

sourcedata <- baseline_corrected[1:length(undetrended_data$undetrended)]

pulse <- data.frame()

for(i in 1:(length(w_poly_peaks))){              # change this from 1:length(w_poly_peaks) to Q1 when averaging quartiles
  
  
  #### Something isn't working in this (see intermediate_poly_Wave[[21]]) source 2 data....
  spline_poly_wave_subset <- CubicInterpSplineAsPiecePoly(    (round(w_poly_peaks[i])-15):(round(w_poly_peaks[i]) +  (source_data_column_length-10)) ,     sourcedata[(round(w_poly_peaks[i])-15):(round(w_poly_peaks[i]) + (source_data_column_length-10))], "natural")
  
  # Now get the y-values to fill the dataframe using the predict function
  
  xxxx <- predict(spline_poly_wave_subset, c(seq((round(w_poly_peaks[i])-15), 
                                                 (round(w_poly_peaks[i])-1), 0.1),   # this is now 141 in length
                                             w_poly_peaks[i],                        # this is now number 142
                                             seq((round(w_poly_peaks[i])+1), 
                                                 (round(w_poly_peaks[i])+(source_data_column_length-10)), 0.1)))  
  xxxx <-  as.data.frame(xxxx)
  xxxx <- cbind(xxxx, c(     seq(   (round(w_poly_peaks[i])-15), (round(w_poly_peaks[i])-1), 0.1), w_poly_peaks[i],  seq(   (round(w_poly_peaks[i])+1), (round(w_poly_peaks[i])+(source_data_column_length-10)), 0.1)))
  colnames(xxxx) <- c('y', 'x') 
  xxxx[, 2] <- xxxx[, 2] - (xxxx$x[1]-1) 
  xxxx$wave <- i 
  # Need to scale so that v-u = 1
  xxxx$y <- xxxx$y/(v_minus_u[i])      ## scaling is an issue - did you calculate u and v?
  # Need to calculate difference between this w point and the first w point
  y_axis_difference <- w_poly_peaks_yval[1] - xxxx$y[142]
  xxxx$y <- xxxx$y + y_axis_difference
  # Need to adjust x values so that all w's line up on x-axis
  x_axis_difference <- w_poly_peaks[1] - xxxx$x[142]
  xxxx$x <- xxxx$x + x_axis_difference
  # Adjust such that w = 0 on x-axis
  xxxx$x <- xxxx$x - xxxx$x[142]   # 92nd element of each wave is w
  
  pulse <- rbind(pulse, xxxx)
}


## Create intermediate splines for converting to non-stacked data-frame:     
intermediate_poly_wave <- list()
for(i in 1:(length(w_poly_peaks))){             
  intermediate_poly_wave[[i]] <- CubicInterpSplineAsPiecePoly(pulse$x[(((i-1)*(
    ((source_data_column_length*10)+33)))+1):((((i-1)*(
      ((source_data_column_length*10)+33)))+1)+(
        ((source_data_column_length*10)+32)))], pulse$y[(((i-1)*(
          ((source_data_column_length*10)+33)))+1):((((i-1)*(
            ((source_data_column_length*10)+33)))+1)+
              ((source_data_column_length*10)+32))], "natural")  
}

pulse2 <- data.frame(-14:(source_data_column_length-11))      ## -9 isn't far back enough - o's are being missed... is - 14?
for(i in 1:(length(w_poly_peaks))){
  
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
  geom_line(size = 1.5, alpha = ((1/length(w_poly_peaks)*10)-(1/length(w_poly_peaks)))) + geom_line(data = average, aes(x, average_wave), size = 1.125, color = "red") + ylim(-0.5, 1.75) +   # ylim will vary based on source data
  theme(legend.position = "none")


## Now data is scaled and y-axis normalized*, create a new polynomial spline for each wave

poly_wave <- list()
for(i in 2:ncol(pulse2)){
  poly_wave[[i-1]] <-CubicInterpSplineAsPiecePoly(1:length(pulse2[, i]), pulse2[, i], "natural")
}






## Find W, O, S, N, D on the new polynomial splines:

osnd <- list()
pulse_width <- c()
x_osnd <- list()
s <- c()
s_yval <- c()
next_o <- c()
next_o_yval <- c()
x <- c(1:(nrow(pulse)))
d_x <- c(1:(nrow(pulse)))
s_sine <- list()
d_sine <- list()

for(i in 2:(ncol(pulse2))){        
  
  sfunction2 <- splinefun(1:(source_data_column_length+4), pulse2[, i], method = "natural")
  deriv1_wave <- sfunction2(seq(1, (source_data_column_length+4)), deriv = 1)
  deriv1_wave_poly <- CubicInterpSplineAsPiecePoly(1:(source_data_column_length+4), deriv1_wave, "natural") 
  
  
  # Find inflexion points on deriv1_wave_poly
  inflexion_points_deriv1_wave_poly <- solve(deriv1_wave_poly, b = 0, deriv = 1)
  inflexion_points_deriv1_wave_poly_yval <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)
  #plot(deriv1_wave_poly)
  #points(inflexion_points_deriv1_wave_poly, inflexion_points_deriv1_wave_poly_yval, pch = 19)
  
  # Find correct threshold using histogram
  # hdat_2<-hist(deriv1_wave)
  quantiles_2 <- quantile(deriv1_wave, probs=c(.025,.95))
  threshold_2 <- quantiles_2[2]
  
  # Identifying peaks of deriv1_poly:
  w_poly_peaks_wave <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)
  w_poly_peaks_wave <- which(w_poly_peaks_wave > threshold_2)     
  w_poly_peaks_wave <- inflexion_points_deriv1_wave_poly[w_poly_peaks_wave]
  

  # Finding the 'notch' (renal wave)    
  notch <- inflexion_points_deriv1_wave_poly[(which(inflexion_points_deriv1_wave_poly_yval == min(inflexion_points_deriv1_wave_poly_yval)))+1]
  notch_poly_yval <- predict(poly_wave[[i-1]], notch)
	
	
  #Find U and V
  
  # Find half the height of w (on derivative y-axis)
  w_half_height_wave <- predict(deriv1_wave_poly, w_poly_peaks_wave[1])/2
  # Find u and v for derivative:
  half_heights_wave_new <- solve(deriv1_wave_poly, b = w_half_height_wave[1])
  half_heights_wave_new_yval <- predict(deriv1_wave_poly, half_heights_wave_new)
  u <- half_heights_wave_new[1]
  v <- half_heights_wave_new[2]
  # Find u and v y-values for original wave:
  u_v_yval_wave <- predict(poly_wave[[i-1]], half_heights_wave_new)
  u_yval <- u_v_yval_wave[1]
  v_yval <- u_v_yval_wave[2]

  
  # Find OSND
  inflexion_points_new <- solve(poly_wave[[i-1]], b = 0, deriv = 1)
  inflexion_points_new_yval <- predict(poly_wave[[i-1]], inflexion_points_new)
  o <- w_poly_peaks_wave[1] - o_w_difference[i-1]
  o_yval <- predict(poly_wave[[i-1]], o)
  
  if(length(inflexion_points_new) >= 4){
  
  #Find the four inflexion points that are 1 to the left and 3 to the right of W
  old_o <- max(which(inflexion_points_new < w_poly_peaks_wave[1]))
  # Find x coords of OSND
  x_osnd[[i-1]] <- inflexion_points_new[old_o:(old_o+3)]
  # Replace old o with new o                 # the new o is already confirmed as the correct one to use from a previous for loop          
  x_osnd[[c(i-1, 1)]] <- o
  
  # Find new S
  s[i-1] <- w_poly_peaks_wave[1] + (2*(v - w_poly_peaks_wave[1]))
  s_yval[i-1] <- predict(poly_wave[[i-1]], s[i-1])
  
  # Continue to define OSND (divergence here between waveforms based on if Paul type or not)
  # if s - w is greater than new s - w, aka if its a Paul type:
  if((inflexion_points_new[old_o+1] - w_poly_peaks_wave[1]) > (s[i-1] - w_poly_peaks_wave[1])){
    
    # Remove false N and D values 
    x_osnd_precursor <- x_osnd[[i-1]]
    if(length(which(complete.cases(x_osnd_precursor) ==1)) != 4){
      x_osnd_precursor <- x_osnd_precursor[-(which(is.na(x_osnd_precursor)))]
    }
    if(length(x_osnd_precursor) == 2){
      x_osnd[[i-1]] <- x_osnd_precursor
    }else{
      false_points <- which(x_osnd_precursor > notch) 
      x_osnd_precursor <- x_osnd_precursor[-false_points]
      x_osnd[[i-1]] <- x_osnd_precursor
    }
    # Find y coords of OSND
    osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
    osnd[[c(i-1, 1)]] <- o_yval
    
    # Remove false values again
    osnd_precursor <- osnd[[i-1]]
    if(length(which(complete.cases(osnd_precursor) ==1)) != 4){
      osnd_precursor <- osnd_precursor[-(which(is.na(osnd_precursor)))]
    }
    if(length(osnd_precursor) == 2){
      osnd[[i-1]] <- osnd_precursor
    }else{
      osnd_precursor <- osnd_precursor[-false_points]
      osnd[[i-1]] <- osnd_precursor
    }
    
    osnd[[c(i-1, 3)]] <- osnd[[c(i-1, 2)]]
    osnd[[c(i-1, 2)]] <- s_yval[i-1]
    x_osnd[[c(i-1, 3)]] <- x_osnd[[c(i-1, 2)]]
    x_osnd[[c(i-1, 2)]] <- s[i-1]
    osnd[[c(i-1, 4)]] <- notch_poly_yval
    x_osnd[[c(i-1, 4)]] <- notch
    

  }else{
    osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
    osnd[[c(i-1, 1)]] <- o_yval
    #osnd[[i-1]] <- osnd[[i-1]] - inflexion_points_new_yval[o]
    
    if((x_osnd[[c(i-1, 4)]] - x_osnd[[c(i-1, 3)]]) < 3 & (osnd[[c(i-1, 3)]] / (osnd[[c(i-1, 2)]] - osnd[[c(i-1, 1)]])) > 0.5){         # these lines are for cases when the canonical waveform has a prominent
     osnd[[c(i-1, 3)]] <- inflexion_points_new_yval[old_o + 4]                                                                         # renal wave such that o+3 and o+4 are inflexion points between the 
     x_osnd[[c(i-1, 3)]] <- inflexion_points_new[old_o + 4]                                                                            # systolic peak and the notch. This shifts them along by two. 
     osnd[[c(i-1, 4)]] <- inflexion_points_new_yval[old_o + 5]                                                                         # A fixed threshold isn't ideal, but a threshold relative to O-S could work. 
     x_osnd[[c(i-1, 4)]] <- inflexion_points_new[old_o + 5]
    }
    
    if((x_osnd[[c(i-1, 4)]] - x_osnd[[c(i-1, 3)]]) < 3 &  x_osnd[[c(i-1, 4)]] < 30){
      osnd[[c(i-1, 3)]] <- inflexion_points_new_yval[old_o + 4]                             # This is a special case for when there are two notches apparent
      x_osnd[[c(i-1, 3)]] <- inflexion_points_new[old_o + 4]                                # between systolic and diastolic waves
      osnd[[c(i-1, 4)]] <- inflexion_points_new_yval[old_o + 5]                             
      x_osnd[[c(i-1, 4)]] <- inflexion_points_new[old_o + 5]
    }
    

    
    if(osnd[[c(i-1, 3)]] < 0 & (osnd[[c(i-1, 4)]] - osnd[[c(i-1, 3)]]) > 1){        # If N is less than 0 and D-N >1, this implies that N has risen above D such 
      # find points either side of the notch                                        # that neither are inflection points and hence the inflection points
      d <- notch + 2                                                                # of next O and S are incorrectly assigned to N and D. 
      n <- notch - 2                                                                # These lines correct for that by finding the inflection point instead, 
      d_yval <- predict(poly_wave[[i-1]], d)                                        # and taking values either side of it that approximate N and D. 
      n_yval <- predict(poly_wave[[i-1]], n) 
      # assign those points
      osnd[[c(i-1, 3)]] <- n_yval
      x_osnd[[c(i-1, 3)]] <- n
      osnd[[c(i-1, 4)]] <- d_yval
      x_osnd[[c(i-1, 4)]] <- d
    }
    
    if((x_osnd[[c(i-1, 3)]] - x_osnd[[c(i-1, 2)]]) < 5){     # despite this being on the canonical waveform side, these lines
      osnd[[c(i-1, 3)]] <- osnd[[c(i-1, 4)]]                 # are here because when an inflection point exists between systolic
      x_osnd[[c(i-1, 3)]] <- x_osnd[[c(i-1, 4)]]             # peaks on non-canonical waveforms, they fulfill the condition
      osnd[[c(i-1, 4)]] <- notch_poly_yval                   # to be processed on the canonical side. 
      x_osnd[[c(i-1, 4)]] <- notch
    }
  }    
  }else{
    
    #define old_o
    old_o <- max(which(inflexion_points_new < w_poly_peaks_wave[1]))
    #define osnd
    x_osnd[[i-1]] <- inflexion_points_new[old_o:(old_o+3)]
    osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
    # Find new S
    s[i-1] <- w_poly_peaks_wave[1] + (2*(v - w_poly_peaks_wave[1]))
    s_yval[i-1] <- predict(poly_wave[[i-1]], s[i-1])
    
    # Distinguish between canonical vs non-canonical
    if((inflexion_points_new[old_o+1] - w_poly_peaks_wave[1]) > (s[i-1] - w_poly_peaks_wave[1])){
      
      # Remove false N and D values 
      x_osnd_precursor <- x_osnd[[i-1]]
      if(length(which(complete.cases(x_osnd_precursor) ==1)) != 4){
        x_osnd_precursor <- x_osnd_precursor[-(which(is.na(x_osnd_precursor)))]
      }
      if(length(x_osnd_precursor) == 2){
        x_osnd[[i-1]] <- x_osnd_precursor
      }else{
        false_points <- which(x_osnd_precursor > notch) 
        x_osnd_precursor <- x_osnd_precursor[-false_points]
        x_osnd[[i-1]] <- x_osnd_precursor
      }
      # Find y coords of OSND
      osnd[[i-1]] <- inflexion_points_new_yval[old_o:(old_o+3)]
      osnd[[c(i-1, 1)]] <- o_yval
      
      # Remove false values again
      osnd_precursor <- osnd[[i-1]]
      if(length(which(complete.cases(osnd_precursor) ==1)) != 4){
        osnd_precursor <- osnd_precursor[-(which(is.na(osnd_precursor)))]
      }
      if(length(osnd_precursor) == 2){
        osnd[[i-1]] <- osnd_precursor
      }else{
        osnd_precursor <- osnd_precursor[-false_points]
        osnd[[i-1]] <- osnd_precursor
      }
      
      osnd[[c(i-1, 3)]] <- osnd[[c(i-1, 2)]]
      osnd[[c(i-1, 2)]] <- s_yval[i-1]
      x_osnd[[c(i-1, 3)]] <- x_osnd[[c(i-1, 2)]]
      x_osnd[[c(i-1, 2)]] <- s[i-1]
      osnd[[c(i-1, 4)]] <- notch_poly_yval
      x_osnd[[c(i-1, 4)]] <- notch
      
    }else{
      x_osnd[[c(i-1, 1)]] <- o                                 # Cases where there are fewer than 4 inflexion points identified 
      x_osnd[[c(i-1, 2)]] <- inflexion_points_new[2]           # can occur if N and D converge and the spline isn't long enough for
      x_osnd[[c(i-1, 3)]] <- notch                             # the next S to register as an inflexion point.
      x_osnd[[c(i-1, 4)]] <- notch
      osnd[[c(i-1, 1)]] <- o_yval
      osnd[[c(i-1, 2)]] <- inflexion_points_new_yval[2]
      osnd[[c(i-1, 3)]] <- notch_poly_yval
      osnd[[c(i-1, 4)]] <- notch_poly_yval
      
    }   
  }
	

  # Find pulse width = width at half the amplitude of the wave
  #max_amp_precursor <- osnd[[i-1]] # second element will give amplitude of the wave
  #half_amp <- max_amp_precursor[2]/2
  #half_amp <- half_amp - inflexion_points_new_yval[o] # correcting for the fact that the spline goes below 0
  # find all x values where y = half_amp
  #xval_width <- solve(poly_wave[[i-1]], b = half_amp)
  # Calculate width from the first two x-values:
  #pulse_width[i-1] <- xval_width[2] - xval_width[1]

  
  # Find the next o_point
  next_o <- x_osnd[[c(i-1, 1)]] + o_difference[i-1]
  next_o_yval <- predict(poly_wave[[i-1]], next_o) 
	
  #Define period as 3*v-u
  period <- 3*(v-u)
  
  #Define phi and replace all values larger than pi with pi, and all values smaller than -pi with -pi. This ensures a single peak
  phi <- ((2*pi)/period)  * (x-(w_poly_peaks_wave[1]+(period/4)))
  phi[phi>(pi)] <- pi
  phi[phi<(-pi)] <- -pi
	
	
  # Finding S peak 
  # y = (S-O)/2 * cos(phi) + S/2
  y <-((osnd[[c(i-1, 2)]] -  osnd[[c(i-1, 1)]])/2) * cos(phi) + (osnd[[c(i-1, 2)]]/2) 
  #the last bit (s/2) only works if O is always 0 - otherwise should be (O+S)/2 
 
  
  #Defining phi for the D peak
  d_phi <- ((2*pi)/period)  * (d_x-(x_osnd[[c(i-1,4)]]))
  d_phi[d_phi>(pi)] <- pi
  d_phi[d_phi<(-pi)] <- -pi
	
  #Finding the sine for the D peak
  d_y <- ((osnd[[c(i-1,4)]]-osnd[[c(i-1, 1)]])/2) * cos(d_phi) + (osnd[[c(i-1,4)]]/2)

  
  
  # Plot back on poly_wave[[i]]:
  plot(poly_wave[[i-1]])
  w_poly_peaks_yval <- predict(poly_wave[[i-1]], w_poly_peaks_wave)
  points(w_poly_peaks_wave[1], w_poly_peaks_yval[1], pch = 19)
  points(x_osnd[[i-1]], osnd[[i-1]], pch = 19, col = "red")
  points(next_o, next_o_yval, pch = 19)
  lines(x, y, col = 'red')
  lines(d_x, d_y, col='purple')
  #points(half_heights_wave_new, u_v_yval_wave, pch = 19)
  #points(notch, notch_poly_yval, pch = 19)
  

  s_sine[[i-1]] <- y 
  d_sine[[i-1]] <- d_y
  
  
}

#Print all osnd's
for(i in 1:length(osnd)){
  osnd_correction <- osnd[[i]]
  osnd_correction <- osnd_correction - osnd_correction[1]
  osnd[[i]] <- osnd_correction
  print(osnd[[i]])
}

resid_test <- list()
resid_peaks_y <- c()
resid_peaks_x <- c()

#Find residual by subtracting the s_sine and the d_sine from the pulse data
#Find the x and y values of the peak of the residual to plot the N peak

for(i in 2:(ncol(pulse))){
  
  plot(pulse$x, pulse[,i], type = 'l')
  lines(1:(nrow(pulse)),s_sine[[i-1]], col='red')
  lines(1:(nrow(pulse)),d_sine[[i-1]], col='purple')
  d_y_resid <- pulse[,i] - s_sine[[i-1]]
  d_y_resid <- d_y_resid - d_sine[[i-1]]
  lines(1:(nrow(pulse)),d_y_resid, col='green')
  #points(x_osnd[[i-1]], osnd[[i-1]], pch = 19, col = "red")
  
  resid_test[[i-1]] <- d_y_resid
  
  resid_peaks_y[i-1] <- findpeaks(d_y_resid[20:40])[1]
  resid_peaks_x[i-1] <- 20+(findpeaks(d_y_resid[20:40])[2])
  points(resid_peaks_x[i-1], resid_peaks_y[i-1], pch = 19)
  
}

#Use residual to find peaks and plot the sine curve for the tidal wave/N peak 

n_x <- c(1:(nrow(pulse)))
n_sine <- list()
n_resid <- list()

for(i in 2:(ncol(pulse))){
  
  plot(pulse$x, resid_test[[i-1]], type='l')
  
  #define phi for notch
  n_phi <- ((2*pi)/period)  * (n_x-(resid_peaks_x[i-1]))
  n_phi[n_phi>(pi)] <- pi
  n_phi[n_phi<(-pi)] <- -pi
  
  #create notch sine curve
  n_y <- ((resid_peaks_y[i-1]-osnd[[c(i-1, 1)]])/2) * cos(n_phi) + (resid_peaks_y[i-1]/2)
  
  lines(1:(nrow(pulse)),n_y, col='green')
  n_sine[[i-1]] <- n_y
  
  resid <- resid_test[[i-1]] - n_y
  lines(pulse$x, resid, col='blue')
  n_resid[[i-1]] <- resid
}


####LUCIE'S PRELIMINARY REFITTING SPLINES BIT THAT ALMOST WORKS?


max_shift=4
sample_rate = 75
avg_period <- 3*(mean(v_vals)-mean(u_vals))
s_x <- c()
s_y <- c()
d_x <- c()
d_y <- c()
n_x <- c()
n_y  <- c()


##Takes in the current estimated values for the S, N and D curves and refits them using the trace derivatives 
or(i in 2:(ncol(pulse)-1)){

  fitted_s <- FALSE
  fitted_n <- FALSE
  fitted_d <- FALSE

  while(!fitted_s | !fitted_n | !fitted_d){
    
    ##define the peaks of the current waveform
    s_peak <- max(s_sine[[i-1]])
    n_peak <- max(n_sine[[i-1]])
    d_peak <- max(d_sine[[i-1]])
    d_peak_x <- match(max(d_sine[[i-1]]), d_sine[[i-1]])
    n_peak_x <- match(max(n_sine[[i-1]]), n_sine[[i-1]])
    s_peak_x <- match(max(s_sine[[i-1]]), s_sine[[i-1]])
    
    #if we've refitted any of the peaks already, define them as -1
    if(fitted_s){
        s_peak = -1
      }
    if(fitted_n){
        n_peak = -1
      }
    if(fitted_d){
        d_peak = -1
      }
    
    if(s_peak > n_peak){
      #subtract the d and n peaks from the original trace
      d_resid <- pulse[,i] - d_sine[[i-1]]
      d_n_resid <- d_resid - n_sine[[i-1]]
      
      #find peaks of the residual as the new S peak
      x_peaks <- findpeaks(d_n_resid)[,2]
      pk_loc <- which(abs(x_peaks-s_peak_x)==min(abs(x_peaks - s_peak_x)))
      s_x[i-1] <- findpeaks(d_n_resid)[pk_loc,2]
      s_y[i-1] <- findpeaks(d_n_resid)[pk_loc,1]

      
      fitted_s <- TRUE

    
    } else if(n_peak > d_peak){
    
      s_resid <- pulse[,i] - s_sine[[i-1]]
      s_d_resid <- s_resid - d_sine[[i-1]]
    
      s_d_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(pulse[, i]), s_d_resid, "natural")
      sfunction3 <- splinefun(1:source_data_column_length, s_d_resid, method = "natural")
      deriv1_sd <- sfunction3(seq(1, source_data_column_length), deriv = 1)
      deriv1_sd_poly <- CubicInterpSplineAsPiecePoly(1:source_data_column_length, deriv1_sd, "natural") 
      current_deriv_n <- predict(deriv1_sd_poly, n_peak_x)
    
        if(current_deriv_n < 0){
          n_x_est <- s_peak_x + findpeaks(deriv1_sd[s_peak_x:n_peak_x])[2]
          n_y_est <- s_d_resid[n_x_est]
        } else if(current_deriv_n > 0){
          n_x_est <- n_peak_x + findpeaks(deriv1_sd[n_peak_x:d_peak_x])[2]
          n_y_est <- s_d_resid[n_x_est]
        } 
        
        n_y[i-1] <- n_y_est
        n_x[i-1] <- n_x_est
        
        if(fitted_s & (!fitted_d)){
        fitted_s = FALSE;
        }

      fitted_n = TRUE;
    
    
    } else{
    
      s_resid <- pulse[,i] - s_sine[[i-1]]
      s_n_resid <- s_resid - n_sine[[i-1]]
      
      s_n_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(pulse[, i]), s_n_resid, "natural")
      sfunction4 <- splinefun(1:source_data_column_length, s_n_resid, method = "natural")
      deriv1_sn <- sfunction4(seq(1, source_data_column_length), deriv = 1)
      deriv1_sn_poly <- CubicInterpSplineAsPiecePoly(1:source_data_column_length, deriv1_sn, "natural") 
      d_y_est <- predict(deriv1_sn_poly, d_peak_x)

      
      if(d_y_est < 0){
        d_y_est <- findpeaks(deriv1_sn[(d_peak_x-1):length(s_n_resid)])[1]
        d_x_est <- d_peak_x + findpeaks(deriv1_sn[(d_peak_x-1):length(s_n_resid)])[2]
        } else if(d_y_est > 0){
        d_y_est <- findpeaks(deriv1_sn[(n_peak_x):d_peak_x])[1]
        d_x_est <- n_peak_x + findpeaks(deriv1_sn[n_peak_x:d_peak_x])[2]
        }
        
      if(d_x_est > (d_peak_x + max_shift)){
        max_deriv_d_x <- findpeaks(deriv1_sn[n_peak_x:length(deriv1_sn)])[1,2]
        max_deriv_d <- findpeaks(deriv1_sn[n_peak_x:length(deriv1_sn)])[1,1]
        d_y_est = 2 * max_deriv_d * sample_rate * avg_period
        d_x_est = max_deriv_d_x + 0.5 * pi * avg_period
      }
      
      d_y[i-1] <- d_y_est
      d_x[i-1] <- d_x_est
      
        
      fitted_d = TRUE
    
    }
  }
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

#filtered <- filter.fft(spline_poly,fc=0.0925,BW=0.0525,n=50)
#plot.fft(filtered)
#spectral<-spec.fft(spline_poly)
#plot(data_downsampled$PPG.PulseOx1)
#plot(baseline_corrected)

powerspectrum<-spectrum(data_downsampled$PPG.PulseOx1)
powerspectrum<-spectrum(baseline_corrected)

LF<-ffilter(data_downsampled$PPG.PulseOx1, f=75, from = 0.04, to = 0.145, bandpass = TRUE) #low bandpass, sampling frequency of 75 Hz (75 times per second)
powerspectrum<-spectrum(LF)
HF<-ffilter(data_downsampled$PPG.PulseOx1, f=75, from = 0.145, to = 0.45, bandpass = TRUE) #high bandpass
powerspectrum<-spectrum(HF)

spectralratio<-sum(LF)/sum(HF) #ratio of LF/HF

###sines sines sines 


pre_phase_test <- list()
post_phase_test <- list()


for(i in 2:((ncol(pulse))-1)){

  
  plot(pulse[,1], pulse[,i], type = 'l')

  #identify the different systolic variables to be entered into the sine eq
  #b is defined as 2*O-S
  syst_b <- (2*pi)/abs(2*(x_osnd[[i]][1]-x_osnd[[i]][2]))
  s_indx <- x_osnd[[i]][2]
  o_indx <- x_osnd[[i]][1]
  s_amp <- (pulse[,i][s_indx] - pulse[,i][o_indx])/2
  
  syst_y <- s_amp*sin(syst_b*(pulse[,1])+0)+(osnd[[i]][2]-s_amp)
  pre_phase_test[[i-1]] <- syst_y
  
  lines(pulse[,1], syst_y, col='red')

  
  #estimate phase from the sine wave
  #this DOESN'T WORK AND I DON'T KNOW WHY
  syst_y_peak <- findpeaks(syst_y) 
  xval <- approx(x=syst_y, y = pulse[,1], xout = syst_y_peak[2])$y
  
  #phase estimation - 360 * time delay / period
  phase_est <- (360*(x_osnd[[i]][2]-xval))/((2*pi)/syst_b)
  
  syst_y_phase <- s_amp*sin(syst_b*(pulse[,1])+phase_est)+(osnd[[i]][2]-s_amp)
  post_phase_test[[i-1]] <- syst_y_phase
  
  lines(pulse[,1], syst_y_phase, col = 'green')
  

}



########

#Plot average trace
averagetrace<-rowMeans(pulse[-1])
means <- data.frame(id=1:length(averagetrace), av=averagetrace)
ggplot(data = pulse_stacked, aes(x = pulse_stacked$x, y = pulse_stacked$values, col = pulse_stacked$wave_ID)) +
  geom_line(size = 1.5)+
  geom_line(data=means, aes(x=id, y=av), color="black")

	
	# Sense rapid increase as a possible beat                                         

       # Empirical checks of whether the current increase in the trace value looks
       # like a beat profile.  There should be a resolvable peak in the gradient,
	     # or else we should wait for more context.  The residual should have
	     # increased, so we don't double count activity from a previous beat.
       
       # If the same beat is already in the history, it will be processed elsewhere.
       # Attempt to match the beat profile to the empirical template.
  
  }


#Refitbeats
# Beats are added to record keeping once the first peak is detected, but some
# analysis is delayed until the third peak can be detected, or a subsequent
# beat prevents the tail of the beat being observed.


RefitBeats <- function()
 {
       # Try to fit fourth peak
       # For all unconfirmed beats in the history, attempt to refit with current trace data.
       # Do not attempt to fit later profiles; We might incorrectly accept a secondary peak.
       # Find valid beats and exclude false positives to an acceptable degree.  The beat must
       # be of sufficient amplitude and enough time must have passed since the previous beat, 
       # but there is some subtlety in combining these two requirements.
       # The fourth peak is added quite late, as the beat is generally accepted once the third 
       # peak is resolved.  Thus, the residual needs to be updated specifically.  The other 
       # peaks are handled during ConfirmBeat.
 }


#UpdatePPGanalysis

UpdatePPGanalysis <- function(). 
{
        #Enter calibration mode if the current interval between beats is too long.
        #Validate peak of last beat if sufficient time has passed
  
  }




#CalcMaxGradient

CalcMaxGradient <- function(i_amplitude, i_HasSkippedBeat, io_Interval)
{
  
  }

#TestMergeBeats
# Check whether the last beat in the history should be confirmed on the basis
# of timings, or a confirmed beat should be replaced with this new beat.

TestMergeBeats <- function()
 {
  
  }

#TestSkippedBeat
# Sometimes, in some individuals, a fraction of heart beats are 'skipped'.
# Physiologically this is probably because the heart tries to beat too early,
# so there's no pressure wave to propagate out to the PPG.  If we assume the
# existence of this skipped beat, the overall heart rate is nearly unaffected.
# We therefore need to detect skipped beats to keep meaningful statistical
# measures, which themselves are used to detect subsequent beats.

TestSkippedBeat <- function()
{
    # Multiple skipped beats may indicate that we've just got the wrong rate...
  }


#UndoBeat
# Reject a beat assignment;
# * Recalculate the residue after subtracting all beats from the trace
# * Recalculate the Bayesian model for the P50 statistic.

UndoBeat <- function()
{
  
  }

#UpdateAverageBeat 

UpdateAverageBeat <- function()
{
  
  }

#CalibrateBeat
#Method to find pulse 

CalibrateBeat <- function()
{
  
  }

#FindPeak

FindPeak <- funtion()
{
    # Look for a well-defined peak
    # If there's no definite maximum, look for a maximum in the first
	  # derivative.  Sometimes the 'secondary' peak rises from the 'primary'
	  # without the gradient decreasing to zero.
  
  }

#HasCalibrationData

HasCalibrationData <- function()
{
  
  }


TestArea <- function()
{
  
  }

# Accept that we have found a new beat if we can adequately fit it with the following basic model; 
# there exists a primary peak of a certain shape, whose echo will follow a certain time later. There is aslo
# a secondary component which (may) be added to this with its own echo.
#Once each component is located, substract it from the trace to make it easier to find the components
FitBeat <- function(FNTBeatParams, io_Beat)
{                                                   
  
  }

FindStartIndex <- function(i_PeakIndex)
{
  
  }

#inputs for this not sure about in R
FitPeaks <- function(i_argW, i_val0, i_start, i_end, TOptional<FVector2D> (o_Peaks)[FitBeatConst::NUM_PEAKS]) 
{
 sfunction<- splinefun(1:500, undetrended_data$PPG.PulseOx1[1:500], method="natural")
  deriv1<-sfunction(1:500, deriv = 1) #first derivative of the spline function
  plot(deriv1, type="l")
  p<-findpeaks(deriv1, nups = 6, minpeakdistance = 40) #finding maximum of first derivative
    plot(deriv1, type = 'l')
        points(p[,2], p[,1], col = 'red', pch = 19)
      peakindex<-p[,2]

  }
  ###This is just an application for fun to visualise the derivs:
  #splinevectorx<-as.numeric(unlist(spline[1]))
  #splinevectory<-as.numeric(unlist(spline[2]))
  #TkSpline(x=splinevectorx, y=splinevectory)
  
  }

SplineFitting <- function(io_Beat)
{
  
  }

TranscribeTrace <- function(o_SplineData, i_Start, i_Length, i_Offset)
{
  
  }

#Try to find events for the current pulse, based on Elgendi, Liang & Ward (2018)
# Events are intended to be unambiguous features in the PPG trace or its derivatives
CalculateEvents <- function(io_Beat)
{
  
  }

#PPG traces must be aligned and scaled if they are to be averaged.  Defining
# 'W' as time of steepest increase as the beat starts, 'U' and 'V' correspond
# to the time when the gradient is half as steep (the gradient may not go to
# zero in a consistent manner).  Assigning 'U' a value of zero, and 'V' a value
# of one (the y scale and offset may vary between beats), the time at which the
# scaled trace has a value of one half is used for time alignment.
FindReferenceTime <- function(Spline)
{
  
  }

# Look for 'O' and 'S' - the beginning and end of the main increase
FindEventsTrace <- function(io_Beat, Spline)
{
  
  }

#Main peak in the first derivative corresponds to the onset of the beat, powered
# by the contraction of the left ventricle
FindEvents1Dev <- function(io_Beat, Spline)
{
  
  }


 

#PPG traces must be aligned and scaled if they are to be averaged.  Defining
#'W' as time of steepest increase as the beat starts, 'U' and 'V' correspond
#to the time when the gradient is half as steep (the gradient may not go to
#zero in a consistent manner).  Assigning 'U' a value of zero, and 'V' a value
#of one (the y scale and offset may vary between beats), the time at which the
#scaled trace has a value of one half is used for time alignment.
FindReferenceTime  <- function(FNTSpline, Spline)
{
}

#Look for 'O' and 'S' - essentially the start and end of the main increase.
FindEventsTrace <- function(FNTBeatParams, io_Beat, const FNTSplin, Spline)
{
}

#Main peak in the first derivative corresponds to the onset of the beat - powered
#by the contraction of the left ventricle.
FindEvents1Dev <- function(FNTBeatParams, io_Beat, const FNTSpline, Spline) 
{
	}

#Look for a series of three maxima (and two minima) in the second derivative
#of the PPG trace.
#If we found too many, reject the smallest ones to try to get to what by eye
#are the `obvious' peaks
#If we found too few, loosen our definition.  Look for inflections, 
FindEvents2Dev<- function(FNTBeatParams, io_Beat, FNTSpline, Spline) 
{
	}


#Having labelled a maximum in the first derivative as 'Z', we expect to find
#the notch ('N') before it, and the diastolic peak ('D') after it.  Many
#people do not have the peak, or it is only sometimes apparent.  In these
#cases, we equate N, Z, and D.

#Typically the first derivative has two minima exist between its main peak
#('W') and Z; these we call 'X' and 'Y'.  Often there's only one real minimum.
#If this is skewed early in the interval, it's probably X, if it's skewed
#late, it's probably Y.  We might be able to label the other event based on
#less strong criteria.
FindEventsNotch<-function(FNTBeatParams, io_Beat, FNTSpline, Spline)
	{
	}

FindZ<-function(FNTSpline, i_Spline, i_argXY, float i_rangeMax)
	{
	}

#Generate a spline for the beat, with various protections
#Spline x-coordinate will be time in units of seconds,
#with zero value at 'W', the maximum first derivative.
#'Truncated' means that we'd use more of the PPG buffer if we had it
GetBeatSpline<-function(FNTBeatParams, i_Beat, FNTSpline, o_Spline) 
	{
	}

GetBeatSpline<-function(FNTBeatParams, i_Beat, FNTSpline, o_Spline, i_AllowTruncation, o_truncated) 
	{
	}

HasEventPeaks<-function(FNTBeatParams, i_Beat)
	{
	}

Stats<-function(i_SampleSize, PublicStats, o_Stats, PrivateStats, o_Private)
	{
	}

CalculateResidual<-function(i_Index)
	{
	}

SubtractBeat<-function(FNTBeatParams, i_Beat)
	{
	}

TestClamp<-function(FNTSpline, i_Spline, Start, MaxFirstDerivative)
	{
	}

SubtractSine<-function(io_Data, Count, Peak, Amplitude, Width, io_Fit)
{
	}

SinePeak<-function(DeltaTime, Amplitude, Width)
	{
	}
