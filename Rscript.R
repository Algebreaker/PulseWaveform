setwd("/home/johanna/Documents/Ninja theory/PulseAnalysis/Data/Craig")
data <- read.table("Source2.csv", header=T, sep=",") #first line of the csv file needs to be deleted

source("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/preproc.R")
source("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/find_w.R")
source("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/find_u_v.R")
source("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/baseline.R")
source("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/find_osnd.R")
source("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/spectrum.R")
source("/Users/luciedaniel-watanabe/Desktop/attempt at pulse analysis/find_wuv.R")



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

#Print all osnd's
for(i in 1:length(osnd_y)){
  osnd_correction <- osnd_y[[i]]
  osnd_correction <- osnd_correction - osnd_correction[1]
  osnd_y[[i]] <- osnd_correction
  print(osnd_y[[i]])
}



###PLEASE DON'T RUN ANY OF THE SINE STUFF FOR NOW, IT WILL NOT WORK!!###
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


####LUCIE'S PRELIMINARY REFITTING SPLINES BIT THAT WORKS?


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

for(i in 2:(ncol(pulse)-1)){

  fitted_s <- FALSE
  fitted_n <- FALSE
  fitted_d <- FALSE

  while(!fitted_s | !fitted_n | !fitted_d){
    
    ##define the peaks of the current waveform
    s_peak <- max(s_sine[[i-1]])
    n_peak <- max(n_sine[[i-1]])
    d_peak <- max(d_sine[[i-1]])
    d_peak_x <- match(d_peak, d_sine[[i-1]])
    n_peak_x <- match(n_peak, n_sine[[i-1]])
    s_peak_x <- match(s_peak, s_sine[[i-1]])
    
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
      
      d_n_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(pulse[, i]), d_n_resid, "natural")
      inflex_test_x <- solve(d_n_resid_spline, b=0, deriv=1)
      inflex_test_y <- predict(d_n_resid_spline, inflex_test_x)
      
      pk_loc <- which(abs(inflex_test_x-s_peak_x)==min(abs(inflex_test_x - s_peak_x)))
      
      s_x[i-1] <- inflex_test_x[pk_loc]
      s_y[i-1] <- inflex_test_y[pk_loc]

      
      fitted_s <- TRUE

    
    } else if(n_peak > d_peak){
    
      s_resid <- pulse[,i] - s_sine[[i-1]]
      s_d_resid <- s_resid - d_sine[[i-1]]
    
      s_d_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(pulse[, i]), s_d_resid, "natural")
      current_deriv_n <- predict(s_d_resid_spline, n_peak_x, deriv=1)
      
      inflex_test_x <- solve(s_d_resid_spline, b=0, deriv=1)
      inflex_test_y <- predict(s_d_resid_spline, inflex_test_x)
        
        if(current_deriv_n < 0){
          subset_x <- inflex_test_x[inflex_test_x > s_peak_x & inflex_test_x < n_peak_x]
          x_loc <- match(subset_x, inflex_test_x)
          subset_y <- inflex_test_y[x_loc]
          pk_loc <- which(abs(subset_x - n_peak_x)==min(abs(subset_x - n_peak_x)))

        } else if(current_deriv_n > 0){
          subset_x <- inflex_test_x[inflex_test_x > n_peak_x & inflex_test_x < d_peak_x]
          x_loc <- match(subset_x, inflex_test_x)
          subset_y <- inflex_test_y[x_loc]
          pk_loc <- which(abs(subset_x - n_peak_x)==min(abs(subset_x - n_peak_x)))

        } 
        
        n_y[i-1] <- subset_y[pk_loc]
        n_x[i-1] <- subset_x[pk_loc]
        
        if(fitted_s & (!fitted_d)){
        fitted_s = FALSE;
        }

      fitted_n = TRUE;
    
    
    } else{
    
      s_resid <- pulse[,i] - s_sine[[i-1]]
      s_n_resid <- s_resid - n_sine[[i-1]]
      
      s_n_resid_spline <-CubicInterpSplineAsPiecePoly(1:length(pulse[, i]), s_n_resid, "natural")
      current_deriv_d <- predict(s_n_resid_spline, d_peak_x, deriv=1)
      
      inflex_test_x <- solve(s_n_resid_spline, b=0, deriv=1)
      inflex_test_y <- predict(s_n_resid_spline, inflex_test_x)

      if(current_deriv_d < 0){
        subset_x <- inflex_test_x[inflex_test_x > n_peak_x & inflex_test_x < d_peak_x]
        x_loc <- match(subset_x, inflex_test_x)
        subset_y <- inflex_test_y[x_loc]
        pk_loc <- which(abs(subset_x - d_peak_x)==min(abs(subset_x - d_peak_x)))
      
      } else if(current_deriv_d > 0){
        subset_x <- inflex_test_x[inflex_test_x > d_peak]
        x_loc <- match(subset_x, inflex_test_x)
        subset_y <- inflex_test_y[x_loc]
        pk_loc <- which(abs(subset_x - d_peak_x)==min(abs(subset_x - d_peak_x)))

      }
      
      d_x_est <- subset_x[pk_loc]
      d_y_est <- subset_y[pk_loc]
        
     # if(d_x_est > (d_peak_x + max_shift)){
        #d_y_est <- inflex_test_y[which.max(inflex_test_y)]
        #d_x_est <- inflex_test_x[which.max(inflex_test_y)]
        #d_y_est = 2 * d_y_est * sample_rate * avg_period
        #d_x_est = d_x_est + 0.5 * pi * avg_period
    #  }
      
      d_y[i-1] <- d_y_est
      d_x[i-1] <- d_x_est
      
        
      fitted_d = TRUE
    
    }
  }
}
  
for(i in 2:(ncol(pulse)-1)){
  plot(pulse$x, pulse[,i], type='l')
  points(s_x[i-1], s_y[i-1], pch = 19, col ='red')
  points(n_x[i-1], n_y[i-1], pch = 19, col = 'blue')
  points(d_x[i-1], d_y[i-1], pch = 19, col = 'yellow')
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

