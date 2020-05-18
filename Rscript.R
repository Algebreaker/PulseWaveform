setwd("/home/johanna/Documents/Ninja theory/PulseAnalysis/Data/Craig")
data <- read.table("Source2.csv", header=T, sep=",") #first line of the csv file needs to be deleted

library(tidyverse)                                  #Have tidyverse packages installed and call tidyverse in library()
library(TeachingDemos)
library(splines2)
library(pracma)
library(SplinesUtils)

data<-data[!(data$PPG.PulseOx1=='NaN'),]

min_beats_per_minute <- 30
max_beats_per_minute <- 240
min_interval <- 0.25
max_interval <- 4
min_amplitude <- 0.5
min_frac_amplitude <- 0.2
sample_rate <- 75
peak_interval <- 9
peak_width <- 2.7   # Should peak width be an array of differing widths per peak?
skipped_beat_fraction <- 1.75

#main function
FNTBioRadioPulse<-function(input)
{
}


# Ignore public functions for now
# Ignore RestoreState for now


#Downsample
# The BioRadio device provides 250 samples per second, but the PPG is only sampled 75 times per second, 
# so we typically have repeated values in a pattern of 3-3-4 repeating. DownSample tries to retrieve the 
# unique values, with an effort to be robust against variation in the repeat pattern and also against 
# genuine repeated values.

list<-rle(data)
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
sfunction <- splinefun(1:1000, undetrended_data$undetrended[1:1000], method = "natural")
spline <- sfunction(seq(1, 1000), deriv = 0)
deriv1 <- sfunction(seq(1, 1000), deriv = 1)
deriv2 <- sfunction(seq(1, 1000), deriv = 2)

spl<-CubicInterpSplineAsPiecePoly(1:50, undetrended_data$undetrended[1:50], "natural")
yval <- solve(spl, b = 78) #returns a vector with all equations of piecewise polynomials that yield this y-value (2.85 in this case)

########## 

# Turning the undetrended data into a piece-wise polynomial spline (non-discrete): 
spline_poly <-CubicInterpSplineAsPiecePoly(1:1000, undetrended_data$undetrended[1:1000], "natural")

## Finding W

# Finding inflexion points on deriv1 requires redefining 1st deriv as a piece-meal spline
deriv1_poly <- CubicInterpSplineAsPiecePoly(1:1000, deriv1, "natural") 
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
# Plot back on spline_poly:
plot(spline_poly)
w_poly_peaks_yval <- predict(spline_poly, w_poly_peaks)
points(w_poly_peaks, w_poly_peaks_yval, pch = 19)

##Finding U and V

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


## Chopping up the original data_undetrended into individual waves:
sourcedata <- undetrended_data$undetrended[1:1000]
pulse <- data.frame(1:91)
for(i in 1:length(w_poly_peaks)){
  if(i == 1){
    pulse <- cbind(pulse, sourcedata[1:91])     # Special case for first wave needed as you cannot specify elements further back than 0
    }
  else 
    pulse <- cbind(pulse, sourcedata[(round(w_poly_peaks[i]) - 15):(round(w_poly_peaks[i]) + 75)])
    colnames(pulse)[i+1] <- paste("wave", i, sep = "_") 
}
colnames(pulse)[1] <- "x"
# Find the minimum of each chopped waveset within first 30 elements (should be approximately* o), 
# and then take the y values down by that much such that o is zero. 
precurser_o <- c()
for(i in 2:ncol(pulse)){
  wave <- pulse[, i]
  wave30 <- wave[1:30]
  precurser_o[i-1] <- min(wave30)
}
for(i in 2:ncol(pulse)){
  pulse[, i] <- pulse[, i] - precurser_o[i-1]
}


## Scale each wave by its own scalar
for(i in 2:ncol(pulse)){                                  
  pulse[, i] <- pulse[, i]/v_minus_u[i-1]        
}


## Now that data is scaled and y-axis normalized*, create a new polynomial spline for each wave
poly_wave <- list()
for(i in 2:ncol(pulse)){
  poly_wave[[i-1]] <-CubicInterpSplineAsPiecePoly(1:length(pulse[, i]), pulse[, i], "natural")
}


## Find U, V, W, O, S, N, D on the new polynomial splines:
osnd <- list()

for(i in 2:ncol(pulse)){
  
  sfunction2 <- splinefun(1:91, pulse[, i], method = "natural")
  deriv1_wave <- sfunction2(seq(1, 91), deriv = 1)
  deriv1_wave_poly <- CubicInterpSplineAsPiecePoly(1:91, deriv1_wave, "natural") 
  
  # Find inflexion points on deriv1_wave_poly
  inflexion_points_deriv1_wave_poly <- solve(deriv1_wave_poly, b = 0, deriv = 1)
  inflexion_points_deriv1_wave_poly_yval <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)
  # plot(deriv1_wave_poly)
  # points(inflexion_points_deriv1_wave_poly, inflexion_points_deriv1_wave_poly_yval, pch = 19)
  
  # Find correct threshold using histogram
  # hdat<-hist(deriv1_wave)
  quantiles<-quantile(deriv1_wave, probs=c(.025,.95))
  threshold<-quantiles[2]
  
  # Identifying peaks of deriv1_poly:
  w_poly_peaks_wave <- predict(deriv1_wave_poly, inflexion_points_deriv1_wave_poly)
  w_poly_peaks_wave <- which(w_poly_peaks_wave > threshold)     
  w_poly_peaks_wave <- inflexion_points_deriv1_wave_poly[w_poly_peaks_wave]
  
  # Plot on deriv1_wave_poly
  # plot(deriv1_wave_poly)
  # w_poly_peaks_yval <- predict(deriv1_wave_poly, w_poly_peaks_wave)
  # points(w_poly_peaks_wave, w_poly_peaks_yval, pch = 19)
	
  #Identifying 'notch' of deriv1_poly (for use on non-canonical waveforms):
  # take minimum inflection point and then find the next one
  notch <- inflexion_points_deriv1_wave_poly[(which(inflexion_points_deriv1_wave_poly_yval == min(inflexion_points_deriv1_wave_poly_yval)))+1]
  notch_yval <- inflexion_points_deriv1_wave_poly_yval[(which(inflexion_points_deriv1_wave_poly_yval == min(inflexion_points_deriv1_wave_poly_yval)))+1]
  plot(deriv1_wave_poly)
  points(notch, notch_yval, pch = 19)
  #Find corresponding y value on poly_wave
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
  #Find the four inflexion points that are 1 to the left and 3 to the right of W
  o <- max(which(inflexion_points_new < w_poly_peaks_wave[1]))
  inflexion_points_new_yval <- predict(poly_wave[[i-1]], inflexion_points_new)
  osnd[[i-1]] <- inflexion_points_new_yval[o:(o+3)]
  osnd[[i-1]] <- osnd[[i-1]] - inflexion_points_new_yval[o]
  
  # Plot back on poly_wave[[i]]:
  plot(poly_wave[[i-1]])
  w_poly_peaks_yval <- predict(poly_wave[[i-1]], w_poly_peaks_wave)
  points(w_poly_peaks_wave[1], w_poly_peaks_yval[1], pch = 19)
  points(inflexion_points_new, inflexion_points_new_yval, pch = 19)
  points(half_heights_wave_new, u_v_yval_wave, pch = 19)
  
}





########

#Plot average trace
averagetrace<-rowMeans(pulse[-1])
means <- data.frame(id=1:length(averagetrace), av=averagetrace)
ggplot(data = pulse_stacked, aes(x = pulse_stacked$x, y = pulse_stacked$values, col = pulse_stacked$wave_ID)) +
  geom_line(size = 1.5)+
  geom_line(data=means, aes(x=id, y=av), color="black")

	
#estimating S values for data without clear systolic peaks
# S value is estimated as W + 2*(V-W)
s_index <- c()

for(i in 1:length(pd1index)){
s_index[i] <- pd1index[i] + (2*(v_index[i]-pd1index[i]))
}
	
#estimating O values 
#O is estimated as U - (W-U)
o_index <- c()

for(i in 1:length(u_index)){
o_index[i] <- u_index[i] - (pd1index[i]-u_index[i])
}

#plotting estimted S and O values on original spline

plot(spline, type = "l")
points(u_index, spline[u_index], col = "red", pch = 19)
points(v_index, spline[v_index], col = "red", pch = 19)
points(pd1index, spline[pd1index], col = 'blue', pch = 19)
points(s_index, spline[s_index], col="green", pch=19)
points(o_index, spline[o_index], col="green", pch=19)

#Fitting a sine curve badly
	
#create a subset of data 
	
subspline <- data.frame(spline[1:551]) 
subspline$x <- seq.int(nrow(subspline))

#plot the subset with the inflexion points 
plot(subspline$x, subspline$spline.1.551., type = 'l')
points(posneginflexionpoints, subspline$spline.1.551.[posneginflexionpoints], col = 'red', pch = 19)

#identify the different systolic variables to be entered into the sine eq
syst_b <- (2*pi)/(abs(2*(posneginflexionpoints[1]-posneginflexionpoints[2])))
s_indx <- posneginflexionpoints[2]
o_indx <- posneginflexionpoints[1]
s_amp <- (subspline$spline.1.551.[s_indx]-subspline$spline.1.551.[o_indx])/2

#specifysystolic peak
syst_y <- s_amp*sin(syst_b*(subspline$x)+3.744)+(84-amp)

#identify the different disatolic variables 
d_indx <- posneginflexionpoints[4]
d_amp <- (subspline$spline.1.551.[d_indx]-subspline$spline.1.551.[o_indx])/2

#specify diastolic peak
d_y <- d_amp*sin(syst_b*(subspline$x)-7.144)+(82-d_amp)

# Really awkardly take the data that we want to use and nothing else
syst_y[300:551]<-NA
d_y[0:200]<-NA
d_y[480:551] <- NA

#plot resulting sines on waveform
lines(subspline$x, syst_y, col='red')
lines(subspline$x, d_y, col='blue')

	
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
