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
spline <- sfunction(seq(1, 1000, 0.1), deriv = 0)
deriv1 <- sfunction(seq(1, 1000, 0.1), deriv = 1)
deriv2 <- sfunction(seq(1, 1000, 0.1), deriv = 2)

# Turning the undetrended data into a piece-wise polynomial spline (non-discrete): 
spline_poly <-CubicInterpSplineAsPiecePoly(1:1000, undetrended_data$undetrended[1:1000], "natural")


###find peak values for individual thresholding
hdat<-hist(deriv1)
quantiles<-quantile(deriv1,probs=c(.025,.95))
threshold<-quantiles[2]
#as density plot
plot(density(deriv1))



# Finding peaks of the first derivative
pd1 <-findpeaks(deriv1, threshold=threshold)     # this will need to be fine tuned

# Plot points of peaks on first derivative
plot(deriv1, type = "l")
points(pd1[,2], pd1[,1], col = 'red', pch = 19)

# Finding second peaks of the first derivative (point "s") by first finding both w and s
pd2 <-findpeaks(deriv1, threshold=0.3)     # this will need to be fine tuned
# Plot points of peaks on first derivative
plot(deriv1, type = "l")
points(pd2[,2], pd2[,1], col = 'red', pch = 19)
#now find Z by cancelling out W
if (pd1[1]==pd2[1]){pd3 = pd2[seq(2, nrow(pd2), 2), ]# keep all even numbered elements in pd2 if the first value in the two vectors is the same (i.e. a w point)
}else { #else keep all odd numbered elements
  pd3 = pd2[seq(1, nrow(pd2), 2), ]
}

# Create vector of all x axis coordinates of peaks of first derivative
pd1index <- pd1[, 2]

# Plot deriv1 peaks on original spline
plot(spline, type = "l")
points(pd1index, spline[pd1index], col = 'red', pch = 19)



## Finding inflexion points on spline_poly
# Find all x values for inflexion points (points on deriv1 that equal 0):
inflexion_points <- solve(spline_poly, b = 0, deriv = 1)

# Find the y values for inflexion points:
inflexion_points_yval <- predict(spline_poly, inflexion_points)

# Plot the y values:
plot(spline_poly)
points(inflexion_points, inflexion_points_yval, pch = 19)




###lucie
subspline <- data.frame(spline[1:80])
subspline$x <- seq.int(nrow(subspline))
plot(subspline$x, subspline$spline.1.80., type = 'l')
points(posneginflexionpoints, subspline$spline.1.80.[posneginflexionpoints], col = 'red', pch = 19)
b <- 2*pi/(abs(4*(posneginflexionpoints[1]-posneginflexionpoints[2])))
indx <- posneginflexionpoints[2]
amp <- subspline$spline.1.80.[indx]
y <- amp*sin(b*(subspline$x)+10)
lines(subspline$x, y)



###### Chopping up and plotting all waveforms V2 ########

## Finding W
# First plot 1st deriv:
plot(spline_poly, deriv = 1)

# Finding inflexion points on deriv1 requires redefining 1st deriv as a piece-wise spline
deriv1_poly <- CubicInterpSplineAsPiecePoly(1:1000, deriv1, "natural") 

# Find inflexion points on deriv1_poly
inflexion_points_deriv1 <- solve(deriv1_poly, b = 0, deriv = 1)
inflexion_points_deriv1_yval <- predict(deriv1_poly, inflexion_points_deriv1)
plot(deriv1_poly)
points(inflexion_points_deriv1, inflexion_points_deriv1_yval, pch = 19)

# Identifying peaks of deriv1_poly:
w_poly_peaks <- predict(deriv1_poly, inflexion_points_deriv1)
w_poly_peaks <- which(w_poly_peaks > 0.5)     # this is still a cutoff specific to this data
w_poly_peaks <- inflexion_points_deriv1[w_poly_peaks]

# Plot back on spline_poly:
plot(spline_poly)
w_poly_peaks_yval <- predict(spline_poly, w_poly_peaks)
points(w_poly_peaks, w_poly_peaks_yval, pch = 19)


## Chopping up and plotting all waveforms:

## Find half the height of w (on derivative y-axis)

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

# Now need to scale according to v-u difference:

# Find individual u and v y-values 

u_yval <- u_v_yval[seq_along(u_v_yval) %%2 != 0] 
v_yval <- u_v_yval[seq_along(u_v_yval) %%2 == 0]   

# Find v-u differences

v_minus_u <- v_yval - u_yval


# You can't really scale spline_poly waves, but you can use the more accurate u, v, and w values
# to scale the original discrete spline correctly. 

## Building a stacked data_frame

pulse <- data.frame()

for(i in 1:length(w_poly_peaks)){
  
  spline_poly_wave_subset <- CubicInterpSplineAsPiecePoly((round((w_poly_peaks[i])-10)):(round((w_poly_peaks[i])+51)), spline[(round((w_poly_peaks[i])-10)):(round((w_poly_peaks[i])+51))], "natural")
  
  # Now get the y-values to fill the dataframe using the predict function
  
  xxxx <- predict(spline_poly_wave_subset, c(seq(round((w_poly_peaks[i])-10), 
                                                 round((w_poly_peaks[i])-1), 0.1), 
                                             w_poly_peaks[i],  
                                             seq(round((w_poly_peaks[i])+1), 
                                                 round((w_poly_peaks[i])+50), 0.1)))    
  
 xxxx <-  as.data.frame(xxxx)
 xxxx <- cbind(xxxx, c(seq(round((w_poly_peaks[i])-10), round((w_poly_peaks[i])-1), 0.1), w_poly_peaks[i],  seq(round((w_poly_peaks[i])+1), round((w_poly_peaks[i])+50), 0.1)))
 colnames(xxxx) <- c('y', 'x') 
 xxxx[, 2] <- xxxx[, 2] - (xxxx$x[1]-1)
 xxxx$wave <- i 
 # Need to scale so that v-u = 1
 xxxx$y <- xxxx$y/(v_minus_u[i])
 # Need to calculate difference between this w point and the first w point
 y_axis_difference <- w_poly_peaks_yval[1] - xxxx$y[92]
 xxxx$y <- xxxx$y + y_axis_difference
 # Need to adjust x values so that all w's line up on x-axis
 x_axis_difference <- w_poly_peaks[1] - xxxx$x[92]
 xxxx$x <- xxxx$x + x_axis_difference
 # Adjust such that w = 0.5, u ~ 0, v ~ 1
 xxxx$y <- xxxx$y - 80.25
 
 
 pulse <- rbind(pulse, xxxx)
 
}

# Plot it
pulse$wave <- as.factor(pulse$wave)
ggplot(data = pulse, aes(x = pulse$x, y = pulse$y, col = pulse$wave)) + geom_line(size = 1.5)

############







########

#Plot average trace
averagetrace<-rowMeans(pulse[-1])
means <- data.frame(id=1:length(averagetrace), av=averagetrace)
ggplot(data = pulse_stacked, aes(x = pulse_stacked$x, y = pulse_stacked$values, col = pulse_stacked$wave_ID)) +
  geom_line(size = 1.5)+
  geom_line(data=means, aes(x=id, y=av), color="black")







#FindNewBeat

FindNewBeat <- function()
{
       ## This is where splines are first fitted 
	
	splined_data <- spline(1:200, undetrended_data[1:200])    # extend subset as required, replace 1:200 with time data?
	plot(splined_data)                                        
	
       # Sense rapid increase as a possible beat                                         # 

	spline <- spline(1:500, undetrended_data$PPG.PulseOx1[1:500])  # (in this case to first 100 values)
	plot(spline, type='l')
	  #plot(deriv1, type='l')
	xcoord<-as.vector(spline$x)
	closestsplines <-c()    #gives you the index of the closes x-coordinate in the spline for each peak
	for (i in 1:length(peakindex))
	{
	  closestsplines[i]<-which(abs(xcoord-(peakindex[i]))==min(abs(xcoord-(peakindex[i]))))
	}
	ycoord<-spline$y
	peakycoord<-ycoord[closestsplines]
	points(spline$x[c(closestsplines)], spline$y[c(closestsplines)], col = 'red', pch = 19)
   
	
	#Fitting a sine curve   
	
	#create a subset of data 
	
subspline <- data.frame(spline[1:80]) 
subspline$x <- seq.int(nrow(subspline))

#plot the subset with the inflexion points 
plot(subspline$x, subspline$spline.1.80., type = 'l')
points(posneginflexionpoints, subspline$spline.1.80.[posneginflexionpoints], col = 'red', pch = 19)

#identify the different variables to be entered into the sine eq
b <- (2*pi)/(abs(2*(posneginflexionpoints[1]-posneginflexionpoints[2])))
s_indx <- posneginflexionpoints[2]
o_indx <- posneginflexionpoints[1]
amp <- (subspline$spline.1.80.[s_indx]-subspline$spline.1.80.[o_indx])/2

#specify sine and plot
y <- amp*sin(b*(subspline$x)+3.4444)+(84-amp)
lines(subspline$x, y)

	
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
