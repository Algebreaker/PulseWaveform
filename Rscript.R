setwd("/home/johanna/Documents/Ninja theory/PulseAnalysis/Data/Craig")
data <- read.table("Source.csv", header=T, sep=",") #first line of the csv file needs to be deleted
library(tidyverse)                                  #Have tidyverse packages installed and call tidyverse in library()
library(TeachingDemos)
library(splines2)
library(pracma)

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


# DownSample 

list<-rle(data$PPG.PulseOx1)

ID <- rep(1:length(list$values), times = list$lengths)
data2 <- cbind(data, ID)       # Cbind a numbered column rather than true/false

data_downsampled <-c()              # Create empty vector

for (i in 1:max(ID)){
  sub.data <- filter(data2, ID == i) # Isolate rows with all the same value as a subset
  if(nrow(sub.data) <= 4){           # If the subset has 4 rows or less, add the first row to the empty vector
    data_downsampled <- rbind(data_downsampled, sub.data[1,])
  }else if(nrow(sub.data) > 4 ){data_downsampled <- rbind(data_downsampled, sub.data[1,], sub.data[5,])} #If the data has more than four rows, add the first and fifth
}
  # This may need to be adjusted if we have datasets where true values are repeated more than twice


# DownSample 
ID <- rep(1:length(list$values), times = list$lengths)     
data2 <- cbind(pulse_full, ID)                                  # Cbind a numbered column rather than true/false

data_downsampled <-c()                                          # Create empty vector

for (i in 1:max(ID)){                                           
  
  sub.data <- filter(data2, ID == i)                            # Isolate rows with all the same value as a subset
  
  if(nrow(sub.data) <= 4){                                      # If the subset has 4 rows or less, add the first row to the empty vector
    
    data_downsampled <- rbind(data_downsampled, sub.data[1,])   
    
  }else if(nrow(sub.data) > 4 ){data_downsampled <- rbind(data_downsampled, sub.data[1,], sub.data[5,])}  #If the data has more than four rows, add the first and fifth
  
}          # This may need to be adjusted if we have datasets where true values are repeated more than twice

#Undetrend
# Analysis of device output indicates that the PPG signal is detrended by application of the following
# formula: OUT[i] = 80 + (OUT[i-1]-80) * 0.96875 + (IN[i] - [IN[i-1]), where the constance 0.96875 is
# an approximation fitted to the data.
# Individual pulse events are more comprehensible if the detrending is not used, so this function 
#removes it by inverting the above function. 
   # can simply reverse the equatino as in the c++ script. 
undetrended <-c()    
for (i in 2:length(data_downsampled$PPG.PulseOx1)+1)   
{
  undetrended[i-1]<-((data_downsampled$PPG.PulseOx1[i]-80) - ((data_downsampled$PPG.PulseOx1[i-1]-80) * 0.96875) + (data_downsampled$PPG.PulseOx1[i-1]))
}
undetrended_data<-cbind(data_downsampled,undetrended)


## create spline + derivatives 
sfunction <- splinefun(1:1000, undetrended_data$undetrended[1:1000], method = "natural")
spline <- sfunction(seq(1, 1000, 0.1), deriv = 0)
deriv1 <- sfunction(seq(1, 1000, 0.1), deriv = 1)
deriv2 <- sfunction(seq(1, 1000, 0.1), deriv = 2)

# Finding peaks of the first derivative
pd1 <-findpeaks(deriv1, nups = 70, minpeakdistance = 10)     # this will need to be fine tuned

# Plot points of peaks on first derivative
plot(deriv1, type = "l")
points(pd1[,2], pd1[,1], col = 'red', pch = 19)

# Create vector of all x axis coordinates of peaks of first derivative
pd1index <- pd1[, 2]

# Plot deriv1 peaks on original spline
plot(spline, type = "l")
points(pd1index, spline[pd1index], col = 'red', pch = 19)


#Finding inflexion points

posneginflexionpoints <- c()

for(i in 1:(length(deriv1)-1)){
  if  (deriv1[i] > 0){
    if (deriv1[i+1] < 0){
      newvector <- c(deriv1[i], deriv1[i+1])
      if ((which(abs(0-newvector) == min(abs(0 - newvector)))) == 1){
        posneginflexionpoints <- c(posneginflexionpoints, i)
      }else {
        posneginflexionpoints <- c(posneginflexionpoints, i+1)
      }
    }
  }else if(deriv1[i] < 0){
    if(deriv1[i+1] > 0){
      newvector <- c(deriv1[i], deriv1[i+1])
      if ((which(abs(0-newvector) == min(abs(0 - newvector)))) == 1){
        posneginflexionpoints <- c(posneginflexionpoints, i)
      }else {
        posneginflexionpoints <- c(posneginflexionpoints, i+1)
      }
    }
  }
}

# Plot inflection points on first derivative
plot(deriv1, type = "l")
points(posneginflexionpoints, deriv1[posneginflexionpoints], col = 'red', pch = 19) # plots inflexion points on 1st deriv

# Plot inflection points on original spline
plot(spline[1:10000], type = "l")
points(posneginflexionpoints, spline[posneginflexionpoints], col = 'red', pch = 19)











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
